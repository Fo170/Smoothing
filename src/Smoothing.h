/*
  Smoothing.h - Bibliothèque de lissage / filtrage pour signaux de capteurs.
  ================================================================================
  Version : 1.1.1

  Header-only, SANS dépendance (ni Arduino ni STL) -> réutilisable sur tout
  projet C++ (ESP8266, ESP32, Arduino, PC).

  RÔLE
  ----
  La classe `Smoother` conserve les N dernières valeurs BRUTES d'un capteur
  (N = fenêtre, <= MAX) et renvoie une valeur « lissée » selon la méthode
  choisie. Elle est prévue pour être instanciée une fois PAR VOIE de mesure.

  UTILISATION TYPIQUE
  -------------------
      Smoother sm;                       // un par capteur/voie

      // ... à chaque acquisition (ex. toutes les 2 s) :
      sm.push(mesureBrute);              // mémorise + met à jour Kalman/cumul

      // ... quand on a besoin de la valeur lissée :
      float y = sm.value(Smoother::SMA,   // méthode
                         8,               // taille de fenêtre (nb d'échantillons)
                         3.0f,            // k (Hampel / z-score / deadband / a-b / Holt)
                         0.15f);          // proportion de trim (tronquée/winsorisée)

  Ordre recommandé : lisser le BRUT puis appliquer la calibration. Ne jamais
  moyenner des valeurs déjà calibrées (un changement de calibration fausserait
  la moyenne).

  LES 5 FAMILLES DE MÉTHODES (enum Type)
  --------------------------------------
  1) MOYENNES MOBILES
       SMA         : moyenne mobile simple
       WMA         : moyenne mobile pondérée (poids croissants)
       EMA         : moyenne exponentielle (poids décroissant avec l'âge)
       DEMA        : double EMA (retard réduit)
       TEMA        : triple EMA (retard minimal)
       MEDIAN      : médiane
       TRIMMED     : moyenne tronquée (retire les extrêmes)
       WINSORIZED  : moyenne winsorisée (borne les extrêmes)
       CA          : moyenne cumulée (tout l'historique)
  2) FILTRES ROBUSTES AUX PARASITES
       HAMPEL      : médiane + MAD, remplace les aberrants
       ZSCORE      : ne garde que les points dans k*sigma
       KALMAN      : filtre de Kalman 1D (optimal bruit + modèle)
       ALPHABETA   : filtre alpha-beta (position + vitesse)
  3) FILTRES DE SIGNAL (DSP)
       GAUSSIAN    : noyau gaussien (lissage doux)
       SAVGOL      : Savitzky-Golay ordre 2 (préserve les pics)
       IIR         : filtre passe-bas du 2e ordre (type Butterworth)
  4) LISSAGE ORIENTÉ TENDANCE
       HOLT        : double lissage exponentiel (niveau + tendance)
  5) ANTI-REBOND (hystérésis)
       DEADBAND    : ne change la sortie que si la variation dépasse un seuil

  Paramètres communs :
    - win   : taille de fenêtre (nb d'échantillons, 1..MAX)
    - k     : seuil générique (Hampel : nb de MAD ; z-score : nb d'écarts-types ;
              deadband : seuil absolu ; alpha-beta/Holt : gain)
    - trim  : proportion d'extrêmes (tronquée/winsorisée), 0..0.4

  NOTE ÉNUMÉRATION : NONE/SMA/EMA/MEDIAN conservent les valeurs 0..3 pour la
  compatibilité des configurations enregistrées ; les nouvelles méthodes sont
  ajoutées à la suite.
*/

#ifndef SMOOTHING_H
#define SMOOTHING_H

#include <math.h>
#include <stdint.h>

class Smoother {
  public:
    // ----------------------------------------------------------------------
    // Méthodes de lissage (voir le descriptif en tête de fichier).
    // ----------------------------------------------------------------------
    enum Type {
      NONE        = 0,   // pas de lissage
      SMA         = 1,   // moyenne mobile simple
      EMA         = 2,   // moyenne mobile exponentielle
      MEDIAN      = 3,   // médiane
      WMA         = 4,   // moyenne mobile pondérée
      DEMA        = 5,   // double EMA
      TEMA        = 6,   // triple EMA
      HAMPEL      = 7,   // filtre de Hampel (médiane + MAD)
      TRIMMED     = 8,   // moyenne tronquée
      WINSORIZED  = 9,   // moyenne winsorisée
      ZSCORE      = 10,  // rejet des points hors k*sigma
      KALMAN      = 11,  // filtre de Kalman 1D
      CA          = 12,  // moyenne cumulée
      // --- familles complémentaires ---
      GAUSSIAN    = 13,  // noyau gaussien (DSP)
      SAVGOL      = 14,  // Savitzky-Golay ordre 2 (DSP, préserve les pics)
      IIR         = 15,  // passe-bas 2e ordre type Butterworth (DSP)
      HOLT        = 16,  // double lissage exponentiel (tendance)
      ALPHABETA   = 17,  // filtre alpha-beta (robuste, position+vitesse)
      DEADBAND    = 18,  // hystérésis / anti-rebond
      TYPE_COUNT  = 19   // nb de méthodes (utile pour l'interface)
    };

    static const uint8_t MAX = 16;

    Smoother() { reset(); }

    // ----------------------------------------------------------------------
    // reset() : vide le buffer et remet à zéro l'état (Kalman, cumul, deadband).
    // ----------------------------------------------------------------------
    void reset() {
      _count = 0;
      _head  = 0;
      _kfX   = NAN;
      _kfP   = 0.0f;
      _kfQ   = 1e-4f;
      _kfR   = 1e-2f;
      _caSum = 0.0;
      _caN   = 0;
      _dbOut = NAN;
    }

    // ----------------------------------------------------------------------
    // push(v) : ajoute une nouvelle mesure BRUTE et met à jour l'état récursif
    // (Kalman, moyenne cumulée).
    // ----------------------------------------------------------------------
    void push(float v) {
      _buf[_head] = v;
      _head = (uint8_t)((_head + 1) % MAX);
      if (_count < MAX) _count++;

      // Filtre de Kalman 1D : prédiction + correction.
      if (isnan(_kfX)) {
        _kfX = v; _kfP = 0.0f;
      } else {
        _kfP += _kfQ;
        float denom = _kfP + _kfR;
        float K = (denom > 0.0f) ? (_kfP / denom) : 0.0f;
        _kfX += K * (v - _kfX);
        _kfP *= (1.0f - K);
      }

      // Moyenne cumulée.
      if (!isnan(v)) { _caSum += v; _caN++; }
    }

    float latest() const {
      return _count ? _buf[(_head - 1 + MAX) % MAX] : NAN;
    }

    void setKalman(float q, float r) { _kfQ = q; _kfR = r; }

    // ----------------------------------------------------------------------
    // value(type, win, k, trim) : renvoie la valeur lissée (const sauf DEADBAND
    // qui maintient un état d'hystérésis via un membre `mutable`).
    // ----------------------------------------------------------------------
    float value(uint8_t type, uint8_t win, float k = 3.0f, float trim = 0.15f) const {
      if (_count == 0) return NAN;
      if (win < 1)      win = 1;
      if (win > _count) win = _count;

      // --- Méthodes à état persistant (sans fenêtre) ---
      if (type == NONE)   return latest();
      if (type == KALMAN) return _kfX;
      if (type == CA)     return _caN ? (float)(_caSum / (double)_caN) : NAN;

      // DEADBAND : hystérésis — la sortie ne bouge que si |x - sortie| > k.
      // x = dernière valeur BRUTE (latest) : pour anti-rebonder une valeur lissée,
      // pousser la valeur lissée dans une 2e instance de Smoother.
      // `mutable` car un filtre à hystérésis est intrinsèquement à état.
      if (type == DEADBAND) {
        float x = latest();
        if (isnan(_dbOut)) _dbOut = x;
        else if (fabsf(x - _dbOut) > k) _dbOut = x;
        return _dbOut;
      }

      // Copie chronologique des `win` dernières valeurs.
      float v[MAX];
      extract(v, win);

      switch (type) {

        // ---- Famille 1 : moyennes mobiles ----
        case SMA: {                        // moyenne simple
          float s = 0.0f;
          for (uint8_t i = 0; i < win; i++) s += v[i];
          return s / (float)win;
        }
        case WMA: {                        // pondérée (récent = plus fort)
          float s = 0.0f, wsum = 0.0f;
          for (uint8_t i = 0; i < win; i++) {
            float wgt = (float)(i + 1);
            s += v[i] * wgt; wsum += wgt;
          }
          return wsum > 0.0f ? s / wsum : NAN;
        }
        case EMA: {                        // exponentielle
          float e[MAX];
          return emaStage(v, win, e);
        }
        case DEMA: {                       // double EMA
          float e1[MAX], e2[MAX];
          float m1 = emaStage(v, win, e1);
          float m2 = emaStage(e1, win, e2);
          return 2.0f * m1 - m2;
        }
        case TEMA: {                       // triple EMA
          float e1[MAX], e2[MAX], e3[MAX];
          float m1 = emaStage(v, win, e1);
          float m2 = emaStage(e1, win, e2);
          float m3 = emaStage(e2, win, e3);
          return 3.0f * m1 - 3.0f * m2 + m3;
        }
        case MEDIAN:
          sortAsc(v, win);
          return medianSorted(v, win);
        case TRIMMED: {
          sortAsc(v, win);
          uint8_t t = trimCount(trim, win);
          float s = 0.0f; uint8_t n = 0;
          for (uint8_t i = t; i < win - t; i++) { s += v[i]; n++; }
          return s / (float)n;
        }
        case WINSORIZED: {
          sortAsc(v, win);
          uint8_t t = trimCount(trim, win);
          float lo = v[t], hi = v[win - 1 - t];
          float s = 0.0f;
          for (uint8_t i = 0; i < win; i++) {
            float x = v[i];
            if (x < lo) x = lo;
            if (x > hi) x = hi;
            s += x;
          }
          return s / (float)win;
        }

        // ---- Famille 2 : filtres robustes ----
        case HAMPEL: {
          float last = v[win - 1];
          sortAsc(v, win);
          float med = medianSorted(v, win);
          float dev[MAX];
          for (uint8_t i = 0; i < win; i++) dev[i] = fabsf(v[i] - med);
          sortAsc(dev, win);
          float mad = medianSorted(dev, win);
          if (mad <= 1e-9f) return med;   // fenêtre quasi constante : jamais de pic
          if (fabsf(last - med) > k * mad) return med;
          return last;
        }
        case ZSCORE: {
          float mean = 0.0f;
          for (uint8_t i = 0; i < win; i++) mean += v[i];
          mean /= (float)win;
          float var = 0.0f;
          for (uint8_t i = 0; i < win; i++) { float d = v[i] - mean; var += d * d; }
          float sd = sqrtf(var / (float)win);
          if (sd <= 1e-9f) return mean;
          float s = 0.0f; uint8_t n = 0;
          for (uint8_t i = 0; i < win; i++)
            if (fabsf(v[i] - mean) <= k * sd) { s += v[i]; n++; }
          return n ? s / (float)n : mean;
        }
        case ALPHABETA: {                  // alpha-beta (position + vitesse)
          float alpha = 2.0f / ((float)win + 1.0f);
          float beta  = clampf(k / 10.0f, 0.02f, 0.8f);
          float x = v[0], vel = 0.0f;
          for (uint8_t i = 1; i < win; i++) {
            float pred = x + vel;
            float res  = v[i] - pred;
            x   = pred + alpha * res;
            vel = vel  + beta  * res;
          }
          return x;
        }

        // ---- Famille 3 : filtres de signal (DSP) ----
        case GAUSSIAN: {                   // noyau gaussien centré sur le récent
          float sigma = (float)win / 3.0f;
          float s = 0.0f, wsum = 0.0f;
          for (uint8_t i = 0; i < win; i++) {
            float d = (float)(win - 1 - i);   // âge (0 = plus récent)
            float wgt = expf(-0.5f * (d / sigma) * (d / sigma));
            s += v[i] * wgt; wsum += wgt;
          }
          return s / wsum;
        }
        case SAVGOL: {                     // Savitzky-Golay ordre 2
          // Ajuste une parabole sur la fenêtre et l'évalue au point le plus
          // récent : lisse en préservant les courbures/pics.
          if (win < 3) return latest();
          float mn = v[0], mx = v[0];
          double S0 = win, S1 = 0, S2 = 0, S3 = 0, S4 = 0, T0 = 0, T1 = 0, T2 = 0;
          for (uint8_t i = 0; i < win; i++) {
            if (v[i] < mn) mn = v[i];
            if (v[i] > mx) mx = v[i];
            double x = i, x2 = x * x;
            S1 += x; S2 += x2; S3 += x2 * x; S4 += x2 * x2;
            T0 += v[i]; T1 += x * v[i]; T2 += x2 * v[i];
          }
          double m[3][4] = {{S0, S1, S2, T0}, {S1, S2, S3, T1}, {S2, S3, S4, T2}};
          for (int col = 0; col < 3; col++) {           // élimination de Gauss
            int piv = col;
            for (int r = col + 1; r < 3; r++)
              if (fabs(m[r][col]) > fabs(m[piv][col])) piv = r;
            if (fabs(m[piv][col]) < 1e-12) return latest();
            if (piv != col)
              for (int kk = 0; kk < 4; kk++) { double t = m[col][kk]; m[col][kk] = m[piv][kk]; m[piv][kk] = t; }
            for (int r = col + 1; r < 3; r++) {
              double f = m[r][col] / m[col][col];
              for (int kk = col; kk < 4; kk++) m[r][kk] -= f * m[col][kk];
            }
          }
          double c2 = m[2][3] / m[2][2];
          double c1 = (m[1][3] - m[1][2] * c2) / m[1][1];
          double c0 = (m[0][3] - m[0][1] * c1 - m[0][2] * c2) / m[0][0];
          double xe = (double)(win - 1);
          float y = (float)(c0 + c1 * xe + c2 * xe * xe);
          // Bornage à la plage de la fenêtre : évite les « overshoots » absurdes
          // d'un ajustement parabolique en présence d'un pic aberrant.
          if (y < mn) y = mn;
          if (y > mx) y = mx;
          return y;
        }
        case IIR: {                        // passe-bas 2e ordre (Butterworth)
          float fc = 1.0f / (2.0f * (float)win);   // fréquence de coupure
          if (fc <= 0.0f || fc >= 0.5f) return latest();
          float w0 = 2.0f * 3.14159265f * fc;
          float cw = cosf(w0), sw = sinf(w0);
          float q  = 0.7071f;
          float al = sw / (2.0f * q);
          float b0 = (1.0f - cw) * 0.5f, b1 = 1.0f - cw, b2 = b0;
          float a0 = 1.0f + al, a1 = -2.0f * cw, a2 = 1.0f - al;
          // Initialisation de l'état sur la première valeur : réduit le
          // transitoire de démarrage (sans le supprimer ; sinon sortie biaisée
          // sur fenêtre courte).
          float x1 = v[0], x2 = v[0], y1 = v[0], y2 = v[0], y = v[0];
          for (uint8_t i = 0; i < win; i++) {
            y = (b0 / a0) * v[i] + (b1 / a0) * x1 + (b2 / a0) * x2
                - (a1 / a0) * y1 - (a2 / a0) * y2;
            x2 = x1; x1 = v[i]; y2 = y1; y1 = y;
          }
          return y;
        }

        // ---- Famille 4 : tendance ----
        case HOLT: {                       // double lissage exponentiel
          float alpha = 2.0f / ((float)win + 1.0f);
          float beta  = clampf(k / 20.0f, 0.01f, 0.5f);
          float level = v[0];
          float trend = (win > 1) ? (v[1] - v[0]) : 0.0f;
          for (uint8_t i = 1; i < win; i++) {
            float prev = level;
            level = alpha * v[i] + (1.0f - alpha) * (level + trend);
            trend = beta * (level - prev) + (1.0f - beta) * trend;
          }
          return level;                    // niveau lissé (tendance suivie)
        }
      }
      return latest();
    }

    // Nom lisible d'une méthode (interface web / logs).
    static const char* typeName(uint8_t t) {
      switch (t) {
        case NONE:       return "Aucun";
        case SMA:        return "Moyenne mobile (SMA)";
        case EMA:        return "Moyenne exponentielle (EMA)";
        case MEDIAN:     return "Médiane";
        case WMA:        return "Moyenne pondérée (WMA)";
        case DEMA:       return "Double EMA (DEMA)";
        case TEMA:       return "Triple EMA (TEMA)";
        case HAMPEL:     return "Hampel (médiane + MAD)";
        case TRIMMED:    return "Moyenne tronquée";
        case WINSORIZED: return "Moyenne winsorisée";
        case ZSCORE:     return "Rejet z-score";
        case KALMAN:     return "Filtre de Kalman 1D";
        case CA:         return "Moyenne cumulée";
        case GAUSSIAN:   return "Filtre gaussien";
        case SAVGOL:     return "Savitzky-Golay (ordre 2)";
        case IIR:        return "Passe-bas IIR (Butterworth)";
        case HOLT:       return "Holt (niveau + tendance)";
        case ALPHABETA:  return "Filtre alpha-beta";
        case DEADBAND:   return "Hystérésis (deadband)";
      }
      return "?";
    }

  private:
    float    _buf[MAX];
    uint8_t  _count;
    uint8_t  _head;
    float    _kfX, _kfP, _kfQ, _kfR;   // état Kalman
    double   _caSum;
    uint32_t _caN;                     // état moyenne cumulée
    mutable float _dbOut;              // état deadband (hystérésis)

    void extract(float* out, uint8_t win) const {
      for (uint8_t i = 0; i < win; i++) {
        int idx = (_head - win + i + MAX * 2) % MAX;
        out[i] = _buf[idx];
      }
    }

    static void sortAsc(float* a, uint8_t n) {
      for (uint8_t i = 0; i < n; i++)
        for (uint8_t j = i + 1; j < n; j++)
          if (a[j] < a[i]) { float t = a[i]; a[i] = a[j]; a[j] = t; }
    }

    static float medianSorted(const float* a, uint8_t n) {
      if (n == 0) return NAN;
      if (n & 1)  return a[n / 2];
      return (a[n / 2 - 1] + a[n / 2]) * 0.5f;
    }

    static float emaStage(const float* src, uint8_t n, float* out) {
      out[0] = src[0];
      float alpha = 2.0f / ((float)n + 1.0f);
      for (uint8_t i = 1; i < n; i++)
        out[i] = src[i] * alpha + out[i - 1] * (1.0f - alpha);
      return out[n - 1];
    }

    static uint8_t trimCount(float trim, uint8_t win) {
      if (trim < 0.0f) trim = 0.0f;
      if (trim > 0.4f) trim = 0.4f;
      uint8_t t = (uint8_t)(trim * (float)win + 0.5f);
      if (t * 2 >= win) t = (uint8_t)((win - 1) / 2);
      return t;
    }

    static float clampf(float v, float lo, float hi) {
      if (v < lo) return lo;
      if (v > hi) return hi;
      return v;
    }
};

#endif // SMOOTHING_H
