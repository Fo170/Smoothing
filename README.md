# Smoothing — lissage et filtrage de signaux de capteurs

Bibliothèque **header-only**, **sans dépendance** (ni Arduino ni STL), pour lisser
ou filtrer les mesures **brutes** d'un capteur. Réutilisable sur ESP8266, ESP32,
Arduino (AVR, SAM, STM32…) et PC.

Auteur : Olivier FOURNET (olivier.fournet@free.fr) — org GitHub Fo170  •  Licence : GPL-3.0  •  Version : 1.1.1

---

## Sommaire

- [Pourquoi / quand l'utiliser](#pourquoi--quand-lutiliser)
- [Installation](#installation)
- [Utilisation rapide](#utilisation-rapide)
- [API détaillée](#api-détaillée)
- [Les 5 familles de méthodes](#les-5-familles-de-méthodes)
  - [1. Moyennes mobiles](#1-moyennes-mobiles)
  - [2. Filtres robustes](#2-filtres-robustes)
  - [3. Filtres de signal (DSP)](#3-filtres-de-signal-dsp)
  - [4. Lissage orienté tendance](#4-lissage-orienté-tendance)
  - [5. Anti-rebond (hystérésis)](#5-anti-rebond-hystérésis)
- [Paramètres : comment les régler](#paramètres--comment-les-régler)
- [Guide de choix par situation](#guide-de-choix-par-situation)
- [Formules complètes](#formules-complètes)
- [Mémoire, CPU et performance](#mémoire-cpu-et-performance)
- [Intégration dans un projet (brut → lissage → calibration)](#intégration-dans-un-projet-brut--lissage--calibration)
- [Limites et pièges](#limites-et-pièges)
- [Exemples](#exemples)

---

## Pourquoi / quand l'utiliser

Un capteur réel est **bruité** : bruit blanc haute fréquence, pics ponctuels
(parasites électriques), dérive lente, valeurs aberrantes. Lisser permet
d'obtenir une mesure **stable** sans relire le capteur plus lentement.

La classe `Smoother` :

- conserve les **N dernières valeurs brutes** dans un **buffer circulaire** de
  taille fixe (mémoire bornée, pas d'allocation dynamique) ;
- renvoie une valeur lissée selon la **méthode** choisie ;
- peut fonctionner en **temps réel** (un appel à `push()` par acquisition, un
  appel à `value()` quand on a besoin de la valeur).

> ⚠️ **Règle d'or** : lisser le **BRUT** puis appliquer la **calibration**
> ensuite. Ne jamais moyenner des valeurs déjà calibrées : un changement de
> calibration fausserait la moyenne (mélange d'anciens et de nouveaux
> coefficients).

---

## Installation

### PlatformIO (recommandé — depuis GitHub)

La bibliothèque est hébergée sur GitHub. Déclarez-la dans le `platformio.ini` de votre projet :

```ini
[env:esp32dev]
platform  = espressif32
board     = esp32dev
framework = arduino
lib_deps  =
    https://github.com/Fo170/Smoothing.git
```

> ⚠️ Utilisez toujours l'**URL git complète** (`https://github.com/Fo170/Smoothing.git`),
> jamais la forme courte `Fo170/Smoothing` (non publiée sur le registre PlatformIO).

Pour **figer une version**, ajoutez un tag (ou une branche) :

```ini
lib_deps =
    https://github.com/Fo170/Smoothing.git#v1.1.1
```

PlatformIO télécharge alors la librairie automatiquement : elle n'est **pas** à copier dans le
dossier `lib/` du projet. Dans le code :

```cpp
#include <Smoothing.h>
```

### Arduino IDE

Copiez `Smoothing.h` dans le dossier de votre sketch ou dans
`libraries/Smoothing/src/`.

### Installation locale (développement / dépannage)

Pour modifier la librairie localement, placez le dossier `Smoothing/` dans `lib/` du projet :

```
lib/Smoothing/    (library.json + src/Smoothing.h)
```

---

## Utilisation rapide

```cpp
#include <Smoothing.h>

Smoother sm;                          // un Smoother PAR voie de mesure

void setup() {
  sm.setKalman(1e-4f, 1e-2f);         // optionnel (si méthode KALMAN)
}

void loop() {
  float brut = lireCapteur();

  sm.push(brut);                      // mémorise + met à jour l'état récursif
  float lisse = sm.value(Smoother::EMA,  // méthode
                         8,              // fenêtre (nb d'échantillons, 1..16)
                         3.0f,           // k (Hampel / z-score / deadband / a-b / Holt)
                         0.15f);         // proportion de trim (Tronquée/Winsorisée)

  // lisse = valeur brute lissée -> à calibrer ensuite
  delay(2000);
}
```

`value()` est **non destructif** (on peut l'appeler plusieurs fois), sauf pour
`DEADBAND` qui maintient un état d'hystérésis (voir §5).

---

## API détaillée

### La classe `Smoother`

```cpp
class Smoother {
  public:
    enum Type { NONE=0, SMA=1, EMA=2, MEDIAN=3, WMA=4, DEMA=5, TEMA=6,
                HAMPEL=7, TRIMMED=8, WINSORIZED=9, ZSCORE=10, KALMAN=11, CA=12,
                GAUSSIAN=13, SAVGOL=14, IIR=15, HOLT=16, ALPHABETA=17,
                DEADBAND=18, TYPE_COUNT=19 };

    static const uint8_t MAX = 16;      // taille max de la fenêtre

    Smoother();                         // état vierge
    void  reset();                      // vide le buffer + remet l'état à zéro
    void  push(float v);                // ajoute une mesure BRUTE
    float latest() const;               // dernière valeur brute
    float value(uint8_t type, uint8_t win,
                float k = 3.0f, float trim = 0.15f) const;   // valeur lissée
    void  setKalman(float q, float r);  // réglage du filtre de Kalman
    static const char* typeName(uint8_t t);   // nom lisible d'une méthode
};
```

| Méthode | Description |
|---------|-------------|
| `Smoother()` | Constructeur : buffer vide, aucun historique. |
| `reset()` | Vide le buffer et réinitialise **tout** l'état (Kalman, moyenne cumulée, deadband). À appeler après un redémarrage de logique ou un changement de capteur. |
| `push(v)` | Ajoute une mesure **brute**. Met aussi à jour l'état de Kalman et de la moyenne cumulée (indépendamment du type sélectionné). À appeler à **chaque acquisition**. |
| `latest()` | Dernière valeur brute (utile pour la méthode `NONE`). |
| `value(type, win, k, trim)` | Renvoie la valeur lissée selon la méthode. `const` (n'altère pas l'historique). |
| `setKalman(q, r)` | Règle le bruit de process `q` et le bruit de mesure `r` (utilisé par `KALMAN`). |
| `typeName(t)` | Chaîne de caractères pour l'affichage (web, logs). |

### Paramètres de `value()`

| Paramètre | Type | Défaut | Rôle |
|-----------|------|--------|------|
| `type` | `Smoother::Type` | — | Méthode de lissage. |
| `win` | `uint8_t` | — | Fenêtre : nombre d'échantillons récents utilisés (borné à `[1..min(MAX, count)]`). |
| `k` | `float` | `3.0` | Seuil générique : **nb de MAD** (Hampel), **nb d'écarts-types** (z-score), **seuil absolu** (deadband), **gain** (alpha-beta, Holt). |
| `trim` | `float` | `0.15` | Proportion d'extrêmes écartés (Tronquée) ou bornés (Winsorisée), bornée à `[0..0.4]`. |

---

## Les 5 familles de méthodes

### 1. Moyennes mobiles

| Méthode | Formule (N = fenêtre, x = valeurs brutes, récent en dernier) | Notes |
|---------|----------------------------------------------------------------|-------|
| `NONE` | `y = x[N-1]` | Réactif, aucune réduction de bruit. |
| `SMA` | `y = (Σ x[i]) / N` | Stable ; **retard ≈ N/2** échantillons. |
| `WMA` | `y = Σ (w_i · x[i]) / Σ w_i`, `w_i = i+1` | Plus réactif que SMA (le récent pèse plus). |
| `EMA` | `α = 2/(N+1)` ; `e ← α·x[i] + (1−α)·e[i−1]` | Lisse en privilégiant le récent. |
| `DEMA` | `2·EMA₁ − EMA₂` | Retard réduit par rapport à l'EMA. |
| `TEMA` | `3·EMA₁ − 3·EMA₂ + EMA₃` | Retard minimal. |
| `MEDIAN` | valeur centrale de x trié | Robuste aux pics isolés. |
| `TRIMMED` | moyenne après retrait de `trim·N` valeurs à chaque extrémité | Rejette les extrêmes, conserve la moyenne. |
| `WINSORIZED` | comme tronquée, mais les extrêmes sont **ramenés** à la valeur du centre puis moyennés | Plus doux que la tronquée. |
| `CA` | `y = (Σ toutes valeurs) / n` | Moyenne de tout l'historique : très stable, très lente. |

### 2. Filtres robustes

| Méthode | Principe | Notes |
|---------|----------|-------|
| `HAMPEL` | médiane `m` + écart absolu médian `MAD` ; si `abs(x_last − m) > k·MAD`, on remplace par `m` | Excellent contre les **pics ponctuels**. |
| `ZSCORE` | moyenne + écart-type ; ne moyenne que les points dans `k·σ` | Rejette les valeurs « anormales ». |
| `KALMAN` | filtre récursif 1D : `P += Q ; K = P/(P+R) ; x += K(v−x) ; P *= (1−K)` | Optimal bruit + modèle ; réglé par `Q`/`R`. |
| `ALPHABETA` | position + vitesse : `pred = x + v ; res = z − pred ; x = pred + α·res ; v = v + β·res` | Robuste, suit une rampe ; `α` déduit de `win`, `β` de `k`. |

### 3. Filtres de signal (DSP)

| Méthode | Principe | Notes |
|---------|----------|-------|
| `GAUSSIAN` | moyenne pondérée par un noyau gaussien `exp(−d²/2σ²)`, `σ = N/3` | Lissage très doux. |
| `SAVGOL` | ajustement d'une **parabole** (moindres carrés) sur la fenêtre, évaluée au point le plus récent | Lisse en **préservant les pics** (contrairement aux moyennes). |
| `IIR` | filtre passe-bas **2e ordre** (biquad type Butterworth, `Q=0.707`), fréquence de coupure liée à `N` | Rejette le bruit haute fréquence. |

### 4. Lissage orienté tendance

| Méthode | Principe | Notes |
|---------|----------|-------|
| `HOLT` | double lissage exponentiel (niveau + tendance) : `level ← α·x + (1−α)(level+trend) ; trend ← β(level−prev) + (1−β)·trend` | Suit une dérive ; `α` déduit de `win`, `β` de `k`. |

### 5. Anti-rebond (hystérésis)

| Méthode | Principe | Notes |
|---------|----------|-------|
| `DEADBAND` | la sortie ne change que si `abs(x − sortie) > k` | Évite le « scintillement ». Seuil `k` **en unités du capteur**. Agit sur la **dernière valeur poussée** (brute) : pour anti-rebonder une valeur **lissée**, la pousser dans une **2ᵉ instance** de `Smoother`. |

---

## Paramètres : comment les régler

- **`win` (fenêtre)** : plus elle est grande, plus le lissage est fort, mais plus
  la réponse est lente. Typiques : `4` (réactif), `8` (équilibré), `16` (très
  lissé). La fenêtre s'adapte automatiquement si l'historique est incomplet.
- **`k`** :
  - `HAMPEL` : nombre de MAD ; `ZSCORE` : nombre d'écarts-types. `3` = standard ; `2` = plus
    strict (rejette plus), `4-5` = plus permissif.
  - `DEADBAND` : seuil **absolu** dans l'unité de la mesure (ex. `0.1` °C). Un
    `k` à 3 bloquerait pratiquement tout : réglez-le à la précision voulue.
  - `ALPHABETA` / `HOLT` : sert de gain (`β`). Plus grand = plus réactif.
- **`trim`** : proportion par extrémité pour `TRIMMED`/`WINSORIZED`. `0.15` ≈
  15 % de chaque côté. Borné à `0.4`.
- **Kalman `Q`/`R`** (`setKalman`) :
  - `Q` (bruit de process) petit → très lisse/lent ; grand → plus réactif.
  - `R` (bruit de mesure) grand → fait plus confiance au modèle (lisse).
  - Défauts `Q=1e-4`, `R=1e-2`.

---

## Guide de choix par situation

| Situation | Méthode conseillée |
|-----------|--------------------|
| Bruit blanc modéré | `SMA`, `EMA`, `WMA`, `Kalman` |
| Pics/parasites ponctuels | `Médiane`, `Hampel`, `Z-score` |
| Valeurs aberrantes fréquentes | `Hampel`, `Tronquée`, `Winsorisée` |
| Besoin de préserver les pics utiles | `Savitzky-Golay` |
| Bruit haute fréquence (DSP) | `Gaussien`, `IIR` |
| Grandeur qui dérive (rampe) | `Holt`, `Alpha-beta`, `Kalman` |
| Affichage/Domoticz qui « saute » | `Hystérésis` (en complément) |
| Comparaison/recherche | `Aucun` |

---

## Formules complètes

Notations : `N` = fenêtre, `x[0]` = plus ancien, `x[N-1]` = plus récent.

- **SMA** : `y = (1/N) Σ x[i]`
- **WMA** : `y = (Σ (i+1)·x[i]) / (Σ (i+1))`
- **EMA** (et étages suivants) : `α = 2/(N+1)`, `e[0] = x[0]`, `e[i] = α·x[i] + (1−α)·e[i−1]`, `y = e[N-1]`
- **DEMA** : `y = 2·E1 − E2`, où `E1 = EMA(x)`, `E2 = EMA(E1)`
- **TEMA** : `y = 3·E1 − 3·E2 + E3`
- **Médiane** : `y = x_sorted[N/2]` (ou moyenne des deux centraux si N pair)
- **Tronquée** : `t = round(trim·N)` ; `y = moyenne(x_sorted[t .. N-1-t])`
- **Winsorisée** : `t = round(trim·N)` ; les valeurs `< x_sorted[t]` sont remontées à `x_sorted[t]`, celles `> x_sorted[N-1-t]` rabaissées, puis `y = moyenne`
- **Cumulée** : `y = Σ_historique / n`
- **Hampel** : `m = médiane(x)`, `MAD = médiane(|x−m|)` ; `y = m` si `|x[N-1]−m| > k·MAD` sinon `x[N-1]`
- **Z-score** : `μ`, `σ` sur la fenêtre (σ **population**, ÷N) ; `y = moyenne` des `x[i]` tels que `|x[i]−μ| ≤ k·σ`
- **Kalman** : voir §2 ; `Q`, `R` réglables
- **Alpha-beta** : voir §2
- **Gaussien** : `w_i = exp(−0.5·((N−1−i)/σ)²)`, `σ = N/3` ; `y = Σ w_i·x[i] / Σ w_i`
- **Savitzky-Golay** : moindres carrés paraboliques `a + b·i + c·i²`, `y = a + b·(N−1) + c·(N−1)²`
- **IIR (Butterworth 2e ordre)** : biquad passe-bas, `f_c = 1/(2N)` échantillons, `Q = 0.7071`
- **Holt** : voir §4
- **Deadband** : `y` ne change que si `|x − y_précédent| > k`

---

## Mémoire, CPU et performance

- **RAM par instance** : buffer `float _buf[16]` (64 o) + `count/head` + état Kalman
  (4 `float`) + moyenne cumulée (`double` + `uint32_t`) + état deadband ≈ **~104 o**
  par voie. Pour 5 voies : ~520 o.
- **Aucune allocation dynamique** : tout est statique dans l'objet.
- **CPU** : l'essentiel est O(win) avec `win ≤ 16`. Les méthodes à tri
  (`Médiane`, `Tronquée`, `Winsorisée`, `Hampel`) font un tri à bulles O(win²)
  (≤ 256 comparaisons) : négligeable.
- **Coût flash** : variable selon les méthodes activées. Les fonctions
  mathématiques (`sqrtf`, `expf`, `sinf`, `cosf`) pèsent surtout en flash, pas en
  RAM.
- **Flottants** : l'ESP8266 n'a pas de FPU ; `float` est natif mais les doubles
  (Savitzky-Golay, moyenne cumulée) sont plus lents. À `2 s` par acquisition,
  c'est négligeable.

---

## Intégration dans un projet (brut → lissage → calibration)

```
   capteur ──► valeur BRUTE ──► Smoother::push()
                                      │
                                      ▼
                              Smoother::value()  ──► moyenne brute
                                      │
                                      ▼
                              calibration (offset/gain) ──► valeur affichée
```

- Un `Smoother` **par voie** (température, humidité, pression, …).
- `push()` à chaque acquisition ; `value()` au moment de calculer la mesure
  publiée.
- La calibration s'applique **après** le lissage.

---

## Limites et pièges

- **Ne pas moyenner du calibré** : toujours lisser le brut.
- **`DEADBAND` est à état** : `value(DEADBAND, …)` modifie l'état interne
  (mutable). Ne l'appeler qu'**une fois par acquisition** ; pour un changement de
  configuration, faire `reset()`.
- **`DEADBAND` agit sur le brut** : il compare la **dernière valeur poussée** au
  seuil, pas la valeur lissée. Pour anti-rebonder un lissage, utiliser un **2ᵉ**
  `Smoother` alimenté par la valeur lissée.
- **`win` > échantillons disponibles** : la fenêtre est automatiquement réduite à
  `count` (utile au démarrage).
- **`k` mal adapté à `DEADBAND`** : le seuil est **absolu** (pas un multiple
  d'écart-type) ; le régler dans l'unité de la mesure.
- **Savitzky-Golay / IIR sur très peu de points** : nécessitent une fenêtre
  suffisante (`≥ 3` pour Savitzky-Golay) ; sinon `latest()` est renvoyé.
- **IIR, DEMA, TEMA, Holt** ont une réponse transitoire après `reset()` : laisser
  quelques échantillons avant d'exploiter la valeur. (`IIR` est initialisé sur la
  première valeur, ce qui **réduit** ce transitoire sans le supprimer.)

---

## Exemples

### Comparer deux méthodes

```cpp
Smoother sma, ham;

void loop() {
  float x = lireCapteur();
  sma.push(x);
  ham.push(x);

  Serial.printf("brut=%.3f  SMA(8)=%.3f  Hampel(k=3)=%.3f\n",
                x,
                sma.value(Smoother::SMA, 8),
                ham.value(Smoother::HAMPEL, 8, 3.0f));
  delay(2000);
}
```

### Filtre de Kalman réglé

```cpp
Smoother kf;
kf.setKalman(1e-4f, 1e-2f);           // lisse et peu de retard

void loop() {
  kf.push(lireCapteur());
  float y = kf.value(Smoother::KALMAN, 1);   // win ignoré pour Kalman
  delay(2000);
}
```

### Anti-rebond d'affichage (après un lissage)

```cpp
Smoother lis;      // 1er : lissage
Smoother aff;      // 2e  : anti-rebond (état propre)

lis.push(x);
float lisse = lis.value(Smoother::EMA, 8);

aff.push(lisse);                                          // on anti-rebonde la valeur LISSÉE
float affiche = aff.value(Smoother::DEADBAND, 1, 0.1f);   // seuil 0,1 unité
```

> `DEADBAND` agit toujours sur la **dernière valeur poussée** (`latest()`). Pour
> anti-rebonder une valeur lissée, il faut donc un **second** `Smoother` alimenté
> par cette valeur lissée.

### Réinitialiser après changement de capteur

```cpp
sm.reset();
sm.setKalman(1e-4f, 1e-2f);
```

---

## Licence

GPL-3.0. Header-only, aucune dépendance externe.
