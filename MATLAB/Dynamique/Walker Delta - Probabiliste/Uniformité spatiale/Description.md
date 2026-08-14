# Description des scripts MATLAB

Ce document décrit les fichiers MATLAB, ainsi que leurs fonctions locales.

---

## `analysis_temp.m`

### Objectif

Simule un réseau LEO dynamique à inclinaison commune fixe, avec des positions initiales uniformes en surface dans la bande de latitude accessible :

\[
|\varphi|\leq i.
\]

Le script :

- tire le nombre de satellites selon un processus de Poisson sur la bande ;
- génère des positions initiales uniformes en surface ;
- reconstruit, pour chaque satellite, une orbite circulaire d’inclinaison imposée passant par sa position initiale ;
- fait évoluer les satellites avec une vitesse angulaire orbitale commune ;
- construit le graphe des liens intersatellites à chaque instant ;
- calcule \(\beta_0\), \(\beta_1\), la plus grande composante et le nombre de liens ;
- étudie la distribution des tailles de composantes ;
- compare les composantes de tailles 1, 2 et 3 à des approximations théoriques ;
- construit la suite zigzag par unions ;
- sauvegarde toutes les données nécessaires au calcul des barcodes.

La suite zigzag construite est

\[
G_1
\longrightarrow
G_1\cup G_2
\longleftarrow
G_2
\longrightarrow
G_2\cup G_3
\longleftarrow
G_3
\longrightarrow\cdots
\]

### Type

Script principal avec une fonction locale :

- `walker_delta_positions`

### Entrées du script

Le script ne possède pas d’arguments d’entrée. Les paramètres sont définis directement dans le fichier.

### Paramètres physiques

| Variable | Valeur | Description |
|---|---:|---|
| `R_earth` | `6371` km | Rayon terrestre. |
| `h` | `550` km | Altitude orbitale. |
| `R` | `R_earth + h` | Rayon orbital. |
| `mu` | `398600` km³/s² | Paramètre gravitationnel terrestre. |
| `omega` | \(\sqrt{\mu/R^3}\) | Vitesse angulaire orbitale. |

### Paramètres du modèle spatial

| Variable | Valeur | Description |
|---|---:|---|
| `inc_deg` | `58` degrés | Inclinaison orbitale commune. |
| `inc` | `deg2rad(inc_deg)` | Inclinaison en radians. |
| `lambda` | `4e-7` sat/km² | Densité spatiale dans la bande accessible. |
| `surface_band` | \(4\pi R^2\sin i\) | Aire de la bande sphérique accessible. |
| `N` | aléatoire | Nombre de satellites tiré selon une loi de Poisson. |

Le nombre de satellites suit

\[
N\sim\mathcal P\!\left(
\lambda\,4\pi R^2\sin i
\right).
\]

### Paramètres temporels et de connexion

| Variable | Valeur | Description |
|---|---:|---|
| `dmax` | `1500` km | Distance maximale de connexion. |
| `dt` | `20` s | Pas temporel. |
| `Tmax` | `12000` s | Durée totale simulée. |
| `time_values` | `0:dt:Tmax` | Instants de simulation. |

### Sorties principales

| Variable | Description |
|---|---|
| `Positions` | Positions temporelles des satellites. |
| `Adjacency` | Matrices d’adjacence temporelles. |
| `beta0` | Nombre de composantes connexes. |
| `beta1_graph` | Nombre cyclomatique. |
| `largest_component` | Taille de la plus grande composante. |
| `num_edges` | Nombre de liens. |
| `component_size_counts` | Nombre de composantes par taille et par instant. |
| `mean_component_count_by_size` | Nombre moyen de composantes par taille. |
| `component_size_fraction` | Distribution normalisée des tailles. |
| `mean_component_size_time` | Taille moyenne d’une composante à chaque instant. |
| `mean_component_size` | Moyenne temporelle correspondante. |
| `N1_emp_time`, `N2_emp_time`, `N3_emp_time` | Petites composantes empiriques temporelles. |
| `N1_emp_mean`, `N2_emp_mean`, `N3_emp_mean` | Moyennes empiriques. |
| `N1_theory`, `N2_theory`, `N3_theory` | Approximations théoriques. |
| `relative_error_N123` | Erreurs relatives. |
| `ZigzagAdjacency` | Matrices d’adjacence du zigzag. |
| `ZigzagLabels` | Labels du zigzag. |
| `beta0_zigzag` | Nombre de composantes sur le zigzag. |
| `beta1_zigzag_graph` | Nombre cyclomatique sur le zigzag. |
| `largest_component_zigzag` | Plus grande composante sur le zigzag. |
| `num_edges_zigzag` | Nombre de liens sur le zigzag. |

### Figures produites

Le script produit notamment :

1. \(\beta_0(t)\) ;
2. \(\beta_1^{\mathrm{graphe}}(t)\) ;
3. fraction de satellites dans la plus grande composante ;
4. nombre de liens ;
5. nombre moyen de composantes par taille ;
6. distribution normalisée des tailles ;
7. comparaison temporelle de \(N_1\), \(N_2\), \(N_3\) ;
8. comparaison moyenne empirique/théorique ;
9. erreurs relatives sur \(N_1\), \(N_2\), \(N_3\).

### Fichier sauvegardé

```matlab
analysis_temp_results.mat
```

---

## `anim_3D.m`

### Objectif

Anime en trois dimensions le réseau LEO Walker-Delta à uniformité spatiale initiale dans la bande accessible.

Le script montre :

- les satellites ;
- les liens intersatellites ;
- un sous-ensemble de plans orbitaux ;
- le temps courant ;
- le nombre de satellites ;
- le nombre de liens.

### Type

Script principal avec une fonction locale :

- `walker_delta_positions`

### Entrées du script

| Variable | Valeur | Description |
|---|---:|---|
| `R_earth` | `6371` km | Rayon terrestre. |
| `h` | `550` km | Altitude. |
| `inc_deg` | `90` degrés | Inclinaison commune. |
| `lambda` | `4e-7` sat/km² | Densité dans la bande accessible. |
| `dmax` | `1500` km | Portée des liens. |
| `dt` | `30` s | Pas temporel. |
| `Tmax` | `6000` s | Durée simulée. |

### Génération initiale

Le script suit la même méthode que `analysis_temp.m` :

1. tirage uniforme de la longitude ;
2. tirage uniforme de \(\sin(\text{latitude})\) dans la bande ;
3. choix aléatoire de la branche ascendante ou descendante ;
4. reconstruction de \(u_0\) ;
5. calcul du RAAN \(\Omega\).

### Sorties principales

Le script ne sauvegarde aucun fichier.

Les principales variables du workspace sont :

| Variable | Description |
|---|---|
| `positions0` | Positions initiales. |
| `positions_t` | Positions au dernier instant. |
| `longitude0`, `latitude0` | Coordonnées initiales tirées. |
| `Omega` | RAAN reconstruits. |
| `u0` | Arguments de latitude initiaux. |
| `A` | Dernière matrice d’adjacence. |
| `E` | Nombre de liens au dernier instant. |
| `sat_handle` | Objet graphique des satellites. |
| `link_handle` | Objet graphique des liens. |

### Affichage des plans orbitaux

Un plan individuel est associé à chaque satellite. Pour conserver une figure lisible, le script en affiche au maximum 40.

### Figure produite

Une animation 3D interactive avec :

- satellites ;
- liens ;
- plans orbitaux en pointillés ;
- titre dynamique.

---

## `barcodes.m`

### Objectif

Calcule le barcode zigzag en homologie \(H_0\) à partir des graphes produits par `analysis_temp.m`.

Le script :

- charge la suite zigzag ;
- calcule les composantes connexes de chaque objet ;
- construit les applications induites en \(H_0\) ;
- calcule les intervalles de persistance ;
- convertit les indices en temps physiques ;
- charge les séries théoriques temporelles de fusion et de rupture ;
- calcule leur moyenne ;
- construit un modèle exponentiel de survie ;
- affiche le barcode et l’histogramme des durées ;
- sauvegarde les résultats.

### Type

Script principal avec neuf fonctions locales :

- `first_existing_file`
- `build_H0_map`
- `zigzag_barcode_from_module_mod2`
- `filtration_quotient_dims`
- `gf2_preimage`
- `gf2_col_basis`
- `gf2_rank`
- `gf2_null`
- `gf2_rref`

### Fichier principal d’entrée

```matlab
analysis_temp_results.mat
```

Variables chargées :

| Variable | Description |
|---|---|
| `ZigzagAdjacency` | Matrices d’adjacence du zigzag. |
| `ZigzagLabels` | Labels entiers et demi-entiers. |
| `time_values` | Instants des graphes réels. |
| `N` | Nombre de satellites. |
| `R` | Rayon orbital. |
| `lambda` | Densité satellitaire. |
| `dmax` | Portée des liens. |
| `dt` | Pas temporel. |
| `omega` | Vitesse angulaire orbitale. |
| `inc_deg` | Inclinaison. |
| `P` | Nombre de plans orbitaux. |

Des valeurs de secours sont prévues lorsque `omega`, `inc_deg` ou `P` sont absents.

### Sorties principales

| Variable | Description |
|---|---|
| `intervals` | Intervalles de persistance. |
| `birth_index` | Indices de naissance. |
| `death_index` | Indices de mort. |
| `birth_time` | Temps de naissance. |
| `death_time` | Temps de mort. |
| `lifetimes` | Durées de vie. |
| `ZigzagTime` | Temps associés aux objets zigzag. |
| `h0_dims` | Dimensions successives de \(H_0\). |
| `Lvals` | Durées positives distinctes. |
| `survival` | Survie empirique. |
| `p_merge_th_t` | Série théorique de fusion. |
| `p_break_th_t` | Série théorique de rupture. |
| `p_merge_mean_th` | Moyenne temporelle de fusion. |
| `p_break_mean_th` | Moyenne temporelle de rupture. |
| `p_change_mean_th` | Moyenne entre les deux probabilités. |
| `dt_model` | Pas temporel retenu. |
| `tau_th` | Temps caractéristique théorique. |
| `survival_th` | Survie exponentielle théorique. |

### Figures produites

1. survie empirique et modèle exponentiel ;
2. barcode zigzag \(H_0\) ;
3. histogramme des durées de vie.

### Fichier sauvegardé

```matlab
barcodes_results.mat
```

---

---

# Sous-dossier `Paramètres`

## `Paramètres/nb_bordure_composante.m`

### Objectif

Estime le nombre moyen de satellites situés sur la bordure d’une composante de taille \(n\), ainsi que le nombre de satellites participant effectivement au balayage.

Le modèle est

\[
N_{\mathrm{bord}}
=
K_{\mathrm{bord}}\sqrt n,
\qquad
K_{\mathrm{bord}}
=
2\sqrt{\pi}\,\ell\sqrt{\lambda},
\]

puis

\[
N_{\mathrm{balayage}}
=
f_{\mathrm{balayage}}N_{\mathrm{bord}}.
\]

### Signature

```matlab
[N_bord, N_balayage, K_bord] = ...
    nb_bordure_composante(n, lambda, ell, facteur_balayage)
```

### Entrées

| Argument | Description |
|---|---|
| `n` | Taille de la composante. |
| `lambda` | Densité surfacique en points/km². |
| `ell` | Épaisseur effective de la bordure en kilomètres. |
| `facteur_balayage` | Fraction de la bordure contribuant au balayage, égale à \(1/\pi\) par défaut. |

### Sorties

| Sortie | Description |
|---|---|
| `N_bord` | Nombre moyen de points en bordure. |
| `N_balayage` | Nombre moyen de points contribuant au balayage. |
| `K_bord` | Coefficient tel que \(N_{\mathrm{bord}}=K_{\mathrm{bord}}\sqrt n\). |

Pour

\[
\ell=\frac{1}{\sqrt{\lambda}}
\quad\text{et}\quad
f_{\mathrm{balayage}}=\frac{1}{\pi},
\]

le code donne

\[
N_{\mathrm{balayage}}
=
\frac{2}{\sqrt{\pi}}\sqrt n,
\qquad
\chi_{\mathrm{merge}}
=
\frac{2}{\sqrt{\pi n}}.
\]

---

## `Paramètres/lambda_eff_th.m`

### Objectif

Calcule une densité satellitaire effective temporelle pour le Walker-Delta à uniformité spatiale initiale.

Pour chaque strate \(m\),

\[
\lambda_m(t)
=
\frac{N_m(t)}{A_m},
\]

avec

\[
A_m
=
4\pi R^2
\left[
\sin(\phi_{m+1})-\sin(\phi_m)
\right].
\]

La densité effective est pondérée par le nombre de satellites :

\[
\lambda_{\mathrm{eff}}(t)
=
\frac{
\sum_m N_m(t)\lambda_m(t)
}{
\sum_m N_m(t)
}.
\]

### Fichier d’entrée

```matlab
Nombre liens/liens_quadrature_results.mat
```

### Variables nécessaires

| Variable | Description |
|---|---|
| `time_values` | Instants de simulation. |
| `R` | Rayon orbital. |
| `dmax` | Portée de communication. |
| `inc` | Inclinaison. |
| `mean_n_sat_strate_t` ou `N_strate_theory_time` | Nombre de satellites par strate et par instant. |

### Sorties principales

| Variable | Description |
|---|---|
| `lambda_eff_t` | Densité effective temporelle. |
| `lambda_eff_mean` | Moyenne temporelle. |
| `lambda_eff_min`, `lambda_eff_max` | Valeurs minimale et maximale. |
| `lambda_strate_t` | Densité locale de chaque strate. |
| `N_strate_t` | Nombre de satellites par strate. |
| `A_strates` | Aires des strates. |
| `lambda_global_band` | Densité globale sur la bande. |
| `phi_in`, `phi_out`, `phi_mid` | Limites et centres des strates. |

### Figures produites

1. densité effective en fonction du temps ;
2. densité moyenne et nombre moyen de satellites par strate ;
3. carte temporelle des densités par strate.

### Fichier sauvegardé

```matlab
lambda_eff_th_results.mat
```

Le message final mentionne `lambda_effective_temp_delta_spatial_results.mat`, mais le fichier réellement créé est `lambda_eff_th_results.mat`.

---

## `Paramètres/plink_t_phi.m`

### Objectif

Calcule et compare la probabilité locale de lien

\[
p_{\mathrm{link}}(t,\phi)
=
P\!\left(d_{12}(t)\le d_{\max}\mid \Phi_1(t)=\phi\right)
\]

dans le modèle Walker Delta à uniformité spatiale initiale.

Le script resimule plusieurs constellations indépendantes et compare la mesure empirique à un modèle théorique fondé sur la loi temporelle de phase orbitale. La latitude conditionnelle est traitée en tenant compte des deux branches orbitales compatibles avec une même latitude.

### Type

Script principal avec plusieurs fonctions locales de génération orbitale, de calcul du noyau géométrique de lien et d’échantillonnage conditionnel.

### Paramètres principaux

| Variable | Description |
|---|---|
| `N` | Nombre de satellites. |
| `R` | Rayon orbital. |
| `dmax` | Distance maximale de connexion. |
| `inc` | Inclinaison orbitale. |
| `omega` | Vitesse angulaire orbitale. |
| `dt`, `Tmax` | Grille temporelle. |
| `n_iterations` | Nombre de réalisations empiriques. |
| `n_phi_bins_emp` | Nombre de tranches de latitude empiriques. |
| `n_phi_bins_th` | Résolution de la grille théorique fine. |

### Théorie locale

À l’instant \(t\), la densité de phase est construite à partir de

\[
f_u(u,t)\propto |\cos(u-\omega t)|.
\]

Pour une latitude donnée, les deux phases compatibles sont pondérées conditionnellement, puis la probabilité de lien est obtenue après intégration sur la phase du second satellite et sur la différence de RAAN.

La valeur théorique comparable à l’empirique est moyennée sur toute la tranche de latitude :

\[
p_{\mathrm{link},b}^{\mathrm{th}}(t)
=
\frac{
\int_b p_{\mathrm{link}}^{\mathrm{th}}(t,\phi)
f_\Phi(t,\phi)\,d\phi
}{
\int_b f_\Phi(t,\phi)\,d\phi
}.
\]

### Sorties principales

| Variable | Description |
|---|---|
| `p_link_emp_iterations` | Probabilité locale empirique pour chaque réalisation. |
| `p_link_emp_mean` | Probabilité locale empirique agrégée. |
| `p_link_emp_iteration_sem` | Erreur standard entre réalisations. |
| `p_link_th_fine` | Probabilité locale théorique sur la grille fine. |
| `p_link_th_on_emp` | Théorie moyennée sur les tranches empiriques. |
| `p_link_th_global` | Probabilité globale théorique. |
| `p_link_emp_global_mean` | Probabilité globale empirique moyenne. |
| `f_phi_th_fine` | Loi théorique de latitude dépendant du temps. |
| `satellite_count_emp_iterations` | Nombre de satellites par tranche et réalisation. |
| `degree_sum_emp_iterations` | Somme des degrés par tranche et réalisation. |

### Fichier sauvegardé

```matlab
plink_t_phi_results.mat
```

---

## `Paramètres/N_t_phi.m`

### Objectif

Calcule le nombre local de satellites

\[
N_b(t)
\]

dans chaque tranche de latitude et le compare à la théorie issue de la loi \(f_\Phi(t,\phi)\).

La densité continue vaut

\[
n_\phi^{\mathrm{th}}(t,\phi)
=
Nf_\Phi(t,\phi),
\]

et, pour une tranche \(b\),

\[
N_b^{\mathrm{th}}(t)
=
N\int_b f_\Phi(t,\phi)\,d\phi.
\]

### Fichier d’entrée

```matlab
plink_t_phi_results.mat
```

### Sorties principales

| Variable | Description |
|---|---|
| `satellite_count_th` | Nombre théorique de satellites par tranche et par instant. |
| `satellite_count_emp_iterations` | Comptages empiriques par réalisation. |
| `satellite_count_emp` | Moyenne empirique. |
| `satellite_count_emp_sem` | Erreur standard empirique. |
| `satellite_density_phi_th_fine` | Densité théorique continue en latitude. |
| `satellite_density_phi_emp` | Densité empirique par radian. |
| `N_total_th` | Nombre total théorique reconstruit. |
| `N_total_emp` | Nombre total empirique moyen. |
| `rmse_grid`, `mae_grid`, `bias_grid` | Diagnostics théorie/empirique. |

### Fichier sauvegardé

```matlab
N_t_phi_results.mat
```

---

## `Paramètres/lambda_t_phi.m`

### Objectif

Calcule la densité surfacique locale dépendant du temps et de la latitude :

\[
\lambda(t,\phi)
\quad [\mathrm{satellites/km^2}].
\]

La théorie continue est

\[
\lambda^{\mathrm{th}}(t,\phi)
=
\frac{
Nf_\Phi(t,\phi)
}{
2\pi R^2\cos\phi
}.
\]

Pour une tranche \(b\),

\[
\lambda_b^{\mathrm{th}}(t)
=
\frac{
N\,P(\Phi(t)\in b)
}{
A_b
},
\]

avec

\[
A_b
=
2\pi R^2
\left[
\sin(\phi_b^+)-\sin(\phi_b^-)
\right].
\]

### Fichier d’entrée

```matlab
plink_t_phi_results.mat
```

### Sorties principales

| Variable | Description |
|---|---|
| `area_bin` | Aire exacte de chaque tranche de latitude. |
| `lambda_th_fine` | Densité locale théorique sur la grille fine. |
| `lambda_bin_th` | Densité théorique moyenne par tranche. |
| `lambda_bin_emp_iterations` | Densité empirique par réalisation. |
| `lambda_bin_emp` | Densité empirique moyenne. |
| `lambda_bin_emp_sem` | Erreur standard empirique. |
| `N_th_reconstructed` | Nombre de satellites reconstruit depuis la densité théorique. |
| `N_emp_reconstructed` | Nombre reconstruit depuis la densité empirique. |
| `rmse_grid`, `mae_grid`, `bias_grid` | Erreurs locales. |

### Fichier sauvegardé

```matlab
lambda_t_phi_results.mat
```

---

## `Paramètres/eta_sweep_t_phi.m`

### Objectif

Calcule le facteur local de redondance spatiale

\[
\eta_{\mathrm{sweep}}(t,\phi)
\]

à partir de la densité locale \(\lambda(t,\phi)\).

Le modèle est

\[
\eta_{\mathrm{sweep}}(t,\phi)
=
\exp\!\left[
-\lambda(t,\phi)
A_{\mathrm{inter}}(d_{\max})
\right],
\]

avec

\[
A_{\mathrm{inter}}(d_{\max})
=
\left(
\frac{2\pi}{3}
-\frac{\sqrt3}{2}
\right)d_{\max}^2.
\]

### Fichier d’entrée

```matlab
lambda_t_phi_results.mat
```

Le script récupère également `dmax` depuis le fichier source de `lambda_t_phi_results.mat` ou depuis `plink_t_phi_results.mat` si nécessaire.

### Particularité statistique

Comme l’exponentielle est non linéaire, la valeur empirique principale est calculée par

\[
\mathbb E_r
\left[
e^{-A_{\mathrm{inter}}\lambda^{(r)}(t,\phi)}
\right]
\]

et non par l’application de l’exponentielle à la densité empirique moyenne.

### Sorties principales

| Variable | Description |
|---|---|
| `eta_sweep_th_fine` | Théorie sur la grille fine. |
| `eta_sweep_bin_th` | Théorie par tranche de latitude. |
| `eta_sweep_emp_iterations` | Valeurs empiriques pour chaque réalisation. |
| `eta_sweep_emp` | Moyenne empirique. |
| `eta_sweep_emp_sem` | Erreur standard empirique. |
| `eta_sweep_from_mean_lambda_emp` | Valeur obtenue à partir de la densité empirique moyenne. |
| `jensen_difference` | Écart dû à la non-linéarité de l’exponentielle. |
| `rmse_grid`, `mae_grid`, `bias_grid` | Diagnostics théorie/empirique. |

### Fichier sauvegardé

```matlab
eta_sweep_t_phi_results.mat
```

---

---

# Sous-dossier `Paramètres/Betti`

## `Paramètres/Betti/N1_t_phi.m`

### Objectif

Calcule et compare le nombre local de satellites isolés

\[
N_1(t,\phi)
\]

dans chaque tranche de latitude.

Un satellite est isolé lorsque son degré est nul. Le calcul théorique tient explicitement compte des deux branches orbitales compatibles avec une latitude donnée.

### Fichier d’entrée

```matlab
../plink_t_phi_results.mat
```

### Théorie

La probabilité correcte d’isolement est

\[
p_{\mathrm{iso}}(t,\phi)
=
w_+(t,\phi)
\left[1-p_+(t,\phi)\right]^{N-1}
+
w_-(t,\phi)
\left[1-p_-(t,\phi)\right]^{N-1},
\]

où \(p_+\) et \(p_-\) sont les probabilités de lien conditionnelles aux deux branches orbitales.

Le nombre théorique d’isolés dans une tranche \(b\) est

\[
N_{1,b}^{\mathrm{th}}(t)
=
N
\int_b
f_\Phi(t,\phi)
p_{\mathrm{iso}}(t,\phi)
\,d\phi.
\]

### Sorties principales

| Variable | Description |
|---|---|
| `isolated_count_th` | Nombre théorique d’isolés par tranche et par instant. |
| `isolated_count_emp_iterations` | Comptages empiriques par réalisation. |
| `isolated_count_emp` | Moyenne empirique. |
| `isolated_count_emp_sem` | Erreur standard empirique. |
| `p_iso_th_fine` | Probabilité théorique fine d’isolement. |
| `p_iso_emp` | Probabilité empirique d’isolement. |
| `isolated_total_th` | Nombre total théorique d’isolés. |
| `isolated_total_emp` | Nombre total empirique moyen. |

### Fichier sauvegardé

```matlab
N1_t_phi_results.mat
```

---

## `Paramètres/Betti/N2_t_phi.m`

### Objectif

Calcule le nombre local de dimères isolés \(N_2(t,\phi)\) par un modèle géométrique, puis le compare à des graphes resimulés.

La théorie utilise

\[
C_{2,b}^{\mathrm{th}}(t)
=
\frac{N(N-1)}{2}
P(\Phi_1\in b)
\,
\mathbb E
\left[
\mathbf 1_{\{X_1\sim X_2\}}
(1-q_2)^{N-2}
\mid
\Phi_1\in b
\right],
\]

où \(q_2\) est la probabilité qu’un troisième satellite soit relié à au moins l’un des deux sommets.

### Fichier d’entrée

```matlab
../plink_t_phi_results.mat
```

### Méthode

L’espérance géométrique est estimée par Monte-Carlo imbriqué. Le second satellite est conditionné à être lié au premier afin d’éviter un événement rare, puis une estimation non biaisée de la probabilité d’exclusion des \(N-2\) satellites restants est utilisée.

### Sorties principales

| Variable | Description |
|---|---|
| `component_count_th` | Nombre théorique de dimères par tranche et par instant. |
| `component_count_emp_iterations` | Comptages empiriques par réalisation. |
| `component_count_emp` | Moyenne empirique. |
| `component_count_emp_sem` | Erreur standard empirique. |
| `component_total_th` | Nombre total théorique de dimères. |
| `component_total_emp` | Nombre total empirique moyen. |
| `mean_q2_eval` | Valeur moyenne locale de \(q_2\). |
| `rmse_grid`, `mae_grid`, `bias_grid` | Diagnostics locaux. |

### Fichier sauvegardé

```matlab
N2_t_phi_results.mat
```

---

## `Paramètres/Betti/N3_t_phi.m`

### Objectif

Calcule le nombre local de trimères isolés \(N_3(t,\phi)\) par un modèle géométrique et le compare à des simulations complètes.

Le modèle théorique est

\[
C_{3,b}^{\mathrm{th}}(t)
=
\frac{N}{3}
P(\Phi_1\in b)
\binom{N-1}{2}
\,
\mathbb E
\left[
\mathbf 1_{\{G_3\ \mathrm{connexe}\}}
(1-q_3)^{N-3}
\mid
\Phi_1\in b
\right].
\]

### Fichier d’entrée

```matlab
../plink_t_phi_results.mat
```

### Méthode

Pour réduire la rareté du tirage :

1. \(X_2\) est conditionné à être lié à \(X_1\) ;
2. \(X_3\) est conditionné à être lié à \(X_1\) ou \(X_2\) ;
3. une correction d’importance géométrique est appliquée ;
4. la probabilité d’exclusion des autres satellites est estimée sur un pool indépendant.

### Sorties principales

| Variable | Description |
|---|---|
| `component_count_th` | Nombre théorique de trimères par tranche. |
| `component_count_emp_iterations` | Comptages empiriques par réalisation. |
| `component_count_emp` | Moyenne empirique. |
| `component_count_emp_sem` | Erreur standard empirique. |
| `component_total_th` | Nombre total théorique de trimères. |
| `component_total_emp` | Nombre total empirique moyen. |
| `mean_q3_eval` | Valeur moyenne locale de \(q_3\). |
| `rmse_grid`, `mae_grid`, `bias_grid` | Diagnostics locaux. |

### Fichier sauvegardé

```matlab
N3_t_phi_results.mat
```

---

## `Paramètres/Betti/betti_t_phi.m`

### Objectif

Construit une approximation locale et temporelle de \(\beta_0(t,\phi)\) à partir des petites composantes :

\[
\beta_0^{\mathrm{th}}(t)
\approx
N_1^{\mathrm{th}}(t)
+
N_2^{\mathrm{th}}(t)
+
N_3^{\mathrm{th}}(t)
+
2.
\]

Les deux composantes macroscopiques sont réparties localement selon la masse résiduelle de satellites n’appartenant ni aux isolés, ni aux dimères, ni aux trimères.

### Fichiers d’entrée

```matlab
N1_t_phi_results.mat
N2_t_phi_results.mat
N3_t_phi_results.mat
```

### Mesure empirique

La valeur empirique de \(\beta_0\) est recalculée directement sur des graphes complets et prend en compte toutes les tailles de composantes.

Une composante de taille \(s\) contribue localement pour \(1/s\) dans la tranche de chacun de ses membres. Ainsi, la somme des contributions locales d’une composante vaut exactement 1.

### Sorties principales

| Variable | Description |
|---|---|
| `beta0_th` | Approximation théorique locale \(\beta_0(t,\phi)\). |
| `beta0_macro_th` | Contribution locale des deux composantes macroscopiques. |
| `beta0_emp_true` | Mesure empirique complète de \(\beta_0(t,\phi)\). |
| `beta0_emp_true_sem` | Erreur standard empirique locale. |
| `beta0_total_th` | \(\beta_0(t)\) théorique global. |
| `beta0_total_emp_true` | \(\beta_0(t)\) empirique global. |
| `N1_th`, `N2_th`, `N3_th` | Contributions théoriques des petites composantes. |
| `rmse_local`, `mae_local`, `bias_local` | Diagnostics théorie/empirique. |

### Fichier sauvegardé

```matlab
betti_t_phi_results.mat
```

---

---

# Sous-dossier `Paramètres/Nombre liens`

## `Paramètres/Nombre liens/strates_delta_spatiales.m`

### Objectif

Découpe la bande accessible du Walker-Delta en strates définies par la distance angulaire

\[
\beta=i-|\varphi|
\]

à la frontière orbitale.

La fonction calcule pour chaque strate :

- son aire réelle ;
- l’aire source associée ;
- le nombre moyen de satellites ;
- la densité locale ;
- le degré moyen local.

### Signature

```matlab
strates = strates_delta_spatiales( ...
    lambda, R, alpha_max, inc, varargin)
```

### Entrées principales

| Argument | Description |
|---|---|
| `lambda` | Densité surfacique initiale dans la bande. |
| `R` | Rayon orbital. |
| `alpha_max` | Angle maximal de liaison. |
| `inc` | Inclinaison commune. |
| `beta_step` | Largeur des strates, égale à `alpha_max/2` par défaut. |
| `beta_max` | Étendue maximale analysée. |
| `verbose` | Active l’affichage de la table. |

### Sorties principales

| Champ | Description |
|---|---|
| `all_table` | Table complète des strates. |
| `active_table` | Table des strates actives. |
| `A_link` | Aire d’une zone de liaison. |
| `beta_step`, `beta_max`, `beta_stop` | Paramètres du découpage. |

---

## `Paramètres/Nombre liens/liens_th.m`

### Objectif

Estime le nombre théorique de liens à partir des strates et du degré moyen local.

Pour une strate contenant en moyenne \(\mu_i\) satellites et de degré moyen \(k_i\),

\[
L_i^{\mathrm{pôle}}
\approx
\frac12\mu_i k_i.
\]

Le résultat est doublé pour prendre en compte les deux pôles.

### Signature

```matlab
[L_max_strates, details_table] = liens_th(strates_polaires)
```

### Sorties

| Sortie | Description |
|---|---|
| `L_max_strates` | Nombre total théorique de liens. |
| `details_table` | Contribution de chaque strate. |

---

## `Paramètres/Nombre liens/liens_emp.m`

### Objectif

Compte empiriquement, pour chaque strate, la contribution des liens incidents.

La contribution d’une strate est calculée par

\[
L_i
=
\frac12
\sum_{v\in\mathrm{strate}\ i}
\deg(v).
\]

Cette définition attribue la moitié d’un lien inter-strates à chacune des deux strates concernées.

### Signature

```matlab
[L_incident, n_sat] = ...
    liens_emp(A, z_t, R, strates_polaires)
```

### Entrées

| Argument | Description |
|---|---|
| `A` | Matrice d’adjacence. |
| `z_t` | Coordonnées verticales des satellites. |
| `R` | Rayon orbital. |
| `strates_polaires` | Structure des strates. |

### Sorties

| Sortie | Description |
|---|---|
| `L_incident` | Contribution en liens de chaque strate. |
| `n_sat` | Nombre de satellites par strate. |

---

## `Paramètres/Nombre liens/liens_sinus_temp.m`

### Objectif

Compare le nombre moyen empirique de liens à une sinusoïde théorique construite à partir d’un minimum et d’un maximum analytiques.

La forme utilisée est

\[
L_{\sin}(t)
=
\frac{L_{\max}+L_{\min}}{2}
-
\frac{L_{\max}-L_{\min}}{2}
\cos(2\omega t).
\]

La période du nombre de liens est donc la demi-période orbitale.

### Sorties principales

| Variable | Description |
|---|---|
| `L_min_theory`, `L_max_theory` | Minimum et maximum théoriques. |
| `L_sin_theory` | Sinusoïde théorique. |
| `mean_edges` | Nombre moyen empirique de liens. |
| `mean_links_internal_strate_t` | Liens internes moyens par strate. |
| `mean_n_sat_strate_t` | Satellites moyens par strate. |
| `rmse_sin`, `mae_sin` | Erreurs du modèle sinusoïdal. |
| `strates_comparison` | Comparaison détaillée par strate. |

### Figures produites

1. nombre de liens empirique et sinusoïde théorique ;
2. satellites par strate au pic ;
3. liens internes par strate au pic.

### Fichier sauvegardé

```matlab
liens_sinus_results.mat
```

---

## `Paramètres/Nombre liens/liens_quadrature_temp.m`

### Objectif

Calcule le nombre théorique temporel de liens par quadrature, à partir de la loi exacte de phase induite par l’uniformité spatiale initiale :

\[
f_u(u,t)
=
\frac{|\cos(u-\omega t)|}{4}.
\]

L’intégration sur la différence de RAAN est traitée analytiquement, puis l’intégrale en phase est évaluée par Gauss-Legendre.

### Sorties principales

| Variable | Description |
|---|---|
| `p_link_theory_time` | Probabilité théorique temporelle de lien. |
| `L_theory_time` | Nombre théorique temporel de liens. |
| `L_min_theory`, `L_max_theory` | Extrema théoriques. |
| `mean_edges` | Moyenne empirique. |
| `N_strate_theory_time` | Nombre théorique de satellites par strate. |
| `mean_links_internal_strate_t` | Liens internes empiriques par strate. |
| `mean_links_incident_strate_t` | Liens incidents empiriques par strate. |
| `comparison_strates` | Comparaison théorie/empirique au pic. |
| `rmse_links`, `mae_links` | Erreurs globales. |

### Figures produites

1. nombre de liens empirique et théorique ;
2. probabilité temporelle de lien ;
3. satellites par strate ;
4. liens internes par strate ;
5. liens incidents par strate.

### Fichier sauvegardé

```matlab
liens_quadrature_results.mat
```

---

## `Paramètres/Nombre liens/liens_t_phi.m`

### Objectif

Calcule le nombre local de liens

\[
E(t,\phi)
\]

à partir des résultats de `plink_t_phi.m`, sans recalculer la probabilité de lien.

Chaque lien est partagé entre les tranches de latitude de ses deux extrémités.

### Fichier d’entrée

```matlab
../plink_t_phi_results.mat
```

### Théorie

Pour une tranche \(b\),

\[
N_b^{\mathrm{th}}(t)
=
N
P(\Phi(t)\in b),
\]

puis

\[
E_b^{\mathrm{th}}(t)
=
\frac12
N_b^{\mathrm{th}}(t)
(N-1)
p_{\mathrm{link},b}^{\mathrm{th}}(t).
\]

### Mesure empirique

Pour chaque réalisation,

\[
E_b^{\mathrm{emp},(r)}(t)
=
\frac12
\sum_{i\in b}
\deg_i^{(r)}(t).
\]

Cette convention garantit

\[
\sum_b E_b(t)
=
E_{\mathrm{total}}(t).
\]

### Sorties principales

| Variable | Description |
|---|---|
| `edges_per_bin_th` | Nombre théorique de liens par tranche. |
| `edges_per_bin_emp_iterations` | Nombre empirique par réalisation. |
| `edges_per_bin_emp` | Moyenne empirique. |
| `edges_per_bin_emp_sem` | Erreur standard empirique. |
| `edges_density_th` | Densité théorique de liens par radian. |
| `edges_density_emp` | Densité empirique de liens par radian. |
| `edges_total_th` | Nombre total théorique de liens. |
| `edges_total_emp` | Nombre total empirique moyen. |
| `rmse_grid`, `mae_grid`, `bias_grid` | Diagnostics locaux. |

### Fichier sauvegardé

```matlab
liens_t_phi_results.mat
```

---

---

# Sous-dossier `Paramètres/Vitesse relative`

## `Paramètres/Vitesse relative/vrel_vs_vrad.m`

### Objectif

Compare la norme de la vitesse relative entre deux satellites à sa composante radiale suivant l’axe de leur séparation.

Pour une paire \((i,j)\),

\[
\mathbf v_{\mathrm{rel}}
=
\mathbf v_j-\mathbf v_i,
\qquad
v_{\mathrm{rad}}
=
\mathbf v_{\mathrm{rel}}\cdot\widehat{\mathbf r}_{ij}.
\]

La fonction peut utiliser :

- toutes les paires ;
- uniquement les liens existants ;
- les paires proches de la frontière \(d_{\max}\).

Pour la fusion, la composante utile est

\[
[-v_{\mathrm{rad}}]_+.
\]

### Signature

```matlab
[t, vrel_mean, vrad_mean, ratio_mean] = ...
    vrel_vs_vrad(pos_all, vel_all, t, dmax, varargin)
```

### Options principales

| Option | Description |
|---|---|
| `Pairs` | `all`, `linked` ou `near_boundary`. |
| `Width` | Largeur autour de \(d_{\max}\). |
| `UseAbsRadial` | Utilise \(|v_{\mathrm{rad}}|\) ou \([-v_{\mathrm{rad}}]_+\). |
| `MakeFigure` | Active les figures. |

### Sorties

| Sortie | Description |
|---|---|
| `vrel_mean` | Norme moyenne de la vitesse relative. |
| `vrad_mean` | Composante radiale moyenne utile. |
| `ratio_mean` | Rapport entre composante radiale et vitesse relative totale. |

### Figures produites

1. vitesse relative totale et composante radiale ;
2. rapport radial/total.

---

## `Paramètres/Vitesse relative/vrel_vs_vrad_emp.m`

### Objectif

Applique `vrel_vs_vrad.m` aux positions temporelles produites par `analysis_temp.m`.

Le script reconstruit les vitesses par différences finies, puis réalise trois diagnostics :

1. sur les liens existants ;
2. sur les paires proches de \(d_{\max}\) ;
3. sur la composante radiale entrante utile pour la fusion.

### Fichier d’entrée

```matlab
../../analysis_temp_results.mat
```

### Sorties principales

| Variable | Description |
|---|---|
| `vrel_linked` | Vitesse relative moyenne sur les liens existants. |
| `vrad_abs_linked` | Composante radiale absolue moyenne. |
| `vrel_bound` | Vitesse relative près de \(d_{\max}\). |
| `vrad_abs_bound` | Composante radiale absolue près de la frontière. |
| `vrad_merge` | Composante entrante utile pour la fusion. |
| `ratio_abs_bound` | Rapport radial absolu/total. |
| `ratio_merge` | Rapport radial entrant/total. |

### Figures produites

Le script produit les figures de diagnostic de `vrel_vs_vrad.m`, puis deux figures de synthèse près de \(d_{\max}\).

### Fichier sauvegardé

```matlab
vrel_vs_vrad_emp_results.mat
```

---

## `Paramètres/Vitesse relative/vitesse_rel_temp.m`

### Objectif

Compare la vitesse relative moyenne empirique des satellites liés à un modèle théorique temporel par quadrature et par strates.

La vitesse relative de deux satellites dont les vecteurs vitesse forment un angle \(\gamma\) est

\[
v_{\mathrm{rel}}
=
2v_{\mathrm{orb}}\sin\left(\frac{\gamma}{2}\right).
\]

Le modèle utilise la densité temporelle de phase

\[
f_u(u,t)
=
\frac{|\cos(u-\omega t)|}{4}.
\]

Il calcule, pour chaque strate, la probabilité qu’une paire soit liée et la moyenne conditionnelle de \(\sin(\gamma/2)\).

### Sorties principales

| Variable | Description |
|---|---|
| `vrel_theory` | Vitesse relative théorique temporelle. |
| `vrel_emp` | Vitesse relative empirique moyenne. |
| `gamma_eff_theory` | Angle effectif théorique. |
| `gamma_eff_emp` | Angle effectif empirique. |
| `L_strate_theory` | Nombre théorique de liens par strate. |
| `L_strate_emp` | Nombre empirique de liens par strate. |
| `s_strate_theory` | Moyenne de \(\sin(\gamma/2)\) par strate. |
| `rmse_vrel`, `mae_vrel` | Erreurs entre théorie et simulation. |

### Figures produites

1. vitesse relative empirique et théorique ;
2. angle effectif conditionné aux liens ;
3. liens par strate ;
4. facteur angulaire théorique par strate.

### Fichier sauvegardé

```matlab
vitesse_rel_temp.mat
```

---

## `Paramètres/Vitesse relative/vrel_t_phi.m`

### Objectif

Calcule et compare la vitesse relative locale conditionnée à l’existence d’un lien :

\[
v_{\mathrm{rel}}^{\mathrm{link}}(t,\phi)
=
\mathbb E
\left[
\|\mathbf v_1-\mathbf v_2\|
\mid
d_{12}\le d_{\max},
\Phi_1=\phi
\right].
\]

La latitude \(\phi\) est celle de la première extrémité. Un lien non orienté est donc pris en compte deux fois, une fois depuis chacune de ses extrémités, ce qui est cohérent avec la définition locale de \(p_{\mathrm{link}}(t,\phi)\).

### Fichier d’entrée

```matlab
../plink_t_phi_results.mat
```

### Théorie

Le calcul utilise deux noyaux géométriques :

- \(G(u_1,u_2)\), probabilité de lien après moyenne sur \(\Delta\Omega\) ;
- \(H(u_1,u_2)\), espérance de \(v_{\mathrm{rel}}\mathbf 1_{\{\mathrm{lien}\}}\).

Ainsi,

\[
v_{\mathrm{rel}}^{\mathrm{th}}(t,\phi)
=
\frac{
\mathbb E[
v_{\mathrm{rel}}\mathbf 1_{\{\mathrm{lien}\}}
\mid
\Phi_1=\phi]
}{
p_{\mathrm{link}}(t,\phi)
}.
\]

La moyenne sur une tranche est pondérée par le nombre attendu de liens orientés, donc par

\[
f_\Phi(t,\phi)
p_{\mathrm{link}}(t,\phi).
\]

### Sorties principales

| Variable | Description |
|---|---|
| `v_rel_th_fine` | Vitesse relative théorique sur la grille fine. |
| `v_rel_th_on_emp` | Théorie moyennée sur les tranches empiriques. |
| `v_rel_emp_iterations` | Vitesse empirique locale par réalisation. |
| `v_rel_emp` | Estimateur empirique agrégé. |
| `v_rel_emp_sem` | Erreur standard empirique. |
| `v_rel_th_global` | Vitesse relative globale théorique conditionnée au lien. |
| `v_rel_emp_global` | Vitesse relative globale empirique. |
| `rmse_grid`, `mae_grid`, `bias_grid` | Diagnostics locaux. |

### Fichier sauvegardé

```matlab
vrel_t_phi_results.mat
```

---

---

# Sous-dossier `Probabilité fusion`

## `Probabilité fusion/calc_p_merge_th.m`

### Objectif

Calcule la probabilité théorique temporelle de fusion :

\[
p_{\mathrm{merge}}(t)
=
1-\exp\!\left[
-2d_{\max}\lambda_{\mathrm{eff}}(t)
\chi_{\mathrm{merge}}(t)
v_{\mathrm{rel}}(t)\Delta t
\right].
\]

Le facteur topologique utilisé est

\[
\chi_{\mathrm{merge}}(t)
=
\frac{N-\beta_{0,\mathrm{geom}}}{E(t)},
\]

avec

\[
\beta_{0,\mathrm{geom}}
=
1+N_1+N_2+N_3.
\]

### Signature

```matlab
[p_merge_t, chi_merge_t, E_t, beta0_geom] = ...
    calc_p_merge_th(lambda_eff_t, vrel_t, dmax, dt, liens_file)
```

### Entrées

| Argument | Description |
|---|---|
| `lambda_eff_t` | Densité effective temporelle. |
| `vrel_t` | Vitesse relative temporelle. |
| `dmax` | Portée maximale. |
| `dt` | Pas temporel. |
| `liens_file` | Fichier contenant les données théoriques de lien. |

### Sorties

| Sortie | Description |
|---|---|
| `p_merge_t` | Probabilité théorique de fusion. |
| `chi_merge_t` | Facteur topologique temporel. |
| `E_t` | Nombre de liens utilisé. |
| `beta0_geom` | Approximation géométrique de \(\beta_0\). |

---

## `Probabilité fusion/pmerge_th.m`

### Objectif

Construit la série temporelle théorique \(p_{\mathrm{merge}}(t)\) à partir de :

- la densité effective ;
- la vitesse relative théorique ;
- le facteur topologique calculé par `calc_p_merge_th.m`.

### Fichiers d’entrée

```matlab
../Paramètres/lambda_eff_th_results.mat
../Paramètres/Vitesse relative/vitesse_rel_temp.mat
../Paramètres/Nombre liens/liens_quadrature_results.mat
```

Le dernier fichier est chargé indirectement par `calc_p_merge_th.m`.

### Sorties principales

| Variable | Description |
|---|---|
| `p_merge_t` | Série théorique temporelle. |
| `p_merge_vrel_emp_t` | Diagnostic utilisant la vitesse relative empirique. |
| `chi_merge_t` | Facteur topologique théorique. |
| `E_t_merge` | Nombre de liens utilisé. |
| `beta0_geom_merge` | Approximation de \(\beta_0\). |
| `p_merge_const` | Référence calculée avec des valeurs moyennes. |
| `chi_merge_const` | Facteur topologique de la référence constante. |

### Figure produite

Une figure compare :

- la théorie utilisant \(v_{\mathrm{rel}}^{\mathrm{th}}\) ;
- le diagnostic utilisant \(v_{\mathrm{rel}}^{\mathrm{emp}}\) ;
- la moyenne temporelle ;
- la référence constante.

### Fichier sauvegardé

```matlab
pmerge_th_results.mat
```

---

## `Probabilité fusion/pmerge_emp.m`

### Objectif

Calcule la probabilité empirique de fusion à partir du barcode zigzag \(H_0\).

Une fusion entre \(t_k\) et \(t_{k+1}\) correspond à une barre mourant à l’indice \(2k-1\).

La normalisation utilisée est

\[
p_{\mathrm{merge}}(k)
=
\frac{
\text{nombre de fusions}
}{
\beta_0(G_k)
}.
\]

### Fichiers d’entrée

```matlab
../barcodes_results.mat
../analysis_temp_results.mat
```

### Sorties principales

| Variable | Description |
|---|---|
| `merge_count` | Nombre de fusions par transition. |
| `merge_count_from_beta0` | Vérification par différence de \(\beta_0\). |
| `beta0_before` | Nombre de composantes avant la transition. |
| `beta0_union` | Nombre de composantes dans le graphe union. |
| `p_merge` | Probabilité empirique instantanée. |
| `p_merge_moving` | Moyenne glissante. |
| `p_merge_mean` | Moyenne globale pondérée. |
| `p_merge_time_mean` | Moyenne temporelle simple. |

### Figures produites

1. probabilité empirique de fusion ;
2. nombre de fusions par transition.

### Fichier sauvegardé

```matlab
pmerge_emp_results.mat
```

---

## `Probabilité fusion/pmerge_temp.m`

### Objectif

Compare temporellement les probabilités théorique et empirique de fusion.

La série théorique est interpolée sur la grille empirique avant le calcul des écarts.

### Fichiers d’entrée

```matlab
pmerge_th_results.mat
pmerge_emp_results.mat
```

### Sorties principales

| Variable | Description |
|---|---|
| `t_common` | Instants communs. |
| `p_th` | Probabilité théorique alignée. |
| `p_emp` | Probabilité empirique brute. |
| `p_emp_smooth` | Probabilité empirique lissée. |
| `p_merge_th_mean` | Moyenne théorique. |
| `p_merge_emp_mean` | Moyenne empirique pondérée. |
| `rmse`, `mae` | Indicateurs d’écart. |

### Figure produite

Une figure superpose les séries théorique et empirique ainsi que leurs moyennes.

### Fichier sauvegardé

```matlab
pmerge_temp_results.mat
```

---

## `Probabilité fusion/chi_temp.m`

### Objectif

Compare temporellement le facteur

\[
\chi_{\mathrm{merge}}(t)
=
\frac{N-\beta_0(t)}{E(t)}
\]

entre les données empiriques et plusieurs approximations théoriques.

Les trois modèles théoriques utilisent un \(\beta_0\) constant évalué à \(t=0\) avec :

1. les isolés ;
2. les isolés et dimères ;
3. les isolés, dimères et trimères.

### Fichiers d’entrée

```matlab
../analysis_temp_results.mat
../Paramètres/Nombre liens/liens_quadrature_results.mat
```

### Sorties principales

| Variable | Description |
|---|---|
| `chi_emp` | Facteur empirique. |
| `chi_th_iso` | Théorie avec isolés. |
| `chi_th_iso_dim` | Théorie avec isolés et dimères. |
| `chi_th_iso_dim_tri` | Théorie jusqu’aux trimères. |
| `beta0_iso` | Approximation constante avec isolés. |
| `beta0_iso_dim` | Approximation avec dimères. |
| `beta0_iso_dim_tri` | Approximation avec trimères. |

### Figure produite

Une figure compare les quatre courbes de \((N-\beta_0)/E\).

### Fichier sauvegardé

```matlab
liens_bridge_temp_results.mat
```

---

## `Probabilité fusion/pmerge_t_phi_th.m`

### Objectif

Calcule la probabilité théorique locale de fusion dépendant simultanément du temps et de la latitude :

\[
p_{\mathrm{merge}}(t,\phi)
=
1-
\exp\!\left[
-2d_{\max}\Delta T
\lambda(t,\phi)
N_C(t,\phi)
v_{\mathrm{rel}}(t,\phi)
\eta_{\mathrm{sweep}}(t,\phi)
\right],
\]

avec

\[
N_C(t,\phi)
=
\frac{
N(t,\phi)
}{
\beta_0(t,\phi)
}.
\]

### Fichiers d’entrée

```matlab
../Paramètres/N_t_phi_results.mat
../Paramètres/lambda_t_phi_results.mat
../Paramètres/eta_sweep_t_phi_results.mat
../Paramètres/Betti/betti_t_phi_results.mat
../Paramètres/Vitesse relative/vrel_t_phi_results.mat
```

### Moyenne globale

La probabilité globale est pondérée par le nombre local de composantes :

\[
p_{\mathrm{merge}}^{\mathrm{global}}(t)
=
\frac{
\sum_b
p_{\mathrm{merge}}(t,b)
\beta_0(t,b)
}{
\sum_b
\beta_0(t,b)
}.
\]

Sous l’hypothèse de fusions binaires,

\[
p_{\mathrm{disp,fusion}}(t,\phi)
=
\frac12
p_{\mathrm{merge}}(t,\phi).
\]

### Sorties principales

| Variable | Description |
|---|---|
| `mean_component_size_th` | Taille moyenne locale \(N/\beta_0\). |
| `mu_merge_th` | Exposant local du modèle de fusion. |
| `p_merge_th` | Probabilité locale théorique de fusion. |
| `p_disp_fusion_th` | Probabilité locale de disparition par fusion. |
| `p_merge_global_th` | Probabilité globale de fusion. |
| `p_disp_fusion_global_th` | Contribution globale à la disparition. |
| `p_merge_from_mean_mu` | Diagnostic obtenu en exponentiant l’exposant moyen. |
| `nonlinearity_gap` | Écart dû à la non-linéarité de l’exponentielle. |

### Fichier sauvegardé

```matlab
pmerge_t_phi_th_results.mat
```

---

---

# Sous-dossier `Probabilité rupture`

## `Probabilité rupture/calc_p_break_th.m`

### Objectif

Implémente la formule théorique linéaire de rupture du modèle Walker-Delta spatial :

\[
p_{\mathrm{break}}(t)
=
\frac{\mathbb E[|E(t)|]}
{\mathbb E[\beta_0(t)]}
\frac{2v_{\mathrm{rel}}(t)\Delta t}
{\pi d_{\max}}
\exp\!\left[
-\lambda_{\mathrm{eff}}(t)
A_{\mathrm{inter}}(d_{\max})
\right].
\]

L’aire d’intersection utilisée est

\[
A_{\mathrm{inter}}(d_{\max})
=
d_{\max}^2
\left(
\frac{2\pi}{3}-\frac{\sqrt3}{2}
\right).
\]

### Signature

```matlab
[p_break_t, details] = calc_p_break_th( ...
    vrel_t, lambda_eff_t, E_edges_t, beta0_t, dmax, dt)
```

### Entrées

| Argument | Description |
|---|---|
| `vrel_t` | Vitesse relative théorique. |
| `lambda_eff_t` | Densité effective. |
| `E_edges_t` | Nombre moyen de liens. |
| `beta0_t` | Nombre moyen de composantes. |
| `dmax` | Portée maximale. |
| `dt` | Pas temporel. |

### Sorties

| Sortie | Description |
|---|---|
| `p_break_t` | Probabilité théorique de rupture, bornée dans \([0,1]\). |
| `details` | Structure contenant les facteurs intermédiaires. |

---

## `Probabilité rupture/pbreak_th.m`

### Objectif

Calcule la série temporelle théorique \(p_{\mathrm{break}}(t)\) en combinant :

- la vitesse relative théorique ;
- la densité effective ;
- le nombre théorique de liens ;
- une approximation constante de \(\beta_0\).

Le calcul principal est délégué à `calc_p_break_th.m`.

### Fichiers d’entrée

```matlab
../Paramètres/Vitesse relative/vitesse_rel_temp.mat
../Paramètres/Nombre liens/liens_quadrature_results.mat
../Paramètres/lambda_eff_th_results.mat
```

### Approximation de \(\beta_0\)

Le script calcule d’abord

\[
\beta_{0,N123}(t)
=
1+N_1(t)+N_2(t)+N_3(t),
\]

puis conserve sa valeur initiale pendant toute la simulation :

\[
\beta_0(t)
=
\beta_{0,N123}(0).
\]

### Sorties principales

| Variable | Description |
|---|---|
| `p_break_t` | Série théorique temporelle. |
| `p_break_vrel_emp_t` | Diagnostic utilisant la vitesse empirique. |
| `p_break_const` | Référence construite avec les valeurs moyennes. |
| `beta0_initial` | Valeur constante de \(\beta_0\). |
| `E_edges_transition` | Nombre de liens sur les transitions. |
| `lambda_transition` | Densité effective sur les transitions. |
| `vrel_transition` | Vitesse relative théorique alignée. |
| `ratio_E_beta0_t` | Nombre moyen de liens par composante. |
| `exp_factor_t` | Facteur d’absence de voisin commun. |
| `q_break_link_t` | Probabilité de rupture d’un lien. |
| `p_break_raw_t` | Valeur avant saturation. |

### Figures produites

Le script compare notamment :

1. les variantes de \(p_{\mathrm{break}}(t)\) ;
2. la densité effective ;
3. les grandeurs topologiques ;
4. \(\beta_0\) variable et constant ;
5. les facteurs multiplicatifs ;
6. la valeur brute et la valeur saturée.

### Fichier sauvegardé

```matlab
pbreak_th_results.mat
```

---

## `Probabilité rupture/pbreak_emp.m`

### Objectif

Calcule la probabilité empirique de rupture depuis le barcode zigzag \(H_0\).

Une rupture entre \(t_k\) et \(t_{k+1}\) correspond à une barre naissant à l’indice \(2k+1\).

La normalisation utilisée est

\[
p_{\mathrm{break}}(k)
=
\frac{
\text{nombre de ruptures}
}{
\beta_0(G_k\cup G_{k+1})
}.
\]

Le dénominateur représente les composantes du graphe union exposées à une séparation.

### Fichiers d’entrée

```matlab
../barcodes_results.mat
../analysis_temp_results.mat
```

### Sorties principales

| Variable | Description |
|---|---|
| `break_count` | Nombre de ruptures par transition. |
| `break_count_from_beta0` | Vérification par différence de \(\beta_0\). |
| `beta0_union` | Nombre de composantes du graphe union. |
| `beta0_after` | Nombre de composantes après la transition. |
| `p_break` | Probabilité empirique instantanée. |
| `p_break_moving` | Moyenne glissante. |
| `p_break_mean` | Moyenne globale pondérée. |
| `p_break_time_mean` | Moyenne temporelle simple. |

### Figures produites

1. probabilité empirique de rupture ;
2. nombre de ruptures par transition.

### Fichier sauvegardé

```matlab
pbreak_emp_results.mat
```

---

## `Probabilité rupture/pbreak_temp.m`

### Objectif

Compare temporellement les probabilités théorique et empirique de rupture.

La série théorique est interpolée sur la grille temporelle empirique avant le calcul des écarts.

### Fichiers d’entrée

```matlab
pbreak_th_results.mat
pbreak_emp_results.mat
```

### Sorties principales

| Variable | Description |
|---|---|
| `t_common` | Instants communs aux deux séries. |
| `p_th` | Probabilité théorique alignée. |
| `p_emp` | Probabilité empirique brute. |
| `p_emp_smooth` | Probabilité empirique lissée. |
| `p_break_th_mean` | Moyenne théorique. |
| `p_break_emp_mean` | Moyenne empirique pondérée. |
| `rmse`, `mae` | Indicateurs d’écart. |

### Figure produite

Une figure superpose :

- la série empirique brute ;
- la série empirique lissée ;
- la série théorique ;
- leurs valeurs moyennes.

### Fichier sauvegardé

```matlab
pbreak_temp_results.mat
```

## `Probabilité rupture/pbreak_t_phi_th.m`

### Objectif

Calcule la probabilité théorique locale de rupture dépendant du temps et de la latitude.

La probabilité conditionnelle à une composante non isolée est

\[
p_{\mathrm{break}}^{\mathrm{noniso}}(t,\phi)
=
1-
\exp\!\left[
-\mu_{\mathrm{break}}(t,\phi)
\right],
\]

avec

\[
\mu_{\mathrm{break}}(t,\phi)
=
\frac{
E(t,\phi)
}{
\beta_0(t,\phi)-N_1(t,\phi)
}
p_{\mathrm{break}}^{\mathrm{lien}}(t,\phi)
p_{\mathrm{bridge,bord}}(t,\phi).
\]

### Fichiers d’entrée

```matlab
../Paramètres/Nombre liens/liens_t_phi_results.mat
../Paramètres/Betti/betti_t_phi_results.mat
../Paramètres/Betti/N1_t_phi_results.mat
../Paramètres/eta_sweep_t_phi_results.mat
../Paramètres/Vitesse relative/vrel_t_phi_results.mat
```

### Rupture d’un lien

La probabilité locale de rupture d’un lien est approchée par

\[
p_{\mathrm{break}}^{\mathrm{lien}}(t,\phi)
\simeq
\frac{2}{\pi}
\frac{
v_{\mathrm{rel}}^{\mathrm{link}}(t,\phi)
\Delta T
}{
d_{\max}
}.
\]

Le script utilise

\[
p_{\mathrm{bridge,bord}}(t,\phi)
=
\eta_{\mathrm{sweep}}(t,\phi).
\]

### Déconditionnement

La probabilité finale est

\[
p_{\mathrm{break}}(t,\phi)
=
\frac{
\beta_0(t,\phi)-N_1(t,\phi)
}{
\beta_0(t,\phi)
}
p_{\mathrm{break}}^{\mathrm{noniso}}(t,\phi).
\]

### Sorties principales

| Variable | Description |
|---|---|
| `p_break_link_th` | Probabilité locale de rupture d’un lien. |
| `nonisolated_component_count_th` | Nombre local de composantes non isolées. |
| `mean_edges_per_nonisolated_component_th` | Nombre moyen local de liens par composante non isolée. |
| `mu_break_nonisolated_th` | Exposant local de rupture. |
| `p_break_nonisolated_th` | Probabilité conditionnelle aux composantes non isolées. |
| `fraction_nonisolated_th` | Fraction locale de composantes non isolées. |
| `p_break_th` | Probabilité locale déconditionnée. |
| `p_break_linear_th` | Approximation linéaire. |
| `p_break_global_th` | Probabilité globale de rupture. |
| `p_break_nonisolated_global_th` | Probabilité globale conditionnelle. |
| `p_break_linear_global_th` | Approximation globale linéaire. |

### Fichier sauvegardé

```matlab
pbreak_t_phi_th_results.mat
```
