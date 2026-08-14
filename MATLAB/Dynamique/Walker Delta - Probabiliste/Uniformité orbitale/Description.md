# Description des scripts MATLAB

Ce document décrit les fichiers MATLAB, ainsi que leurs fonctions locales.

---

## `analysis_temp.m`

### Objectif

Simule un réseau LEO dynamique dans lequel chaque satellite possède :

- une inclinaison orbitale commune ;
- un RAAN `Omega` aléatoire ;
- une phase initiale `u0` aléatoire.

Le script construit le graphe des liens intersatellites à chaque instant, calcule les grandeurs topologiques principales, étudie la distribution des tailles de composantes, compare les nombres de composantes de tailles 1, 2 et 3 à des approximations théoriques, puis construit une suite zigzag par unions.

La suite utilisée est

\[
G_1 \to G_1\cup G_2 \leftarrow G_2 \to G_2\cup G_3 \leftarrow G_3 \to \cdots
\]

### Type

Script principal avec une fonction locale :

- `walker_delta_positions`

### Entrées du script

| Variable | Valeur | Description |
|---|---:|---|
| `R_earth` | `6371` km | Rayon terrestre. |
| `h` | `550` km | Altitude orbitale. |
| `mu` | `398600` km³/s² | Paramètre gravitationnel terrestre. |
| `lambda` | `4e-7` sat/km² | Intensité du processus de Poisson. |
| `inc_deg` | `58` degrés | Inclinaison orbitale commune. |
| `dmax` | `1500` km | Distance maximale d’un lien. |
| `dt` | `20` s | Pas temporel. |
| `Tmax` | `12000` s | Durée totale simulée. |

Le nombre de satellites suit

\[
N\sim\mathcal P\!\left(\lambda 4\pi R^2\right).
\]

### Modèle orbital

La phase évolue selon

\[
u(t)=u_0+\omega t,
\qquad
\omega=\sqrt{\mu/R^3}.
\]

Les positions sont données par

\[
x=R(\cos\Omega\cos u-\sin\Omega\sin u\cos i),
\]

\[
y=R(\sin\Omega\cos u+\cos\Omega\sin u\cos i),
\]

\[
z=R\sin u\sin i.
\]

### Sorties principales

| Variable | Description |
|---|---|
| `Positions` | Positions des satellites à chaque instant. |
| `Adjacency` | Matrices d’adjacence temporelles. |
| `beta0` | Nombre de composantes connexes. |
| `beta1_graph` | Nombre cyclomatique `E - N + beta0`. |
| `largest_component` | Taille de la plus grande composante. |
| `num_edges` | Nombre de liens. |
| `component_size_counts` | Nombre de composantes de chaque taille à chaque instant. |
| `mean_component_count_by_size` | Nombre moyen temporel de composantes de chaque taille. |
| `mean_component_fraction_by_size` | Distribution moyenne vue depuis une composante. |
| `mean_satellite_fraction_by_size` | Distribution moyenne vue depuis un satellite. |
| `mean_component_size` | Taille moyenne temporelle d’une composante. |
| `N1_emp_mean`, `N2_emp_mean`, `N3_emp_mean` | Moyennes empiriques des isolés, dimères et trimères. |
| `N1_theory`, `N2_theory`, `N3_theory` | Prédictions théoriques correspondantes. |
| `relative_error_N123` | Erreurs relatives associées. |
| `ZigzagAdjacency` | Matrices d’adjacence de la suite zigzag. |
| `ZigzagLabels` | Labels entiers et demi-entiers. |
| `beta0_zigzag` | Nombre de composantes sur le zigzag. |
| `beta1_zigzag_graph` | Nombre cyclomatique sur le zigzag. |
| `largest_component_zigzag` | Taille de la plus grande composante sur le zigzag. |
| `num_edges_zigzag` | Nombre de liens sur le zigzag. |

### Figures produites

Le script produit notamment :

1. `beta0(t)` ;
2. `beta1_graph(t)` ;
3. fraction de satellites dans la plus grande composante ;
4. nombre de liens ;
5. distributions des tailles de composantes ;
6. comparaison empirique/théorique de `N1`, `N2`, `N3` ;
7. métriques topologiques sur la suite zigzag.

### Fichier sauvegardé

```matlab
analysis_temp_results.mat
```

---

## `anim_3D.m`

### Objectif

Anime en trois dimensions le réseau LEO à inclinaison fixe et RAAN aléatoires.

Le script affiche :

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
| `h` | `550` km | Altitude orbitale. |
| `lambda` | `5e-7` sat/km² | Intensité du processus de Poisson. |
| `inc_deg` | `53` degrés | Inclinaison commune. |
| `dmax` | `1500` km | Portée maximale. |
| `dt` | `30` s | Pas temporel. |
| `Tmax` | `6000` s | Durée simulée. |

### Sorties

Le script ne retourne aucune valeur et ne sauvegarde aucun fichier.

Principales variables :

| Variable | Description |
|---|---|
| `positions0` | Positions initiales. |
| `positions_t` | Positions au dernier instant. |
| `Omega` | RAAN aléatoires. |
| `u0` | Phases initiales. |
| `A` | Dernière matrice d’adjacence. |
| `E` | Nombre de liens au dernier instant. |
| `sat_handle` | Objet graphique des satellites. |
| `link_handle` | Objet graphique des liens. |

---

## `barcodes.m`

### Objectif

Calcule le barcode zigzag en homologie \(H_0\) à partir des graphes produits par `analysis_temp.m`.

Le script :

- charge les graphes zigzag ;
- calcule les composantes de chaque objet ;
- construit les applications induites en \(H_0\) ;
- décompose le module zigzag en intervalles ;
- convertit les indices en temps physiques ;
- charge les probabilités théoriques de fusion et de rupture ;
- calcule `p_disp_th` et le temps caractéristique `tau_th` ;
- compare la survie empirique à un modèle exponentiel ;
- affiche le barcode et l’histogramme des durées.

### Type

Script principal avec huit fonctions locales :

- `build_H0_map`
- `zigzag_barcode_from_module_mod2`
- `filtration_quotient_dims`
- `gf2_preimage`
- `gf2_col_basis`
- `gf2_rank`
- `gf2_null`
- `gf2_rref`

### Fichiers d’entrée

```matlab
analysis_temp_results.mat
Probabilité fusion/pmerge_th_results.mat
Probabilité rupture/pbreak_th_results.mat
```

Pour la fusion, le script cherche `p_merge_mean`, puis `p_merge_delta`.

Pour la rupture, il utilise `p_break_delta`.

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
| `p_merge_th` | Probabilité théorique de fusion. |
| `p_break_th` | Probabilité théorique de rupture. |
| `p_disp_th` | Probabilité théorique de disparition. |
| `tau_th` | Temps caractéristique théorique. |
| `survival_th` | Courbe de survie théorique. |

### Figures produites

1. survie empirique et modèle exponentiel ;
2. barcode zigzag `H0` ;
3. histogramme des durées de vie.

### Fichier sauvegardé

```matlab
barcodes_results.mat
```

---

# Sous-dossier `Paramètres`

## `Paramètres/prob_routage_emp.m`

### Objectif

Calcule la probabilité temporelle de routage multi-sauts dans une constellation Walker-Delta à inclinaison fixe et RAAN aléatoires.

À chaque instant, la probabilité de routage est calculée à partir des tailles \(s_c\) des composantes connexes :

\[
P_{\mathrm{routing}}(t)
=
\frac{
\sum_c s_c(s_c-1)
}{
N(N-1)
}.
\]

Cette quantité est la probabilité que deux satellites distincts tirés uniformément appartiennent à la même composante.

### Type

Script principal avec une fonction locale :

- `walker_delta_positions`

### Entrées du script

| Variable | Valeur | Description |
|---|---:|---|
| `R_earth` | `6371` km | Rayon terrestre. |
| `h` | `550` km | Altitude. |
| `lambda` | `2e-7` sat/km² | Densité satellitaire. |
| `dmax` | `1500` km | Portée des liens. |
| `inc_deg` | `53` degrés | Inclinaison commune. |
| `dt` | `10` s | Pas temporel. |
| `Tmax` | `6000` s | Durée simulée. |

Le nombre de satellites suit une loi de Poisson sur la sphère orbitale.

### Sorties principales

| Variable | Description |
|---|---|
| `P_routing_time` | Probabilité de routage à chaque instant. |
| `num_edges_time` | Nombre de liens à chaque instant. |
| `largest_component_time` | Fraction de satellites dans la plus grande composante. |
| `time_values` | Instants simulés. |
| `N` | Nombre de satellites tiré. |
| `Omega`, `u0` | Paramètres orbitaux aléatoires. |

### Figures produites

1. probabilité de routage en fonction du temps ;
2. fraction de satellites dans la composante géante ;
3. nombre de liens temporel.

---

## `Paramètres/ponts_emp.m`

### Objectif

Compte les ponts exacts dans chaque composante connexe de chaque graphe temporel.

Un pont est une arête dont la suppression augmente le nombre de composantes connexes.

Le script calcule :

- le nombre de ponts par composante ;
- le nombre moyen de ponts par composante ;
- la fraction réelle de liens qui sont des ponts ;
- les moyennes conditionnelles aux composantes non triviales ;
- l’évolution temporelle de ces grandeurs.

### Type

Script principal avec une fonction locale :

- `find_bridges_tarjan`

### Fichier d’entrée

```matlab
analysis_temp_results.mat
```

Variables utilisées :

| Variable | Description |
|---|---|
| `Adjacency` | Matrices d’adjacence temporelles. |
| `time_values` | Instants de simulation. |
| `N` | Nombre de satellites. |
| `dmax` | Portée des liens, si disponible. |
| `dt` | Pas temporel, si disponible. |

### Sorties principales

| Variable | Description |
|---|---|
| `component_table` | Tableau de toutes les composantes observées. |
| `n_components_t` | Nombre de composantes à chaque instant. |
| `n_nontrivial_components_t` | Nombre de composantes de taille au moins 2. |
| `n_components_with_bridge_t` | Nombre de composantes possédant au moins un pont. |
| `n_edges_t` | Nombre de liens. |
| `n_bridges_t` | Nombre de ponts. |
| `mean_bridges_per_component_t` | Ponts moyens par composante. |
| `mean_bridges_per_nontrivial_component_t` | Ponts moyens par composante non triviale. |
| `bridge_fraction_t` | Fraction de liens ponts. |
| `mean_edges_per_component_t` | Liens moyens par composante. |

Les agrégats globaux incluent notamment

\[
\overline B
=
\frac{\text{nombre total de ponts}}
{\text{nombre total de composantes}},
\]

et

\[
\chi_{\mathrm{bridge}}^{\mathrm{emp}}
=
\frac{\text{nombre total de ponts}}
{\text{nombre total de liens}}.
\]

### Figures produites

1. nombre moyen de ponts par composante ;
2. fraction réelle de liens critiques ;
3. distribution du nombre de ponts par composante ;
4. nombre de ponts en fonction de la taille des composantes.

### Fichier sauvegardé

```matlab
ponts_emp_results.mat
```

---

## `Paramètres/plink_th.m`

### Objectif

Calcule la probabilité théorique globale de lien dans le modèle Walker-Delta à inclinaison fixe.

Les paramètres orbitaux sont

\[
\Omega\sim\mathcal U(0,2\pi),
\qquad
u\sim\mathcal U(0,2\pi).
\]

L’intégration sur la différence de RAAN \(\Delta\Omega\) est résolue analytiquement. Il reste une quadrature déterministe bidimensionnelle sur \((u_1,u_2)\).

### Type

Script principal avec une fonction locale :

- `gauss_legendre_interval`

### Entrées du script

| Variable | Valeur | Description |
|---|---:|---|
| `R_earth` | `6371` km | Rayon terrestre. |
| `h` | `550` km | Altitude. |
| `inc_deg` | `90` degrés | Inclinaison. |
| `dmax` | `1500` km | Distance maximale. |
| `quad_orders` | `[50 100 150 200 300 400]` | Ordres de quadrature testés. |

### Sorties principales

| Variable | Description |
|---|---|
| `p_quad` | Valeurs obtenues pour chaque ordre de quadrature. |
| `p_link_delta` | Valeur finale Walker-Delta. |
| `p_link_sphere` | Référence uniforme sur la sphère. |
| `E_delta` | Nombre moyen d’arêtes pour le modèle Delta. |
| `k_delta` | Degré moyen Delta. |
| `E_sphere` | Nombre moyen d’arêtes pour le modèle uniforme sphérique. |
| `k_sphere` | Degré moyen sphérique. |

### Figure produite

Une figure de convergence de la quadrature en fonction de l’ordre.

### Fichier sauvegardé

```matlab
plink_results.mat
```

---

## `Paramètres/distribution_latitude.m`

### Objectif

Vérifie numériquement la distribution spatiale induite par des RAAN et phases orbitales uniformes à inclinaison fixe.

Le script compare :

- la distribution empirique de la latitude ;
- la densité théorique induite par \(u\) uniforme ;
- une distribution uniforme sur la bande accessible.

La latitude vérifie

\[
\sin\varphi
=
\sin i\,\sin u.
\]

### Type

Script principal avec une fonction locale :

- `walker_delta_positions`

### Entrées du script

| Variable | Valeur | Description |
|---|---:|---|
| `R` | `1` | Rayon normalisé. |
| `inc_deg` | `53` degrés | Inclinaison. |
| `M` | `1e6` | Nombre de positions simulées. |
| `nbins` | `120` | Nombre de classes d’histogramme. |

### Densités comparées

La densité théorique de latitude est

\[
f_\Phi(\varphi)
=
\frac{\cos\varphi}
{\pi\sqrt{\sin^2 i-\sin^2\varphi}}.
\]

La densité uniforme sur la bande en surface serait

\[
f_{\mathrm{bande}}(\varphi)
=
\frac{\cos\varphi}{2\sin i}.
\]

Pour \(s=\sin\varphi\), la densité induite est

\[
f_S(s)
=
\frac{1}
{\pi\sqrt{\sin^2i-s^2}}.
\]

### Sorties principales

| Variable | Description |
|---|---|
| `latitude` | Latitudes simulées. |
| `longitude` | Longitudes simulées. |
| `sin_latitude` | Valeurs de \(\sin\varphi\). |
| `counts_lat` | Densité empirique de latitude. |
| `counts_s` | Densité empirique de \(\sin\varphi\). |
| `counts_lon` | Densité empirique de longitude. |
| `h_band`, `p_band` | Résultat du test de Kolmogorov-Smirnov contre l’uniformité de bande. |

### Figures produites

1. distribution des latitudes ;
2. distribution de \(\sin(\text{latitude})\) ;
3. distribution des longitudes ;
4. nuage 3D de positions.

---

## `Paramètres/densite_phi.m`

### Objectif

Calcule et valide la densité surfacique locale Walker-Delta en fonction de la latitude.

Sous les hypothèses

\[
\Omega\sim\mathcal U(0,2\pi),
\qquad
u\sim\mathcal U(0,2\pi),
\]

la latitude satisfait

\[
\sin\varphi=\sin i\,\sin u.
\]

### Entrées du script

| Variable | Valeur | Description |
|---|---:|---|
| `R_earth` | `6371` km | Rayon terrestre. |
| `h` | `550` km | Altitude. |
| `inc_deg` | `90` degrés | Inclinaison. |
| `N` | `204` | Nombre de satellites. |
| `n_realizations` | `10000` | Nombre de constellations simulées. |
| `n_bins` | `80` | Nombre de bandes de latitude. |

### Sorties principales

| Variable | Description |
|---|---|
| `phi` | Grille continue de latitude. |
| `f_phi_theory` | Densité théorique de latitude. |
| `lambda_theory` | Densité locale théorique continue. |
| `lambda_theory_bins` | Densité théorique moyenne par bande. |
| `lambda_empirical` | Densité empirique moyenne par bande. |
| `lambda_band_mean` | Densité moyenne sur la bande. |
| `lambda_sphere_mean` | Densité moyenne sphérique. |
| `rmse` | Erreur quadratique moyenne théorie/simulation. |
| `relative_l1_error` | Erreur \(L^1\) relative pondérée par les aires. |

### Figures produites

1. densité locale théorique et empirique ;
2. facteur de concentration locale ;
3. densité de probabilité de latitude.

### Fichier sauvegardé

```matlab
densite_phi_results.mat
```

---

## `Paramètres/betti_th.m`

### Objectif

Estime de manière améliorée les nombres moyens de composantes isolées de tailles 1, 2 et 3 dans un modèle Walker-Delta à uniformité orbitale.

Le script calcule :

- \(N_1\) par quadrature déterministe locale ;
- \(N_2\) par quadrature imbriquée 3D + 2D ;
- \(N_3\) par quadrature imbriquée 5D + 2D ;
- une approximation de \(\beta_0\) ;
- une validation par simulations de graphes statiques complets.

### Type

Script principal avec trois fonctions locales :

- `walker_delta_position`
- `sample_walker_delta_unit`
- `gauss_legendre_interval`

### Entrées du script

| Variable | Valeur | Description |
|---|---:|---|
| `R` | `6371 + 550` km | Rayon orbital. |
| `inc_deg` | `90` degrés | Inclinaison. |
| `dmax` | `1500` km | Portée de connexion. |
| `N` | `204` | Nombre de satellites. |
| `n_graph_realizations` | `300` | Nombre de graphes simulés pour la validation. |

### Sorties principales

| Variable | Description |
|---|---|
| `phi_vals` | Latitudes utilisées pour l’étude locale. |
| `p_link_phi` | Probabilité locale de lien. |
| `N1_th_local` | Estimation théorique des isolés. |
| `N2_th_geom` | Estimation théorique des dimères. |
| `N3_th_geom` | Estimation théorique des trimères. |
| `beta0_th_123` | Approximation de \(\beta_0\). |
| `N1_emp`, `N2_emp`, `N3_emp` | Valeurs empiriques sur les graphes simulés. |
| `beta0_emp` | Nombre de composantes simulé. |

### Figures produites

1. comparaison théorie/simulation de \(N_1\), \(N_2\), \(N_3\) ;
2. distribution empirique de \(\beta_0\) et approximation ;
3. probabilité locale de lien en fonction de la latitude.

### Fichier sauvegardé

```matlab
betti_results.mat
```

---

## `Paramètres/vrel_vrad_temp.m`

### Objectif

Étudie la **vitesse relative et la vitesse radiale sortante des liens intersatellites en fonction de la latitude** dans le modèle Walker Delta à uniformité orbitale.

Les grandeurs sont calculées sur les liens existants puis moyennées temporellement par tranche de latitude. Pour un lien \((i,j)\),

\[
v_{\mathrm{rel},ij}
=
\|\mathbf v_j-\mathbf v_i\|,
\]

et

\[
v_{\mathrm{rad},ij}
=
(\mathbf v_j-\mathbf v_i)\cdot
\frac{\mathbf r_j-\mathbf r_i}
{\|\mathbf r_j-\mathbf r_i\|}.
\]

La composante pertinente pour les ruptures est

\[
v_{\mathrm{rad,out}}
=
\max(v_{\mathrm{rad}},0).
\]

Le script distingue également les liens appartenant à la couche de rupture définie par

\[
d_{ij}+v_{\mathrm{rad,out}}\,dt\ge d_{\max}.
\]

La latitude associée à un lien est celle de son milieu géométrique reprojeté sur la sphère.

### Type

Script principal avec fonctions locales de calcul et de traitement des moyennes.

### Fichier d’entrée

```matlab
analysis_temp_results.mat
```

Le fichier doit notamment contenir :

```matlab
Positions, Adjacency, dt, dmax
```

Le script récupère également `R`, `mu`, `inc` ou `inc_deg` lorsqu’ils sont disponibles.

### Modèle théorique local

Pour les liens proches de la rupture, le modèle utilisé est

\[
v_{\mathrm{rel}}^{\Delta}(\phi)
=
2v_{\mathrm{orb}}
\frac{\sqrt{\sin^2 i-\sin^2\phi}}
{\cos\phi},
\]

avec

\[
v_{\mathrm{orb}}
=
\sqrt{\frac{\mu}{R}}.
\]

La projection radiale sortante conditionnée par la couche de rupture est approchée par

\[
v_{\mathrm{rad,out}}^{\Delta}(\phi)
=
\frac{\pi}{4}
v_{\mathrm{rel}}^{\Delta}(\phi).
\]

### Sorties principales

| Variable | Description |
|---|---|
| `vorb_time_mean_lat` | Vitesse orbitale moyenne temporelle par latitude. |
| `vrel_time_mean_lat` | Vitesse relative moyenne des liens par latitude. |
| `vrad_out_time_mean_lat` | Composante radiale sortante moyenne par latitude. |
| `vrel_border_time_mean_lat` | Vitesse relative moyenne des liens proches de la rupture. |
| `vrad_out_border_time_mean_lat` | Vitesse radiale sortante moyenne des liens proches de la rupture. |
| `vrel_border_theory_lat` | Modèle théorique local de vitesse relative au bord. |
| `vrad_out_border_theory_lat` | Modèle théorique local de vitesse radiale sortante au bord. |
| `ratio_vrel_emp_theory_lat` | Rapport empirique/théorie pour \(v_{\mathrm{rel}}\). |
| `ratio_vrad_emp_theory_lat` | Rapport empirique/théorie pour \(v_{\mathrm{rad,out}}\). |
| `results_by_latitude` | Tableau récapitulatif par tranche de latitude. |

Le script calcule à la fois des moyennes temporelles, où chaque instant non vide possède le même poids, et des moyennes agrégées, où chaque observation individuelle possède le même poids.

### Fichier sauvegardé

```matlab
vrel_vrad_emp_delta_latitude_results.mat
```

---

## `Paramètres/eta_sweep_phi.m`

### Objectif

Calcule et valide le **facteur local de redondance spatiale**

\[
\eta_{\mathrm{sweep}}^\Delta(\phi),
\]

qui représente la fraction de la zone nouvellement balayée par un satellite qui n’était pas déjà couverte par les autres satellites de sa composante.

Le modèle théorique utilise directement la densité surfacique locale :

\[
\eta_{\mathrm{sweep}}^{\mathrm{th}}(\phi)
=
\exp\!\left[
-\lambda(\phi)A_{\mathrm{intersection}}
\right],
\]

avec

\[
A_{\mathrm{intersection}}
=
\left(
\frac{2\pi}{3}-\frac{\sqrt 3}{2}
\right)d_{\max}^2.
\]

Le script compare ce modèle à une estimation construite avec la densité empirique ainsi qu’à une mesure directe de la zone réellement nouvelle.

### Type

Script principal avec plusieurs fonctions locales, notamment pour :

- rechercher les fichiers d’entrée ;
- échantillonner une calotte sphérique ;
- détecter les ponts exacts du graphe.

### Fichiers d’entrée

```matlab
densite_phi_results.mat
analysis_temp_results.mat
```

### Sorties principales

| Variable | Description |
|---|---|
| `phi_vals` | Centres des tranches de latitude. |
| `lambda_phi_th` | Densité surfacique locale théorique. |
| `lambda_phi_emp_density` | Densité surfacique locale empirique. |
| `eta_sweep_phi_th` | Modèle théorique de \(\eta_{\mathrm{sweep}}(\phi)\). |
| `eta_sweep_phi_from_empirical_density` | Modèle utilisant la densité locale empirique. |
| `eta_sweep_phi_emp_direct` | Mesure directe de \(\eta_{\mathrm{sweep}}(\phi)\). |
| `eta_sweep_phi_emp_mean_per_sat` | Moyenne locale donnant le même poids à chaque satellite testé. |
| `p_bridge_bord_phi_emp_true` | Probabilité empirique \(P(\text{pont}\mid\text{rupture},\phi)\). |
| `p_no_common_neighbor_given_break_phi` | Probabilité de ne pas avoir de voisin commun conditionnellement à une rupture. |
| `ratio_direct_th_phi` | Rapport mesure directe/modèle théorique. |
| `rmse_direct_vs_th` | RMSE entre mesure directe et théorie. |

Le script mesure donc également la vraie probabilité locale

\[
p_{\mathrm{bridge,bord}}^{\mathrm{emp}}(\phi)
=
P(\text{pont à }t
\mid
\text{lien rompu entre }t\text{ et }t+dt,\phi).
\]

### Fichier sauvegardé

```matlab
eta_sweep_phi_results.mat
```

---

## `Paramètres/loi_distances_phi.m`

### Objectif

Vérifie empiriquement la **loi des distances des liens**, globalement et en fonction de la latitude, puis la compare au modèle local uniforme

\[
f_D(d)
=
\frac{2d}{d_{\max}^2},
\qquad
0\le d\le d_{\max},
\]

et

\[
F_D(d)
=
\left(\frac{d}{d_{\max}}\right)^2.
\]

Le script étudie en particulier la densité de probabilité au voisinage de \(d_{\max}\), quantité directement impliquée dans le calcul des ruptures de liens.

### Type

Script principal avec une fonction locale :

- `safe_divide`

### Fichier d’entrée

```matlab
analysis_temp_results.mat
```

Variables nécessaires :

```matlab
Positions, Adjacency, dmax
```

### Calcul local

Les liens sont répartis selon la latitude du milieu géométrique normalisé de leurs deux extrémités.

Pour une couronne de largeur \(\Delta d\) proche de \(d_{\max}\), le script estime

\[
P(D\ge d_{\max}-\Delta d\mid \text{lien},\phi),
\]

puis

\[
f_{D\mid\text{lien},\phi}(d_{\max})
\approx
\frac{
P(D\ge d_{\max}-\Delta d\mid \text{lien},\phi)
}{
\Delta d
}.
\]

Le modèle uniforme donne

\[
f_D(d_{\max})
=
\frac{2}{d_{\max}}.
\]

### Sorties principales

| Variable | Description |
|---|---|
| `all_link_distances` | Ensemble des distances de liens extraites. |
| `all_link_latitudes` | Latitude associée à chaque lien. |
| `pdf_emp_global` | PDF empirique globale des distances. |
| `cdf_emp_global` | CDF empirique globale. |
| `pdf_uniform_local` | PDF théorique \(2d/d_{\max}^2\). |
| `pdf_emp_by_latitude` | PDF empirique des distances pour chaque tranche de latitude. |
| `boundary_probability_emp_by_latitude` | Probabilité locale d’être proche de \(d_{\max}\). |
| `boundary_density_emp_by_latitude` | Densité locale estimée au bord. |
| `boundary_density_ratio_emp_theory` | Rapport entre densité empirique et modèle \(2/d_{\max}\). |
| `ks_distance` | Écart maximal entre CDF empirique et CDF théorique. |
| `rmse_pdf_global` | RMSE de la PDF globale. |

### Figures produites

Le script produit notamment :

1. la PDF globale empirique et théorique ;
2. la CDF globale ;
3. la PDF des distances selon la latitude ;
4. la densité locale au voisinage de \(d_{\max}\) ;
5. la probabilité d’appartenir à la couronne de bord ;
6. le nombre de liens échantillonnés par latitude.

### Fichier sauvegardé

```matlab
verif_loi_distances_liens_results.mat
```

---

# Sous-dossier `Probabilité fusion`

## `Probabilité fusion/chi_temp.m`

### Objectif

Calcule et compare le facteur correctif topologique

\[
\chi_{\mathrm{merge}}
=
\frac{N-\beta_0}{E_\Delta}
\]

pour le modèle Walker-Delta à uniformité orbitale.

La quantité \(N-\beta_0\) correspond au nombre d’arêtes d’une forêt couvrante des composantes. Elle représente donc le nombre minimal de liens topologiquement indispensables pour relier les sommets à l’intérieur des composantes.

Le nombre moyen total de liens est

\[
E_\Delta
=
\binom N2p_{\mathrm{link}}^\Delta.
\]

Le script compare :

- une distribution empirique obtenue à partir des réalisations de \(\beta_0\) ;
- la moyenne empirique ;
- l’approximation théorique de \(\beta_0\) incluant les composantes jusqu’aux trimères.

### Type

Script principal avec deux fonctions locales :

- `get_field_or_nan`
- `get_optional_scalar`

### Fichiers d’entrée

Le script cherche un niveau au-dessus :

```matlab
../plink_results.mat
../betti_results.mat
```

### Variables nécessaires

Depuis `plink_results.mat` :

| Variable | Description |
|---|---|
| `N` | Nombre de satellites. |
| `p_link_delta` | Probabilité moyenne de lien du modèle Delta. |
| `E_delta` | Nombre moyen de liens, lorsqu’il est sauvegardé. |

Depuis `betti_results.mat` :

| Variable | Description |
|---|---|
| `N` | Nombre de satellites. |
| `beta0_emp` | Réalisations empiriques de \(\beta_0\). |
| `beta0_th_123` | Approximation théorique jusqu’aux trimères. |

Le script vérifie que les deux fichiers utilisent la même valeur de \(N\).

### Sorties principales

| Variable | Description |
|---|---|
| `E_delta` | Nombre moyen théorique de liens. |
| `beta0_emp_mean` | Moyenne empirique de \(\beta_0\). |
| `beta0_emp_std` | Écart-type empirique. |
| `beta0_th_123` | Approximation théorique de \(\beta_0\). |
| `n_bridge_equiv_emp` | Liens indispensables équivalents pour chaque réalisation. |
| `n_bridge_equiv_th` | Valeur théorique correspondante. |
| `chi_bridge_emp` | Valeurs empiriques du facteur correctif. |
| `chi_bridge_emp_mean` | Moyenne empirique du facteur. |
| `chi_bridge_emp_std` | Écart-type du facteur. |
| `chi_bridge_th_123` | Facteur théorique jusqu’aux trimères. |
| `T` | Tableau récapitulatif. |

### Figures produites

1. distribution empirique de \(\beta_0\) ;
2. distribution empirique de \(\chi_{\mathrm{merge}}\) ;
3. comparaison empirique/théorique de \(\beta_0\) ;
4. comparaison empirique/théorique du facteur correctif ;
5. profil local de \(p_{\mathrm{link}}(\phi)\), lorsqu’il est disponible.

---

## `Probabilité fusion/pmerge_emp.m`

### Objectif

Calcule la probabilité empirique de fusion \(p_{\mathrm{merge}}\) à partir des résultats du barcode zigzag en homologie \(H_0\).

Pour une transition

\[
G_k
\longrightarrow
G_k\cup G_{k+1}
\longleftarrow
G_{k+1},
\]

une fusion correspond à une barre qui meurt à l’indice \(2k-1\), lors du passage de \(G_k\) vers le graphe union.

La probabilité est normalisée par le nombre de composantes exposées dans \(G_k\) :

\[
p_{\mathrm{merge}}(k)
=
\frac{\text{nombre de fusions à la transition }k}
{\beta_0(G_k)}.
\]

### Fichiers d’entrée

Le script cherche d’abord :

```matlab
../barcodes_results.mat
../analysis_temp_results.mat
```

Il teste ensuite leur présence dans son propre dossier pour assurer une compatibilité avec une autre organisation des fichiers.

Une compatibilité est également prévue avec l’ancien fichier :

```matlab
leo_H0_zigzag_barcodes_delta.mat
```

### Variables chargées

Depuis le fichier de barcode :

| Variable | Description |
|---|---|
| `death_index` | Indices de mort des barres. |
| `ZigzagTime` | Temps associés aux objets zigzag. |
| `h0_dims` | Dimensions de \(H_0\), donc nombres de composantes. |

Depuis le fichier d’analyse, lorsqu’il existe :

| Variable | Description |
|---|---|
| `time_values` | Instants des graphes réels. |
| `dt` | Pas temporel. |
| `inc_deg` | Inclinaison. |
| `N` | Nombre de satellites. |

### Sorties principales

| Variable | Description |
|---|---|
| `merge_count` | Nombre de fusions à chaque transition. |
| `merge_count_from_beta0` | Comptage obtenu par différence de \(\beta_0\). |
| `beta0_before` | Nombre de composantes avant la fusion. |
| `beta0_union` | Nombre de composantes dans le graphe union. |
| `p_merge` | Probabilité empirique instantanée. |
| `p_merge_moving` | Moyenne glissante. |
| `merge_count_moving` | Moyenne glissante du nombre de fusions. |
| `p_merge_mean` | Moyenne globale pondérée. |
| `p_merge_time_mean` | Moyenne temporelle simple. |
| `max_merge_error` | Écart maximal entre les deux méthodes de comptage. |
| `time_transition` | Temps au milieu des transitions. |

La moyenne globale pondérée est

\[
\overline p_{\mathrm{merge}}
=
\frac{\sum_k\text{merge\_count}(k)}
{\sum_k\beta_0(G_k)}.
\]

### Figures produites

1. probabilité empirique instantanée, moyenne glissante et moyenne globale ;
2. nombre de fusions entre graphes consécutifs.

### Fichier sauvegardé

```matlab
pmerge_emp_results.mat
```

Le fichier est sauvegardé dans le même dossier que le script.

---

## `Probabilité fusion/pmerge_temp.m`

### Objectif

Compare l’évolution temporelle de la probabilité empirique de fusion à la moyenne théorique calculée par `pmerge_th.m`.

Le script superpose :

- \(p_{\mathrm{merge}}^{\mathrm{emp}}(t)\) ;
- sa moyenne glissante ;
- la moyenne empirique globale ;
- la moyenne orbitale théorique.

### Fichiers d’entrée

```matlab
pmerge_emp_results.mat
pmerge_th_results.mat
```

### Variables chargées

Depuis le fichier empirique :

| Variable | Description |
|---|---|
| `time_transition` | Temps associés aux transitions. |
| `p_merge` | Probabilité empirique instantanée. |
| `p_merge_moving` | Moyenne glissante empirique. |
| `p_merge_mean` | Moyenne empirique globale pondérée. |
| `p_merge_time_mean` | Moyenne temporelle simple. |
| `moving_window` | Taille de la fenêtre glissante. |

Depuis le fichier théorique :

| Variable | Description |
|---|---|
| `p_merge_mean` | Moyenne orbitale théorique de fusion. |

La variable théorique est renommée :

```matlab
p_merge_th = Sth.p_merge_mean;
```

### Sorties

| Variable | Description |
|---|---|
| `p_merge_emp` | Courbe empirique instantanée. |
| `p_merge_emp_moving` | Moyenne glissante empirique. |
| `p_merge_emp_mean` | Moyenne empirique globale. |
| `p_merge_emp_time_mean` | Moyenne temporelle simple. |
| `p_merge_th` | Moyenne orbitale théorique. |

### Figure produite

Une figure comparant \(p_{\mathrm{merge}}^{\mathrm{emp}}(t)\) et \(p_{\mathrm{merge}}^{\mathrm{th}}\).

### Affichage console

Le script affiche :

- la valeur théorique ;
- la moyenne empirique globale ;
- la moyenne temporelle simple ;
- l’écart absolu ;
- le rapport théorie/empirique.

---

## `Probabilité fusion/pmerge_phi_th.m`

### Objectif

Calcule la **probabilité théorique locale de fusion**

\[
p_{\mathrm{merge}}^\Delta(\phi)
\]

pour le modèle Walker Delta à uniformité orbitale, puis reconstruit la probabilité globale par une moyenne pondérée par la distribution des composantes connexes.

La formule locale utilisée est

\[
p_{\mathrm{merge}}^\Delta(\phi)
=
1-\exp\!\left[
-4d_{\max}v_{\mathrm{orb}}dt\,
\lambda(\phi)\,
\frac{N^\Delta(\phi)}
{B_0^\Delta(\phi)}
\frac{\sqrt{\sin^2i-\sin^2\phi}}
{\cos\phi}
\eta_{\mathrm{sweep}}(\phi)
\right].
\]

Les termes locaux sont :

- \(N^\Delta(\phi)\) : densité de satellites par radian de latitude ;
- \(B_0^\Delta(\phi)\) : densité locale de composantes ;
- \(\lambda(\phi)\) : densité surfacique locale ;
- \(\eta_{\mathrm{sweep}}(\phi)\) : facteur local de redondance spatiale.

### Type

Script principal avec plusieurs fonctions locales de chargement, de validation et de reconstruction empirique de \(\beta_0(\phi)\).

### Fichiers d’entrée

```matlab
N_phi_results.mat
betti_phi_results.mat
densite_phi_results.mat
eta_sweep_phi_results.mat
analysis_temp_results.mat
```

Le dernier fichier est notamment utilisé pour récupérer `dt`, `mu` et reconstruire un \(\beta_0(\phi)\) empirique.

### Vitesse relative locale

Le modèle utilise

\[
v_{\mathrm{rel}}^\Delta(\phi)
=
2v_{\mathrm{orb}}
\frac{\sqrt{\sin^2i-\sin^2\phi}}
{\cos\phi},
\qquad
v_{\mathrm{orb}}
=
\sqrt{\frac{\mu}{R}}.
\]

### Moyenne globale

La probabilité globale est pondérée par la loi de latitude des composantes :

\[
p_{\mathrm{merge}}^\Delta
=
\frac{
\int
p_{\mathrm{merge}}^\Delta(\phi)
B_0^\Delta(\phi)\,d\phi
}{
\int
B_0^\Delta(\phi)\,d\phi
}.
\]

Numériquement, chaque tranche est donc pondérée par le nombre théorique de composantes fourni par `betti_phi_results.mat`.

Le script conserve aussi, pour comparaison, une moyenne pondérée par la loi des satellites.

### Version corrigée

Une seconde courbe utilise :

- le \(\beta_0(\phi)\) empirique reconstruit depuis `analysis_temp_results.mat` ;
- la mesure empirique directe de \(\eta_{\mathrm{sweep}}(\phi)\).

Elle fournit une version corrigée permettant d’identifier les écarts dus au modèle topologique et au facteur de redondance spatiale.

### Sorties principales

| Variable | Description |
|---|---|
| `phi_vals` | Grille de latitude. |
| `satellites_density_phi` | \(N^\Delta(\phi)\). |
| `betti0_density_phi` | \(B_0^\Delta(\phi)\) théorique. |
| `lambda_phi` | Densité surfacique locale utilisée. |
| `eta_sweep_phi` | Facteur local de redondance spatiale. |
| `geometry_factor` | Facteur \(\sqrt{\sin^2i-\sin^2\phi}/\cos\phi\). |
| `v_rel_phi` | Vitesse relative locale. |
| `merge_exponent_phi` | Exposant local du modèle de fusion. |
| `p_merge_phi_th` | Probabilité théorique locale de fusion. |
| `component_probability_bin` | Loi discrète de latitude des composantes. |
| `p_merge_th` | Probabilité théorique globale de fusion. |
| `p_disp_fusion_th` | Contribution de la fusion à la disparition d’une barre \(H_0\). |
| `p_merge_phi_th_corrected` | Probabilité locale avec corrections empiriques. |
| `p_merge_th_corrected` | Probabilité globale corrigée. |
| `correction_betti_phi` | Facteur correctif associé à \(\beta_0(\phi)\). |
| `correction_eta_phi` | Facteur correctif associé à \(\eta_{\mathrm{sweep}}(\phi)\). |

### Fichier sauvegardé

```matlab
pmerge_phi_th_results.mat
```

---

# Sous-dossier `Probabilité rupture`

## `Probabilité rupture/chi_bridge_th.m`

### Objectif

Calcule une approximation analytique de la fraction de liens localement critiques dans un Walker-Delta à uniformité orbitale.

Le script calcule aussi :

- le nombre moyen de liens par composante ;
- le nombre moyen de ponts par composante ;
- une approximation de \(\beta_0\) incluant les isolés, dimères et trimères.

### Fichiers d’entrée

```matlab
../Paramètres/betti_results.mat
../Paramètres/plink_results.mat
```

### Paramètres principaux

| Variable | Valeur | Description |
|---|---:|---|
| `R_earth` | `6371` km | Rayon terrestre. |
| `h` | `550` km | Altitude. |
| `N` | `204` | Nombre de satellites. |
| `i_deg` | `58` degrés | Inclinaison commune. |
| `dmax` | `1500` km | Portée de communication. |
| `epsilon_lat_deg` | `0.10` degré | Régularisation près des latitudes de retournement. |
| `n_phi` | `4000` | Résolution en latitude. |
| `n_r` | `1500` | Résolution sur la longueur des liens. |

### Sorties principales

Les résultats sont regroupés dans la structure `results`.

| Champ | Description |
|---|---|
| `phi` | Grille de latitude régularisée. |
| `lambda_orb` | Densité orbitale locale renormalisée. |
| `chi_bridge_local` | Fraction locale de liens sans voisin commun. |
| `chi_bridge` | Fraction orbitale moyenne. |
| `p_link_delta` | Probabilité globale de lien. |
| `E_edges` | Nombre moyen total de liens. |
| `E_N1`, `E_N2`, `E_N3` | Petites composantes théoriques utilisées. |
| `C_macro` | Nombre de composantes macroscopiques supposé. |
| `beta0_used` | Approximation de \(\beta_0\). |
| `mean_links_per_component` | Nombre moyen de liens par composante. |
| `mean_bridges_per_component` | Nombre moyen de ponts par composante. |

### Figures produites

1. densité orbitale locale régularisée ;
2. fraction locale de liens sans voisin commun.

### Fichier sauvegardé

```matlab
chi_bridge_results.mat
```

---

## `Probabilité rupture/pbreak_emp.m`

### Objectif

Calcule la probabilité empirique de rupture topologique à partir du barcode zigzag \(H_0\).

Pour la transition

\[
G_k
\longrightarrow
G_k\cup G_{k+1}
\longleftarrow
G_{k+1},
\]

une rupture correspond à une barre qui naît dans \(G_{k+1}\), à l’indice zigzag \(2k+1\).

Le code utilise la normalisation

\[
p_{\mathrm{break}}(k)
=
\frac{\text{nombre de ruptures à la transition }k}
{\beta_0(G_{k+1})}.
\]

### Type

Script principal sans fonction locale.

### Fichiers d’entrée

Le script cherche en priorité :

```matlab
../barcodes_results.mat
../analysis_temp_results.mat
```

Il prévoit aussi une recherche dans son dossier courant et une compatibilité avec l’ancien fichier :

```matlab
leo_H0_zigzag_barcodes_delta.mat
```

### Variables chargées

Depuis le barcode :

| Variable | Description |
|---|---|
| `birth_index` | Indices de naissance des barres. |
| `ZigzagTime` | Temps associés aux objets zigzag. |
| `h0_dims` | Dimensions successives de \(H_0\). |

Depuis l’analyse temporelle :

| Variable | Description |
|---|---|
| `time_values` | Instants des graphes réels. |
| `dt` | Pas temporel. |
| `inc_deg` | Inclinaison. |
| `N` | Nombre de satellites. |

### Sorties principales

| Variable | Description |
|---|---|
| `break_count` | Nombre de ruptures à chaque transition. |
| `break_count_from_beta0` | Comptage obtenu par différence de \(\beta_0\). |
| `beta0_union` | Nombre de composantes dans le graphe union. |
| `beta0_after` | Nombre de composantes dans \(G_{k+1}\). |
| `p_break` | Probabilité empirique instantanée. |
| `p_break_moving` | Moyenne glissante. |
| `break_count_moving` | Moyenne glissante du nombre de ruptures. |
| `p_break_mean` | Moyenne globale pondérée. |
| `p_break_time_mean` | Moyenne temporelle simple. |
| `max_break_error` | Erreur maximale entre les deux comptages. |
| `time_transition` | Temps au milieu des transitions. |

La moyenne globale pondérée est

\[
\overline p_{\mathrm{break}}
=
\frac{
\sum_k\text{break\_count}(k)
}{
\sum_k\beta_0(G_{k+1})
}.
\]

### Figures produites

1. probabilité empirique instantanée, moyenne glissante et moyenne globale ;
2. nombre de ruptures entre graphes consécutifs.

### Fichier sauvegardé

```matlab
pbreak_emp_results.mat
```

Le fichier est sauvegardé dans le même dossier que le script.

---

## `Probabilité rupture/pbreak_temp.m`

### Objectif

Compare temporellement la probabilité empirique de rupture à la valeur théorique constante calculée par `pbreak_th.m`.

Le script affiche :

- \(p_{\mathrm{break}}^{\mathrm{emp}}(t)\) ;
- sa moyenne glissante ;
- sa moyenne empirique globale ;
- \(p_{\mathrm{break}}^\Delta\).

### Type

Script principal sans fonction locale.

### Fichiers d’entrée

```matlab
pbreak_emp_results.mat
pbreak_th_results.mat
```

### Variables chargées

Depuis le fichier empirique :

| Variable | Description |
|---|---|
| `time_transition` | Temps des transitions. |
| `p_break` | Probabilité empirique instantanée. |
| `p_break_moving` | Moyenne glissante. |
| `p_break_mean` | Moyenne empirique globale pondérée. |
| `p_break_time_mean` | Moyenne temporelle simple. |
| `moving_window` | Taille de la fenêtre glissante. |

Depuis le fichier théorique :

| Variable | Description |
|---|---|
| `p_break_delta` | Probabilité théorique probabiliste. |
| `p_break_delta_linear` | Approximation théorique linéaire. |

### Sorties

| Variable | Description |
|---|---|
| `p_break_emp` | Courbe empirique instantanée. |
| `p_break_emp_moving` | Moyenne glissante empirique. |
| `p_break_emp_mean` | Moyenne empirique globale. |
| `p_break_emp_time_mean` | Moyenne temporelle simple. |
| `p_break_th` | Valeur théorique probabiliste. |
| `p_break_th_linear` | Approximation théorique linéaire. |

### Figure produite

Une figure superposant les valeurs empiriques et théoriques de \(p_{\mathrm{break}}\).

### Affichage console

Le script affiche :

- la valeur théorique probabiliste ;
- l’approximation linéaire ;
- la moyenne empirique globale ;
- la moyenne temporelle simple ;
- l’écart absolu ;
- le rapport théorie/empirique.

---

## `Probabilité rupture/pbreak_th.m`

### Objectif

Calcule la probabilité théorique de rupture d’une composante dans le modèle Walker-Delta.

Le calcul distingue :

1. la rupture d’un lien quelconque ;
2. la rupture conditionnelle d’un pont ;
3. le passage du niveau du pont au niveau de la composante.

### Type

Script principal sans fonction locale.

### Fichiers d’entrée

```matlab
../analysis_temp_results.mat
chi_bridge_results.mat
```

Le premier fournit notamment \(N\). Le second fournit les facteurs topologiques et les profils orbitaux.

### Paramètres principaux

| Variable | Valeur | Description |
|---|---:|---|
| `R_earth` | `6371` km | Rayon terrestre. |
| `h` | `550` km | Altitude. |
| `d_max` | `1500` km | Portée de communication. |
| `Delta_t` | `20` s | Pas temporel. |
| `v_rel` | \(v_{\mathrm{orb}}/\sqrt2\) | Vitesse relative moyenne. |
| `v_rad` | \(v_{\mathrm{rel}}/\pi\) | Projection radiale retenue. |

### Sorties principales

| Variable | Description |
|---|---|
| `v_orb` | Vitesse orbitale. |
| `v_rel` | Vitesse relative moyenne. |
| `v_rad` | Vitesse radiale retenue. |
| `chi_bridge_delta` | Fraction moyenne de liens critiques. |
| `p_bridge_at_boundary` | Probabilité qu’un lien à la frontière soit un pont. |
| `bridge_boundary_factor` | Facteur conditionnel frontière. |
| `mean_links_per_component` | Liens moyens par composante. |
| `mean_bridges_per_component` | Ponts moyens par composante. |
| `q_break_link` | Rupture d’un lien quelconque. |
| `q_break_bridge` | Rupture conditionnelle d’un pont. |
| `p_break_delta` | Probabilité théorique probabiliste. |
| `p_break_delta_linear` | Approximation linéaire. |
| `p_break_delta_linear_simplified` | Forme linéaire simplifiée. |
| `p_break_vs_bridge_factor` | Évolution en fonction du nombre de ponts. |
| `p_break_vs_dt` | Évolution en fonction du pas temporel. |

### Figures produites

1. \(p_{\mathrm{break}}^\Delta\) en fonction du nombre moyen de ponts par composante ;
2. \(p_{\mathrm{break}}^\Delta\) en fonction du pas temporel.

### Fichier sauvegardé

```matlab
pbreak_th_results.mat
```

## `Probabilité rupture/pbreak_phi_th.m`

### Objectif

Calcule la **probabilité théorique locale de rupture**

\[
p_{\mathrm{break}}^\Delta(\phi)
\]

dans le modèle Walker Delta à uniformité orbitale, en combinant :

1. la probabilité locale de rupture d’un lien ;
2. le nombre moyen de liens par composante non isolée ;
3. la probabilité qu’un lien situé au bord soit un pont ;
4. la fraction locale de composantes non isolées.

Le modèle local est

\[
p_{\mathrm{break}}^\Delta(\phi)
=
\frac{B_0(\phi)-N_1(\phi)}{B_0(\phi)}
\left[
1-
\exp\!\left(
-\mu_{\mathrm{break}}^{\mathrm{non-isolé}}(\phi)
\right)
\right],
\]

avec

\[
\mu_{\mathrm{break}}^{\mathrm{non-isolé}}(\phi)
=
\frac{E(\phi)}
{B_0(\phi)-N_1(\phi)}
\,p_{\mathrm{break}}^{\mathrm{lien}}(\phi)
\,p_{\mathrm{bridge,bord}}(\phi).
\]

### Type

Script principal avec plusieurs fonctions locales pour l’interpolation, la reconstruction empirique de \(\beta_0(\phi)\) et le calcul des versions corrigées.

### Fichiers d’entrée

```matlab
edges_phi_results.mat
N1_phi_results.mat
betti_phi_results.mat
densite_phi_results.mat
eta_sweep_phi_results.mat
analysis_temp_results.mat
```

### Rupture locale d’un lien

La vitesse relative locale est

\[
v_{\mathrm{rel}}^\Delta(\phi)
=
2v_{\mathrm{orb}}
\frac{\sqrt{\sin^2i-\sin^2\phi}}
{\cos\phi},
\]

et la probabilité de rupture d’un lien est modélisée par

\[
p_{\mathrm{break}}^{\mathrm{lien}}(\phi)
=
\frac{2}{\pi}
\frac{
v_{\mathrm{rel}}^\Delta(\phi)\,dt
}{
d_{\max}
}.
\]

### Probabilité qu’un lien de bord soit un pont

La version théorique utilise

\[
p_{\mathrm{bridge,bord}}^\Delta(\phi)
\simeq
\exp\!\left[
-\lambda(\phi)
A_{\mathrm{intersection}}(d_{\max})
\right],
\]

avec

\[
A_{\mathrm{intersection}}(d_{\max})
=
\left(
\frac{2\pi}{3}-\frac{\sqrt3}{2}
\right)d_{\max}^2.
\]

Deux corrections empiriques sont également comparées :

- \(\eta_{\mathrm{sweep}}^{\mathrm{emp}}(\phi)\) utilisé comme approximation de \(p_{\mathrm{bridge,bord}}(\phi)\) ;
- la vraie probabilité empirique \(P(\text{pont}\mid\text{rupture},\phi)\).

### Moyenne globale

La probabilité globale est obtenue par une moyenne pondérée par le nombre de composantes :

\[
p_{\mathrm{break}}^\Delta
=
\frac{
\sum_b
p_{\mathrm{break}}^\Delta(\phi_b)\,
B_{0,b}
}{
\sum_b B_{0,b}
}.
\]

### Versions comparées

Le script calcule :

- la forme probabiliste théorique ;
- son approximation linéaire ;
- une version corrigée avec \(\beta_0(\phi)\) empirique et \(\eta_{\mathrm{sweep}}^{\mathrm{emp}}(\phi)\) ;
- une version corrigée avec \(\beta_0(\phi)\) empirique et la vraie \(P(\text{pont}\mid\text{rupture},\phi)\).

### Sorties principales

| Variable | Description |
|---|---|
| `phi_vals` | Grille de latitude. |
| `edges_density_phi` | Densité locale théorique d’arêtes. |
| `N1_density_phi` | Densité locale de composantes isolées. |
| `betti0_density_phi` | Densité locale théorique de composantes. |
| `nonisolated_density_phi` | Densité locale de composantes non isolées. |
| `lambda_phi` | Densité surfacique locale. |
| `v_rel_phi` | Vitesse relative locale. |
| `p_break_link_phi` | Probabilité locale de rupture d’un lien. |
| `p_bridge_bord_phi` | Modèle théorique de probabilité de pont au bord. |
| `mean_links_per_nonisolated_component_phi` | Nombre moyen de liens par composante non isolée. |
| `mu_break_nonisolated_phi` | Nombre moyen d’événements de rupture critiques par composante non isolée. |
| `p_break_phi_th` | Probabilité locale théorique de rupture. |
| `p_break_phi_th_linear` | Approximation locale linéaire. |
| `p_break_th` | Probabilité globale théorique. |
| `p_break_phi_th_corrected_eta` | Correction locale utilisant \(\eta_{\mathrm{sweep}}^{\mathrm{emp}}\). |
| `p_break_phi_th_corrected_true` | Correction locale utilisant la vraie probabilité de pont. |
| `p_break_th_corrected_eta` | Probabilité globale corrigée avec \(\eta_{\mathrm{sweep}}^{\mathrm{emp}}\). |
| `p_break_th_corrected_true` | Probabilité globale corrigée avec \(P(\text{pont}\mid\text{rupture},\phi)\). |
| `correction_betti_phi` | Facteur correctif lié au modèle de \(\beta_0\). |
| `correction_p_bridge_eta_phi` | Facteur correctif lié à \(\eta_{\mathrm{sweep}}\). |
| `correction_p_bridge_true_phi` | Facteur correctif lié à la vraie probabilité de pont. |

### Fichier sauvegardé

```matlab
pbreak_phi_th_results.mat
```
