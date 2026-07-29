# Description des scripts MATLAB

Ce document décrit les fichiers MATLAB, ainsi que leurs fonctions locales.

---

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

---

# Sous-dossier `Paramètres`

## `Paramètres/vitesse_rad_emp.m`

### Objectif

Calcule empiriquement la vitesse radiale relative des liens intersatellites à partir des positions temporelles sauvegardées dans `analysis_temp_results.mat`.

Pour un lien \((i,j)\), la vitesse radiale exacte est

\[
v_{\mathrm{rad},ij}
=
\frac{
(\mathbf r_j-\mathbf r_i)\cdot(\mathbf v_j-\mathbf v_i)
}{
\|\mathbf r_j-\mathbf r_i\|
}.
\]

Le signe indique :

- \(v_{\mathrm{rad}}>0\) : les satellites s’éloignent ;
- \(v_{\mathrm{rad}}<0\) : les satellites se rapprochent.

Le script calcule aussi la partie sortante

\[
(v_{\mathrm{rad}})_+
=
\max(v_{\mathrm{rad}},0),
\]

qui est la quantité pertinente pour les ruptures de liens.

### Fichier d’entrée

```matlab
analysis_temp_results.mat
```

Variables chargées :

| Variable | Description |
|---|---|
| `Positions` | Positions des satellites à chaque instant. |
| `Adjacency` | Matrices d’adjacence temporelles. |
| `time_values` | Instants de simulation. |
| `N` | Nombre de satellites. |
| `R` | Rayon orbital. |
| `dmax` | Distance maximale de connexion. |
| `dt` | Pas temporel. |
| `inc_deg` | Inclinaison, lorsqu’elle existe. |

### Reconstruction des vitesses

Les vitesses individuelles sont reconstruites par différences finies :

- différence avant au premier instant ;
- différence arrière au dernier instant ;
- différence centrée aux instants intermédiaires.

### Sorties principales

| Variable | Description |
|---|---|
| `mean_vrad_signed_t` | Vitesse radiale signée moyenne à chaque instant. |
| `mean_vrad_abs_t` | Valeur absolue moyenne. |
| `mean_vrad_out_t` | Partie positive moyenne sur tous les liens. |
| `mean_vrad_in_t` | Partie entrante moyenne. |
| `mean_vrad_out_boundary_t` | Vitesse sortante moyenne près de `dmax`. |
| `all_vrad_signed` | Toutes les vitesses radiales signées observées. |
| `all_distances` | Distances associées aux liens. |
| `fraction_outgoing` | Fraction de liens avec \(v_{\mathrm{rad}}>0\). |
| `fraction_boundary` | Fraction de liens proches de `dmax`. |
| `q_break_flux_all_links` | Estimation par flux de frontière. |
| `q_break_crossing_estimate` | Estimation directe des franchissements de `dmax`. |

L’approximation par flux est

\[
q_{\mathrm{break}}
\approx
\frac{2\,dt}{d_{\max}}
\mathbb E[(\dot D)_+\mid \text{lien}].
\]

L’estimation directe compte les liens pour lesquels

\[
D+\dot D\,dt>d_{\max}.
\]

### Figures produites

1. vitesse radiale sortante moyenne au cours du temps ;
2. vitesse sortante des liens proches de `dmax` ;
3. histogramme des vitesses radiales signées ;
4. vitesse radiale en fonction de la longueur du lien.

### Fichier sauvegardé

```matlab
vitesse_rad_emp_results.mat
```

---

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

## `Probabilité fusion/pmerge_th.m`

### Objectif

Calcule la probabilité théorique de fusion dans un modèle Walker-Delta à uniformité orbitale.

### Type

Script principal avec une fonction locale :

- `gauss_legendre`

### Fichiers d’entrée

Le script charge :

```matlab
../analysis_temp_results.mat
../Paramètres/plink_results.mat
../Paramètres/betti_results.mat
```

### Paramètres principaux

| Variable | Valeur | Description |
|---|---:|---|
| `R_earth` | `6371` km | Rayon terrestre. |
| `h` | `550` km | Altitude. |
| `inc_deg` | `58` degrés | Inclinaison utilisée. |
| `d_max` | `1500` km | Portée de connexion. |
| `Delta_t` | `20` s | Pas temporel du modèle. |
| `v_rel` | \(v_{\mathrm{orb}}/\sqrt2\) | Vitesse relative moyenne. |
| `n_quad` | `500` | Ordre de quadrature en phase. |

Le nombre \(N\) est chargé depuis `analysis_temp_results.mat`.

### Sorties principales

| Variable | Description |
|---|---|
| `v_orb` | Vitesse orbitale. |
| `v_rel` | Vitesse relative retenue. |
| `p_link_delta` | Probabilité moyenne de lien. |
| `E_edges_delta` | Nombre moyen d’arêtes. |
| `beta0_delta` | Approximation de \(\beta_0\). |
| `chi_delta_raw` | Facteur correctif avant bornage. |
| `chi_delta` | Facteur correctif utilisé. |
| `phi` | Grille de latitude. |
| `lambda_delta_phi` | Densité locale. |
| `p_merge_phi` | Probabilité locale de fusion. |
| `lambda_band_mean` | Densité moyenne sur la bande. |
| `p_merge_mean` | Moyenne orbitale théorique. |
| `p_merge_mean_density` | Approximation utilisant la densité moyenne. |

### Figures produites

1. probabilité locale de fusion en fonction de la latitude ;
2. densité locale utilisée dans le calcul.

### Fichier sauvegardé

```matlab
pmerge_th_results.mat
```

---

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