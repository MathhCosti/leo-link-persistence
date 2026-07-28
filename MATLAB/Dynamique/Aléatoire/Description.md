# Description des scripts MATLAB

Ce document décrit les fichiers MATLAB, ainsi que leurs fonctions locales.

---

---

## `analysis_temp.m`

### Objectif

Simule un réseau LEO dynamique avec des directions tangentielles aléatoires, construit le graphe des liens à chaque instant et calcule plusieurs métriques topologiques.

Le script :

- tire le nombre de satellites selon un processus de Poisson ;
- génère des positions initiales uniformes sur la sphère ;
- attribue à chaque satellite une direction tangentielle aléatoire ;
- fait évoluer les satellites sur des grands cercles ;
- construit les graphes temporels ;
- calcule \(\beta_0\), \(\beta_1\), la plus grande composante et le nombre de liens ;
- construit la suite zigzag par unions ;
- détecte les ponts exacts des graphes union ;
- estime une probabilité empirique de rupture des ponts ;
- sauvegarde les résultats pour `barcodes.m`.

La suite zigzag utilisée est

\[
G_1 \to G_1\cup G_2 \leftarrow G_2 \to G_2\cup G_3 \leftarrow G_3 \to \cdots
\]

### Type

Script principal avec une fonction locale :

- `find_bridges_tarjan`

### Entrées du script

Les paramètres sont définis directement dans le fichier.

| Variable | Valeur | Description |
|---|---:|---|
| `R_earth` | `6371` km | Rayon terrestre. |
| `h` | `550` km | Altitude orbitale. |
| `mu` | `398600` km³/s² | Paramètre gravitationnel terrestre. |
| `lambda` | `4e-7` sat/km² | Intensité du processus de Poisson. |
| `dmax` | `1500` km | Distance maximale d’un lien. |
| `dt` | `3` s | Pas temporel. |
| `Tmax` | `1500` s | Durée totale simulée. |

Le nombre de satellites suit

\[
N\sim\mathcal P\!\left(\lambda 4\pi R^2\right).
\]

### Modèle de mouvement

Les positions initiales sont uniformes sur la sphère. Un vecteur aléatoire est projeté sur le plan tangent :

\[
\mathbf v=\mathbf a-(\mathbf a\cdot\mathbf r_0)\mathbf r_0.
\]

Après normalisation, la position est

\[
\mathbf r(t)=R\left[\mathbf r_0\cos(\omega t)+\mathbf v\sin(\omega t)\right],
\qquad
\omega=\sqrt{\mu/R^3}.
\]

### Sorties principales

| Variable | Description |
|---|---|
| `Positions` | Positions des satellites à chaque instant. |
| `Adjacency` | Matrices d’adjacence temporelles. |
| `beta0` | Nombre de composantes connexes. |
| `beta1_graph` | Nombre cyclomatique \(E-N+\beta_0\). |
| `largest_component` | Taille de la plus grande composante. |
| `num_edges` | Nombre de liens. |
| `ZigzagAdjacency` | Matrices d’adjacence de la suite zigzag. |
| `ZigzagLabels` | Labels entiers ou demi-entiers des objets zigzag. |
| `beta0_zigzag` | Nombre de composantes dans la suite zigzag. |
| `beta1_zigzag_graph` | Nombre cyclomatique dans la suite zigzag. |
| `n_bridges_union` | Nombre de ponts dans chaque graphe union. |
| `n_removed_bridges_union` | Nombre de ponts retirés à l’étape suivante. |
| `mean_bridges_per_exposed_component` | Nombre moyen de ponts par composante exposée. |
| `q_break_given_bridge_global` | Probabilité empirique qu’un pont disparaisse. |
| `p_break_from_bridges` | Probabilité de rupture reconstruite. |
| `p_break_from_bridges_linear` | Approximation linéaire correspondante. |

### Ponts exposés

Le script calcule

\[
\overline B_{\mathrm{exposé}}
=
\frac{\sum_k B_k}{\sum_k\beta_0(U_k)}
\]

et

\[
q_{\mathrm{break}\mid\mathrm{bridge}}
=
\frac{\text{ponts retirés}}{\text{ponts exposés}}.
\]

La probabilité reconstruite est

\[
p_{\mathrm{break}}
=
1-\left(1-q_{\mathrm{break}\mid\mathrm{bridge}}\right)^{\overline B_{\mathrm{exposé}}}.
\]

### Figures produites

1. \(\beta_0(t)\) ;
2. \(\beta_1(t)\) du graphe ;
3. fraction de satellites dans la plus grande composante ;
4. nombre de liens ;
5. nombre moyen de ponts par composante exposée ;
6. probabilité conditionnelle de rupture d’un pont.

### Fichier sauvegardé

```matlab
analysis_temp_results.mat
```

---

## `anim_3D.m`

### Objectif

Anime en trois dimensions le même réseau LEO à directions tangentielles aléatoires.

Le script affiche en temps réel :

- les satellites ;
- les liens intersatellites ;
- l’instant courant ;
- le nombre de satellites ;
- le nombre de liens.

### Entrées du script

| Variable | Valeur | Description |
|---|---:|---|
| `R_earth` | `6371` km | Rayon terrestre. |
| `h` | `550` km | Altitude orbitale. |
| `lambda` | `4e-7` sat/km² | Intensité du processus de Poisson. |
| `dmax` | `1500` km | Portée des liens. |
| `dt` | `30` s | Pas temporel. |
| `Tmax` | `6000` s | Durée simulée. |

### Sorties

Le script ne retourne pas de valeur et ne sauvegarde aucun fichier.

Les principales variables du workspace sont :

| Variable | Description |
|---|---|
| `positions0` | Positions initiales. |
| `r0` | Vecteurs radiaux unitaires. |
| `v` | Directions tangentielles. |
| `positions_t` | Positions au dernier instant. |
| `A` | Dernière matrice d’adjacence. |
| `E` | Nombre de liens au dernier instant. |
| `sat_handle` | Objet graphique des satellites. |
| `link_handle` | Objet graphique des liens. |

---

## `barcodes.m`

### Objectif

Calcule le barcode zigzag en homologie \(H_0\) à partir des données produites par `analysis_temp.m`.

Le script :

- charge la suite zigzag ;
- calcule les composantes connexes ;
- construit les applications induites en \(H_0\) ;
- décompose le module zigzag en intervalles ;
- calcule les temps de naissance, de mort et les durées ;
- supprime la barre globale persistante ;
- compare la survie empirique à un modèle exponentiel ;
- calcule \(p_{\mathrm{merge}}\), \(p_{\mathrm{break}}\) et \(p_{\mathrm{disp}}\) ;
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

### Fichier d’entrée

```matlab
analysis_temp_results.mat
```

### Sorties principales

| Variable | Description |
|---|---|
| `intervals` | Indices de naissance et de mort. |
| `birth_index` | Indices de naissance. |
| `death_index` | Indices de mort. |
| `birth_time` | Temps de naissance. |
| `death_time` | Temps de mort. |
| `lifetimes` | Durées de vie. |
| `ZigzagTime` | Temps associés aux objets zigzag. |
| `h0_dims` | Dimensions successives de \(H_0\). |
| `p_merge` | Probabilité corrigée de fusion. |
| `q_break` | Probabilité corrigée de rupture. |
| `p_death` | Probabilité totale de disparition. |
| `tau_th` | Temps caractéristique théorique. |
| `p_disp_emp_t` | Probabilité empirique temporelle de disparition. |
| `p_disp_emp_global` | Estimateur global de \(p_{\mathrm{disp}}\). |

### Extraction empirique de \(p_{\mathrm{disp}}\)

Pour chaque intervalle temporel, le script calcule

\[
p_{\mathrm{disp}}^{\mathrm{emp}}(t_k)
=
\frac{\text{barres mourant avant }t_{k+1}}
{\text{barres vivantes à }t_k}.
\]

Il produit :

- `p_disp_emp_t` ;
- `p_disp_emp_mean` ;
- `p_disp_emp_global` ;
- `p_disp_emp_from_mean_lifetime`.

### Figures produites

1. survie empirique et modèle exponentiel ;
2. évolution temporelle de \(p_{\mathrm{disp}}\) ;
3. barcode zigzag \(H_0\) ;
4. histogramme des durées de vie.

### Fichier sauvegardé

```matlab
barcodes_results.mat
```

---

---

# Sous-dossier `Probabilité fusion`

## `Probabilité fusion/chi_temp_lambda.m`

### Objectif

Étudie le facteur correctif topologique utilisé dans le modèle de fusion en fonction de la densité satellitaire \(\lambda\).

Le facteur étudié est

\[
\chi_{\mathrm{merge}}
=
\frac{N-\beta_0}{|E|}.
\]

Il mesure la fraction des arêtes qui participent effectivement à une réduction du nombre de composantes connexes. Dans une forêt, chaque arête fusionne deux composantes et le facteur vaut 1. Lorsque des cycles apparaissent, une partie des arêtes devient topologiquement redondante et le facteur diminue.

Le script compare :

- la moyenne empirique du ratio calculé réalisation par réalisation ;
- le ratio construit à partir des moyennes ;
- une approximation théorique sparse ;
- une approximation théorique fondée sur les satellites isolés.

### Type

Script principal avec trois fonctions locales :

- `compute_betti_0_1`
- `rank_mod2`
- `edge_key`

### Entrées du script

| Variable | Type | Description |
|---|---:|---|
| `Re` | réel | Rayon terrestre en kilomètres. |
| `h` | réel | Altitude orbitale. |
| `R` | réel | Rayon orbital. |
| `d_max` | réel | Distance maximale de connexion. |
| `alpha_max` | réel | Angle maximal de connexion. |
| `lambda_scaled_values` | vecteur réel | Densités testées en satellites par \(10^6\ \mathrm{km}^2\). |
| `n_iter` | entier | Nombre de réalisations Monte-Carlo. |

### Sorties du script

| Variable | Description |
|---|---|
| `N_values` | Nombre de satellites associé à chaque densité. |
| `Betti0_all` | Valeurs simulées de \(\beta_0\). |
| `Betti1_graph_all` | Valeurs simulées de \(\beta_1\) du graphe. |
| `Betti1_complex_all` | Valeurs simulées de \(\beta_1\) du complexe de clique. |
| `E_all` | Nombre d’arêtes pour chaque réalisation. |
| `Frac_bridge_all` | Valeurs de \((N-\beta_0)/|E|\) pour chaque réalisation. |
| `Betti0` | Moyenne empirique de \(\beta_0\). |
| `E_mean` | Nombre moyen empirique d’arêtes. |
| `Frac_bridge` | Moyenne des ratios réalisation par réalisation. |
| `Frac_bridge_ratio_means` | Ratio \((N-\mathbb E[\beta_0])/\mathbb E[|E|]\). |
| `Frac_bridge_theory_sparse` | Facteur théorique sparse. |
| `Frac_bridge_theory_isolated` | Facteur théorique fondé sur les isolés. |

### Figures produites

1. \(\beta_0\) moyen en fonction de \(\lambda\) ;
2. \(\beta_1^{\mathrm{graphe}}\) moyen ;
3. nombre moyen d’arêtes ;
4. facteur correctif \((N-\beta_0)/|E|\).

Le script affiche également un tableau récapitulatif dans la console.

---

## `Probabilité fusion/pmerge_emp.m`

### Objectif

Calcule la probabilité empirique de fusion topologique \(p_{\mathrm{merge}}\) à partir du barcode zigzag \(H_0\).

Pour une transition

\[
G_k\longrightarrow U_k\longleftarrow G_{k+1},
\qquad
U_k=G_k\cup G_{k+1},
\]

une fusion correspond à une barre qui meurt à l’indice zigzag associé à \(G_k\), soit \(2k-1\).

La probabilité empirique est

\[
p_{\mathrm{merge}}(k)
=
\frac{\text{nombre de fusions à la transition }k}
{\beta_0(G_k)}.
\]

### Fichiers d’entrée

```matlab
barcodes_results.mat
analysis_temp_results.mat
```

Le premier fichier est obligatoire. Le second permet de récupérer les instants réels et le pas temporel.

### Variables chargées

Depuis `barcodes_results.mat` :

| Variable | Description |
|---|---|
| `death_index` | Indices de mort des barres. |
| `ZigzagTime` | Temps associés aux objets zigzag. |
| `h0_dims` | Dimensions de \(H_0\), donc nombres de composantes. |

Depuis `analysis_temp_results.mat` :

| Variable | Description |
|---|---|
| `time_values` | Instants des graphes réels. |
| `dt` | Pas temporel. |

### Sorties du script

| Variable | Description |
|---|---|
| `merge_count` | Nombre de fusions à chaque transition. |
| `merge_count_from_beta0` | Nombre de fusions obtenu par différence de \(\beta_0\). |
| `beta0_before` | Nombre de composantes avant la fusion. |
| `beta0_union` | Nombre de composantes dans le graphe union. |
| `p_merge` | Probabilité empirique instantanée. |
| `p_merge_moving` | Moyenne glissante de \(p_{\mathrm{merge}}\). |
| `merge_count_moving` | Moyenne glissante du nombre de fusions. |
| `p_merge_mean` | Moyenne globale pondérée. |
| `p_merge_time_mean` | Moyenne temporelle simple. |
| `time_transition` | Temps au milieu des transitions. |
| `moving_window` | Taille de la fenêtre glissante. |

La moyenne globale pondérée est

\[
\overline p_{\mathrm{merge}}
=
\frac{\sum_k\text{merge\_count}(k)}
{\sum_k\beta_0(G_k)}.
\]

### Figures produites

1. probabilité empirique instantanée, moyenne glissante et moyenne globale ;
2. nombre de fusions par transition.

### Fichier sauvegardé

```matlab
pmerge_emp_results.mat
```

---

## `Probabilité fusion/pmerge_temp.m`

### Objectif

Compare temporellement la probabilité empirique de fusion à la valeur théorique constante calculée par `pmerge_th.m`.

Le script affiche :

- \(p_{\mathrm{merge}}^{\mathrm{emp}}(t)\) ;
- sa moyenne glissante ;
- la moyenne empirique globale ;
- \(p_{\mathrm{merge}}^{\mathrm{th}}\).

### Type

Script principal sans fonction locale.

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
| `p_merge_moving` | Moyenne glissante. |
| `p_merge_mean` | Moyenne empirique globale pondérée. |
| `p_merge_time_mean` | Moyenne temporelle simple. |
| `moving_window` | Taille de la fenêtre glissante. |

Depuis le fichier théorique :

| Variable | Description |
|---|---|
| `p_merge_th` | Probabilité théorique probabiliste. |
| `p_merge_th_linear` | Approximation théorique linéaire. |

### Sorties

| Variable | Description |
|---|---|
| `p_merge_emp` | Courbe empirique instantanée. |
| `p_merge_emp_moving` | Moyenne glissante empirique. |
| `p_merge_emp_mean` | Moyenne empirique globale. |
| `p_merge_emp_time_mean` | Moyenne temporelle simple. |
| `p_merge_th` | Valeur théorique probabiliste. |
| `p_merge_th_linear` | Approximation linéaire. |

### Figure produite

Une figure comparant les probabilités empirique et théorique de fusion.

### Affichage console

Le script affiche :

- la valeur théorique probabiliste ;
- la valeur théorique linéaire ;
- la moyenne empirique globale ;
- la moyenne temporelle simple ;
- l’écart absolu théorie/empirique ;
- le rapport théorie/empirique.

---

## `Probabilité fusion/pmerge_th.m`

### Objectif

Calcule théoriquement la probabilité de fusion \(p_{\mathrm{merge}}\) dans le modèle aléatoire à directions tangentielles.

Le calcul reprend la partie théorique de `barcodes.m` et combine :

1. une aire balayée pendant un pas temporel ;
2. une loi de Poisson spatiale ;
3. un facteur topologique \(\chi_{\mathrm{merge}}\).

### Type

Script principal sans fonction locale.

### Fichier d’entrée

```matlab
analysis_temp_results.mat
```

Variables utilisées :

| Variable | Description |
|---|---|
| `N` | Nombre de satellites. |
| `R` | Rayon orbital. |
| `lambda` | Densité satellitaire. |
| `dmax` | Distance maximale de connexion. |
| `dt` | Pas temporel. |

### Sorties du script

| Variable | Description |
|---|---|
| `v_orb` | Vitesse orbitale. |
| `v_rel` | Vitesse relative moyenne. |
| `A_sweep` | Aire balayée brute. |
| `p_merge_raw` | Probabilité brute de fusion. |
| `alpha_max` | Angle maximal de connexion. |
| `p_link` | Probabilité de lien. |
| `E_theory` | Nombre moyen d’arêtes. |
| `N1_theory` | Nombre moyen d’isolés. |
| `N2_theory` | Nombre moyen de dimères. |
| `N3_theory` | Nombre moyen de trimères. |
| `beta0_theory` | Approximation du nombre de composantes. |
| `chi_merge` | Facteur topologique de fusion. |
| `A_sweep_corrected` | Aire balayée corrigée. |
| `p_merge_th` | Probabilité théorique probabiliste. |
| `p_merge_th_linear` | Approximation linéaire. |

### Fichier sauvegardé

```matlab
pmerge_th_results.mat
```

---

---

# Sous-dossier `Probabilité rupture`

## `Probabilité rupture/chi_bridge_emp.m`

### Objectif

Calcule empiriquement les grandeurs topologiques utilisées pour corriger la probabilité de rupture \(p_{\mathrm{break}}\).

Le script exécute d’abord `analysis_temp.m`, puis analyse tous les graphes temporels produits afin d’estimer :

- le nombre moyen de liens par composante ;
- la proportion de liens qui sont des ponts ;
- le nombre moyen de ponts, ou liens délétères, par composante ;
- une estimation alternative fondée sur les moments de la taille d’une composante.

Le facteur recommandé est

\[
F_{\mathrm{liens}}
=
\frac{\text{nombre total de liens}}
{\text{nombre total de composantes}}.
\]

La proportion empirique de liens critiques est

\[
p_{\mathrm{link},A}
=
\frac{\text{nombre total de ponts}}
{\text{nombre total de liens}},
\]

et le nombre moyen de liens délétères par composante vaut

\[
F_{\mathrm{délétère}}
=
F_{\mathrm{liens}}\,p_{\mathrm{link},A}
=
\frac{\text{nombre total de ponts}}
{\text{nombre total de composantes}}.
\]

### Type

Script principal avec une fonction locale :

- `count_bridges_tarjan`

### Dépendance principale

Le script cherche et exécute :

```matlab
../analysis_temp.m
```

### Entrées du script

Le script ne possède pas d’argument d’entrée. Les données sont produites par `analysis_temp.m`.

### Sorties du script

| Variable | Description |
|---|---|
| `E_NA` | Taille moyenne d’une composante observée. |
| `E_NA2` | Moment d’ordre deux de la taille d’une composante. |
| `p_link` | Probabilité moyenne empirique qu’une paire soit liée. |
| `p_link_A` | Proportion de liens qui sont des ponts. |
| `N_liens_moyen_emp` | Nombre moyen réel de liens par composante. |
| `N_liens_deleteres` | Nombre moyen de ponts par composante. |
| `N_liens_moyen_formule` | Estimation du nombre moyen de liens obtenue par la formule utilisant les moments. |
| `N_deleteres_formule` | Estimation correspondante du nombre de liens délétères. |
| `n_components_t` | Nombre de composantes à chaque instant. |
| `n_edges_t` | Nombre de liens à chaque instant. |
| `n_bridges_t` | Nombre de ponts à chaque instant. |
| `p_link_t` | Probabilité empirique de lien à chaque instant. |

### Figures produites

1. nombre moyen de liens par composante en fonction du temps ;
2. nombre moyen de ponts par composante en fonction du temps.

### Fichier sauvegardé

```matlab
chi_bridge_emp_results.mat
```

---

## `Probabilité rupture/pbreak_emp.m`

### Objectif

Calcule la probabilité empirique de rupture topologique \(p_{\mathrm{break}}\) à partir du barcode zigzag en homologie \(H_0\).

Pour une transition

\[
G_k\longrightarrow U_k\longleftarrow G_{k+1},
\qquad
U_k=G_k\cup G_{k+1},
\]

une rupture correspond à une barre qui naît à l’indice associé à \(G_{k+1}\), soit l’indice zigzag \(2k+1\).

La probabilité empirique est définie par

\[
p_{\mathrm{break}}(k)
=
\frac{\text{nombre de ruptures à la transition }k}
{\beta_0(U_k)}.
\]

### Fichiers d’entrée

```matlab
barcodes_results.mat
analysis_temp_results.mat
```

Le premier fichier est obligatoire. Le second est utilisé pour récupérer `time_values` et `dt` lorsqu’il est disponible.

### Variables chargées

Depuis `barcodes_results.mat` :

| Variable | Description |
|---|---|
| `birth_index` | Indices de naissance des barres. |
| `ZigzagTime` | Temps associés aux objets zigzag. |
| `h0_dims` | Dimensions de \(H_0\), donc nombres de composantes. |

Depuis `analysis_temp_results.mat` :

| Variable | Description |
|---|---|
| `time_values` | Instants des graphes réels. |
| `dt` | Pas temporel. |

### Sorties du script

| Variable | Description |
|---|---|
| `break_count` | Nombre de ruptures à chaque transition. |
| `break_count_from_beta0` | Nombre de ruptures obtenu par différence de \(\beta_0\). |
| `beta0_union` | Nombre de composantes dans les graphes union. |
| `beta0_after` | Nombre de composantes dans les graphes suivants. |
| `p_break` | Probabilité empirique instantanée. |
| `p_break_moving` | Moyenne glissante de \(p_{\mathrm{break}}\). |
| `break_count_moving` | Moyenne glissante du nombre de ruptures. |
| `p_break_mean` | Moyenne globale pondérée. |
| `p_break_time_mean` | Moyenne temporelle simple. |
| `time_transition` | Temps au milieu de chaque transition. |
| `moving_window` | Taille de la fenêtre glissante. |

La moyenne globale est

\[
\overline p_{\mathrm{break}}
=
\frac{\sum_k \text{break\_count}(k)}
{\sum_k\beta_0(U_k)}.
\]

### Figures produites

1. probabilité empirique instantanée et moyenne glissante ;
2. nombre de ruptures par transition.

### Fichier sauvegardé

```matlab
pbreak_emp_results.mat
```

---

## `Probabilité rupture/pbreak_temp.m`

### Objectif

Compare l’évolution temporelle empirique de \(p_{\mathrm{break}}\) à la valeur théorique constante calculée par `pbreak_th.m`.

Le script superpose :

- la courbe empirique instantanée ;
- sa moyenne glissante ;
- la moyenne empirique globale ;
- la valeur théorique probabiliste.

### Fichiers d’entrée

```matlab
pbreak_emp_results.mat
pbreak_th_results.mat
```

### Variables chargées

Depuis le fichier empirique :

| Variable | Description |
|---|---|
| `time_transition` | Temps associés aux transitions. |
| `p_break` | Probabilité empirique instantanée. |
| `p_break_moving` | Moyenne glissante empirique. |
| `p_break_mean` | Moyenne empirique globale pondérée. |
| `p_break_time_mean` | Moyenne temporelle simple. |
| `moving_window` | Taille de la fenêtre glissante. |

Depuis le fichier théorique :

| Variable | Description |
|---|---|
| `p_break_th` | Probabilité théorique sous forme probabiliste. |
| `p_break_th_linear` | Approximation linéaire. |

### Sorties

Ce script ne crée pas de nouvelle variable théorique. Il laisse dans le workspace les données chargées et les variables renommées :

| Variable | Description |
|---|---|
| `p_break_emp` | Courbe empirique instantanée. |
| `p_break_emp_moving` | Moyenne glissante empirique. |
| `p_break_emp_mean` | Moyenne globale empirique. |
| `p_break_emp_time_mean` | Moyenne temporelle simple. |
| `p_break_th` | Valeur théorique probabiliste. |
| `p_break_th_linear` | Valeur théorique linéaire. |

### Figure produite

Une figure comparant \(p_{\mathrm{break}}^{\mathrm{emp}}(t)\) et \(p_{\mathrm{break}}^{\mathrm{th}}\).

### Affichage console

Le script affiche :

- la valeur théorique probabiliste ;
- l’approximation linéaire ;
- la moyenne empirique globale ;
- la moyenne temporelle simple ;
- l’écart absolu ;
- le rapport théorie/empirique lorsque la moyenne empirique est non nulle.

---

## `Probabilité rupture/pbreak_th.m`

### Objectif

Calcule théoriquement la probabilité de rupture \(p_{\mathrm{break}}\) dans le modèle aléatoire à vecteurs tangentiels.

Le modèle combine :

1. une probabilité qu’un lien individuel sorte de la zone de connexion ;
2. une approximation de la fraction de liens critiques ;
3. un nombre moyen de ponts par composante ;
4. une composition probabiliste des ruptures possibles.

### Fichier d’entrée

Le script charge le fichier situé un niveau au-dessus :

```matlab
../analysis_temp_results.mat
```

Variables utilisées :

| Variable | Description |
|---|---|
| `N` | Nombre de satellites. |
| `R` | Rayon orbital. |
| `lambda` | Densité satellitaire. |
| `dmax` | Distance maximale de connexion. |
| `dt` | Pas temporel. |

### Sorties du script

| Variable | Description |
|---|---|
| `v_orb` | Vitesse orbitale. |
| `v_rel` | Vitesse relative moyenne utilisée. |
| `alpha_max` | Angle maximal de connexion. |
| `p_link` | Probabilité de lien. |
| `E_theory` | Nombre moyen d’arêtes. |
| `N1_theory` | Nombre moyen d’isolés. |
| `N2_theory` | Nombre moyen de dimères. |
| `N3_theory` | Nombre moyen de trimères. |
| `beta0_theory` | Approximation du nombre moyen de composantes. |
| `chi_bridge` | Fraction analytique de liens critiques. |
| `mean_links_per_component` | Nombre moyen de liens par composante. |
| `mean_bridges_per_component` | Nombre moyen de ponts par composante. |
| `p_break_link` | Probabilité de rupture d’un lien individuel. |
| `p_break_th` | Probabilité théorique probabiliste. |
| `p_break_th_linear` | Approximation linéaire. |

### Fichier sauvegardé

```matlab
pbreak_th_results.mat
```