# Description des scripts et fonctions MATLAB

Ce document décrit les fichiers MATLAB, ainsi que leurs fonctions locales.

---

## `constellation.m`

### Objectif

Génère une réalisation d’un **processus ponctuel de Poisson homogène sur une sphère orbitale LEO**, puis affiche les positions des satellites en trois dimensions.

Le nombre de satellites est tiré selon

\[
N\sim\mathcal P\!\left(\lambda\,4\pi R^2\right),
\]

où

\[
R=R_{\mathrm{Earth}}+h.
\]

Les positions sont ensuite tirées uniformément sur la sphère de rayon \(R\).

### Type

Script principal sans fonction locale.

### Entrées du script

| Variable | Type | Description |
|---|---:|---|
| `R_earth` | réel | Rayon terrestre en kilomètres. |
| `h` | réel | Altitude orbitale en kilomètres. |
| `R` | réel | Rayon orbital \(R=R_{\mathrm{Earth}}+h\). |
| `lambda` | réel positif | Intensité du processus de Poisson en satellites/km². |
| `surface_sphere` | réel | Surface de la sphère orbitale. |

### Sorties du script

| Variable | Description |
|---|---|
| `N` | Nombre de satellites tiré. |
| `u` | Variables uniformes utilisées pour générer la colatitude. |
| `phi` | Longitudes des satellites. |
| `theta` | Colatitudes des satellites. |
| `x`, `y`, `z` | Coordonnées cartésiennes des satellites en kilomètres. |

### Figure produite

Le script produit une figure 3D contenant les satellites sous forme de points.

Le titre indique le nombre de satellites effectivement généré.

---

## `graphe_3D.m`

### Objectif

Construit le graphe géométrique des liens intersatellites à partir de positions déjà présentes dans le workspace, puis affiche les satellites et les liens en trois dimensions.

Deux satellites sont reliés lorsque leur distance euclidienne vérifie

\[
D_{ij}\leq d_{\max}.
\]

### Type

Script principal sans fonction locale.

### Dépendances

Le script suppose que les variables suivantes existent déjà dans le workspace :

| Variable | Description |
|---|---|
| `x` | Coordonnées cartésiennes selon l’axe \(x\). |
| `y` | Coordonnées cartésiennes selon l’axe \(y\). |
| `z` | Coordonnées cartésiennes selon l’axe \(z\). |

Il est donc destiné à être exécuté après `constellation.m` ou après tout autre script générant ces coordonnées.

### Entrées du script

| Variable | Type | Description |
|---|---:|---|
| `dmax` | réel positif | Distance maximale de connexion en kilomètres. |
| `x`, `y`, `z` | vecteurs réels | Coordonnées cartésiennes des satellites. |

### Sorties du script

| Variable | Description |
|---|---|
| `positions` | Matrice `N x 3` des positions. |
| `D` | Matrice complète des distances entre satellites. |
| `A` | Matrice d’adjacence logique. |
| `G` | Objet MATLAB `graph`. |
| `row`, `col` | Indices des extrémités des arêtes. |
| `E` | Nombre total d’arêtes. |
| `Xlinks`, `Ylinks`, `Zlinks` | Coordonnées utilisées pour tracer les segments des liens. |

### Figure produite

Le script produit une figure 3D contenant :

- les liens intersatellites sous forme de segments ;
- les satellites sous forme de points ;
- les axes cartésiens en kilomètres.

---

---

# Sous-dossier `Betti`

## `Betti/betti_alpha.m`

### Objectif

Étudie les nombres de Betti \(\beta_0\) et \(\beta_1\) en fonction de l’angle maximal de connexion \(\alpha_{\max}\), pour un nombre de satellites \(N\), une altitude \(h\) et une distance maximale \(d_{\max}\) fixés.

Le graphe de liens est construit avec la double condition

\[
D_{ij}\leq d_{\max}
\]

et

\[
\alpha_{ij}\leq \alpha_{\max}.
\]

L’angle effectivement utilisé est donc

\[
\alpha_{\mathrm{eff}}
=
\min(\alpha_{\max},\alpha_{d_{\max}}),
\]

où

\[
\alpha_{d_{\max}}
=
2\arcsin\left(\frac{d_{\max}}{2R}\right).
\]

Le script compare les valeurs simulées à plusieurs approximations théoriques de \(\beta_0\) et du premier nombre de Betti du graphe.

### Type

Script principal avec trois fonctions locales :

- `compute_betti_0_1`
- `rank_mod2`
- `edge_key`

### Entrées du script

| Variable | Type | Description |
|---|---:|---|
| `N` | entier | Nombre de satellites. |
| `Re` | réel | Rayon terrestre en kilomètres. |
| `h` | réel | Altitude orbitale en kilomètres. |
| `R` | réel | Rayon orbital \(R=R_e+h\). |
| `d_max` | réel | Distance maximale de connexion en kilomètres. |
| `alpha_values` | vecteur réel | Valeurs de \(\alpha_{\max}\) testées, en radians. |
| `n_iter` | entier | Nombre de réalisations Monte-Carlo. |

### Sorties du script

#### Résultats simulés

| Variable | Description |
|---|---|
| `Betti0_all` | Valeurs de \(\beta_0\) pour toutes les réalisations et tous les angles. |
| `Betti1_graph_all` | Valeurs de \(\beta_1\) du graphe. |
| `Betti1_complex_all` | Valeurs de \(\beta_1\) du complexe de clique. |
| `Betti0` | Moyenne Monte-Carlo de \(\beta_0\). |
| `Betti1_graph` | Moyenne Monte-Carlo de \(\beta_1^{\mathrm{graphe}}\). |
| `Betti1_complex` | Moyenne Monte-Carlo de \(\beta_1^{\mathrm{complexe}}\). |
| `Betti0_std` | Écart-type de \(\beta_0\). |
| `Betti1_graph_std` | Écart-type de \(\beta_1^{\mathrm{graphe}}\). |
| `Betti1_complex_std` | Écart-type de \(\beta_1^{\mathrm{complexe}}\). |

#### Résultats théoriques

| Variable | Description |
|---|---|
| `E_theory` | Nombre moyen théorique d’arêtes. |
| `I_theory` | Nombre moyen théorique de satellites isolés. |
| `Beta0_theory_sparse` | Approximation sparse de \(\beta_0\). |
| `Beta0_theory_isolated` | Approximation par satellites isolés. |
| `Beta0_theory_dimers_geom` | Approximation incluant les dimères. |
| `Beta0_theory_trimers_geom` | Approximation incluant les trimères. |
| `Beta1_graph_theory_sparse` | Approximation de \(\beta_1^{\mathrm{graphe}}\) issue du modèle sparse. |
| `Beta1_graph_theory_isolated` | Approximation issue du modèle par isolés. |
| `Beta1_graph_theory_dimers_geom` | Approximation incluant les dimères. |
| `Beta1_graph_theory_trimers_geom` | Approximation incluant les trimères. |
| `Beta1_complex_theory_ER` | Approximation Erdős-Rényi de \(\beta_1\) du complexe de clique. |

### Figures produites

1. \(\beta_0\) moyen en fonction de \(\alpha_{\max}\) ;
2. \(\beta_1^{\mathrm{graphe}}\) moyen en fonction de \(\alpha_{\max}\).

---

## `Betti/betti_dmax.m`

### Objectif

Étudie \(\beta_0\), \(\beta_1^{\mathrm{graphe}}\) et \(\beta_1^{\mathrm{complexe}}\) en fonction de la distance maximale de connexion \(d_{\max}\), pour \(N\), \(h\) et \(\alpha_{\max}\) fixés.

La contrainte angulaire peut saturer le balayage en distance. Le seuil de distance équivalent est

\[
d_{\alpha_{\max}}
=
2R\sin\left(\frac{\alpha_{\max}}{2}\right).
\]

Au-delà de cette distance, augmenter \(d_{\max}\) ne change plus le graphe si \(\alpha_{\max}\) reste la contrainte la plus restrictive.

### Type

Script principal avec trois fonctions locales :

- `compute_betti_0_1`
- `rank_mod2`
- `edge_key`

### Entrées du script

| Variable | Type | Description |
|---|---:|---|
| `N` | entier | Nombre de satellites. |
| `Re` | réel | Rayon terrestre en kilomètres. |
| `h` | réel | Altitude orbitale. |
| `R` | réel | Rayon orbital. |
| `alpha_max` | réel | Angle maximal fixé, en radians. |
| `dmax_values` | vecteur réel | Distances maximales testées, en kilomètres. |
| `n_iter` | entier | Nombre de réalisations Monte-Carlo. |

### Sorties du script

| Variable | Description |
|---|---|
| `Betti0_all` | \(\beta_0\) pour toutes les réalisations et toutes les distances. |
| `Betti1_graph_all` | \(\beta_1\) du graphe. |
| `Betti1_complex_all` | \(\beta_1\) du complexe de clique. |
| `Betti0`, `Betti1_graph`, `Betti1_complex` | Moyennes Monte-Carlo. |
| `Betti0_std`, `Betti1_graph_std`, `Betti1_complex_std` | Écarts-types Monte-Carlo. |
| `Beta0_theory_sparse` | Approximation sparse. |
| `Beta0_theory_isolated` | Approximation par isolés. |
| `Beta0_theory_dimers_geom` | Approximation incluant les dimères. |
| `Beta0_theory_trimers_geom` | Approximation incluant les trimères. |
| `Beta1_graph_theory_sparse` | Approximation sparse de \(\beta_1^{\mathrm{graphe}}\). |
| `Beta1_graph_theory_isolated` | Approximation par isolés. |
| `Beta1_graph_theory_dimers_geom` | Approximation incluant les dimères. |
| `Beta1_graph_theory_trimers_geom` | Approximation incluant les trimères. |
| `d_alpha_max` | Distance correspondant à la contrainte angulaire. |

### Figures produites

1. \(\beta_0\) moyen en fonction de \(d_{\max}\) ;
2. \(\beta_1^{\mathrm{graphe}}\) moyen en fonction de \(d_{\max}\).

### Fichier sauvegardé

Aucun fichier `.mat` n’est sauvegardé.

---

## `Betti/betti_lambda.m`

### Objectif

Étudie les nombres de Betti en fonction de la densité satellitaire \(\lambda\).

La densité est exprimée en satellites par \(10^6\ \mathrm{km}^2\), puis convertie en satellites par kilomètre carré :

\[
\lambda
=
\frac{\lambda_{\mathrm{scaled}}}{10^6}.
\]

Le nombre de satellites utilisé est déterministe dans ce script :

\[
N
=
\operatorname{round}\left(
\lambda\,4\pi R^2
\right),
\]

et non tiré selon une loi de Poisson.

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
| `lambda_scaled_values` | Densités étudiées. |
| `N_values` | Nombre de satellites associé à chaque densité. |
| `Betti0_all` | Valeurs simulées de \(\beta_0\). |
| `Betti1_graph_all` | Valeurs simulées de \(\beta_1^{\mathrm{graphe}}\). |
| `Betti1_complex_all` | Valeurs simulées de \(\beta_1^{\mathrm{complexe}}\). |
| `Betti0`, `Betti1_graph`, `Betti1_complex` | Moyennes Monte-Carlo. |
| `Betti0_std`, `Betti1_graph_std`, `Betti1_complex_std` | Écarts-types. |
| `Beta0_theory_sparse` | Approximation sparse. |
| `Beta0_theory_isolated` | Approximation par isolés. |
| `Beta0_theory_dimers_geom` | Approximation incluant les dimères. |
| `Beta0_theory_trimers_geom` | Approximation incluant les trimères. |
| `Beta1_graph_theory_sparse` | Approximation sparse de \(\beta_1\). |
| `Beta1_graph_theory_isolated` | Approximation par isolés. |
| `Beta1_graph_theory_dimers_geom` | Approximation incluant les dimères. |
| `Beta1_graph_theory_trimers_geom` | Approximation incluant les trimères. |

### Figures produites

1. \(\beta_0\) moyen en fonction de \(\lambda\) ;
2. \(\beta_1^{\mathrm{graphe}}\) moyen en fonction de \(\lambda\).

### Fichier sauvegardé

Aucun fichier `.mat` n’est sauvegardé.

---

---

# Sous-dossier `Percolation`

## `Percolation/percolation_alpha.m`

### Objectif

Étudie la probabilité de **percolation finie** en fonction de l’angle maximal de connexion \(\alpha_{\max}\), pour un nombre de satellites \(N\) fixé.

L’événement de percolation est défini par

\[
\frac{C_{\max}}{N}\geq \eta,
\]

où \(C_{\max}\) est la taille de la plus grande composante connexe.

Le script compare :

1. l’estimation Monte-Carlo de la probabilité de percolation ;
2. un intervalle de confiance de Hoeffding ;
3. plusieurs bornes analytiques ;
4. une approximation sigmoïde construite à partir du seuil de percolation et du seuil de connexité.

### Type

Script principal utilisant :

- `percolation_alpha_sweep`
- `plot_percolation_curve`

### Entrées du script

| Variable | Type | Description |
|---|---:|---|
| `N` | entier | Nombre de satellites. |
| `numTests` | entier | Nombre de simulations Monte-Carlo par valeur de \(\alpha_{\max}\). |
| `eta` | réel dans \((0,1]\) | Fraction minimale de satellites dans la plus grande composante. |
| `delta` | réel dans \((0,1)\) | Niveau d’erreur de l’intervalle de Hoeffding. |
| `alpha_vals` | vecteur réel | Valeurs de \(\alpha_{\max}\) testées, en radians. |

### Sorties du script

| Variable | Description |
|---|---|
| `res` | Structure contenant l’estimation Monte-Carlo, les bornes et l’approximation par seuils. |
| `alpha_perc` | Seuil angulaire de percolation. |
| `alpha_conn` | Seuil angulaire de connexité. |
| `eps_H` | Demi-largeur de l’intervalle de Hoeffding. |
| `p_perc` | Probabilité de lien correspondant au seuil de percolation. |
| `p_conn` | Probabilité de lien correspondant au seuil de connexité. |

---

## `Percolation/percolation_alpha_sweep.m`

### Objectif

Calcule la probabilité de percolation finie lorsque \(N\) est fixé et que \(\alpha_{\max}\) varie.

### Signature

```matlab
res = percolation_alpha_sweep(N, alpha_vals, numTests, eta, delta)
```

### Entrées

| Argument | Type | Description |
|---|---:|---|
| `N` | entier | Nombre de satellites. |
| `alpha_vals` | vecteur réel | Valeurs testées de l’angle maximal. |
| `numTests` | entier | Nombre de simulations par valeur de l’angle. |
| `eta` | réel | Fraction minimale définissant la percolation. |
| `delta` | réel | Niveau d’erreur utilisé dans Hoeffding. |

### Bornes calculées

La fonction calcule :

#### Intervalle de Hoeffding

\[
\widehat P_{\mathrm{perc}}
\pm
\sqrt{\frac{\log(2/\delta)}{2\,\texttt{numTests}}}.
\]

#### Borne supérieure par satellites non isolés

\[
P(C_{\max}\geq m)
\leq
\frac{\mathbb E[N_{\mathrm{non\,isolés}}]}{m}.
\]

#### Borne supérieure par nombre d’arêtes

\[
P(C_{\max}\geq m)
\leq
\frac{\mathbb E[|E|]}{m-1}.
\]

#### Borne supérieure par arbres couvrants

Une union bound est appliquée sur les sous-ensembles de \(m\) sommets et les arbres couvrants possibles.

#### Borne inférieure par étoile

\[
P(C_{\max}\geq m)
\geq
P\bigl(\deg(1)\geq m-1\bigr).
\]

### Sortie

La fonction retourne une structure `res`.

| Champ | Description |
|---|---|
| `P_perc_hat` | Probabilité empirique de percolation. |
| `Hoeffding_low` | Borne inférieure de Hoeffding. |
| `Hoeffding_up` | Borne supérieure de Hoeffding. |
| `Bound_upper_nonisolated` | Borne supérieure par satellites non isolés. |
| `Bound_upper_edges` | Borne supérieure par nombre d’arêtes. |
| `Bound_upper_tree` | Borne supérieure par arbres couvrants. |
| `Bound_upper_math` | Minimum des bornes supérieures disponibles. |
| `Bound_lower_star` | Borne inférieure par étoile. |
| `P_threshold_approx` | Approximation sigmoïde fondée sur les seuils. |
| `threshold_percolation` | Valeur critique de \(\alpha_{\max}\) pour la percolation. |
| `threshold_connectivity` | Valeur critique de \(\alpha_{\max}\) pour la connexité. |

---

## `Percolation/percolation_dmax.m`

### Objectif

Étudie la probabilité de percolation finie en fonction de la portée maximale \(d_{\max}\), pour \(N\) et \(h\) fixés.

La distance de corde est convertie en angle maximal par

\[
\alpha_{\max}
=
2\arcsin\left(\frac{d_{\max}}{2r}\right),
\qquad
r=R_E+h.
\]

### Type

Script principal utilisant :

- `percolation_alpha_sweep`
- `plot_percolation_curve`

### Entrées du script

| Variable | Type | Description |
|---|---:|---|
| `R_E` | réel | Rayon terrestre en kilomètres. |
| `N` | entier | Nombre de satellites. |
| `h` | réel | Altitude orbitale en kilomètres. |
| `numTests` | entier | Nombre de simulations par point. |
| `eta` | réel | Seuil relatif de percolation. |
| `delta` | réel | Niveau d’erreur de Hoeffding. |
| `dmax_vals` | vecteur réel | Portées maximales testées en kilomètres. |

### Sorties du script

| Variable | Description |
|---|---|
| `res` | Structure renvoyée par `percolation_alpha_sweep`. |
| `dmax_perc` | Portée critique de percolation. |
| `dmax_conn` | Portée critique de connexité. |
| `kbar_perc` | Degré moyen au seuil de percolation. |
| `kbar_conn` | Degré moyen au seuil de connexité. |
| `eps_H` | Demi-largeur de Hoeffding. |

### Figures produites

Une figure de probabilité de percolation en fonction de \(d_{\max}\), enrichie des seuils théoriques.

---

## `Percolation/percolation_h.m`

### Objectif

Étudie la probabilité de percolation finie en fonction de l’altitude \(h\), pour \(N\) et \(d_{\max}\) fixés.

Lorsque \(h\) augmente, le rayon orbital augmente et l’angle couvert par une portée fixe diminue.

### Type

Script principal utilisant :

- `percolation_alpha_sweep`
- `plot_percolation_curve`

### Entrées du script

| Variable | Type | Description |
|---|---:|---|
| `R_E` | réel | Rayon terrestre. |
| `N` | entier | Nombre de satellites. |
| `dmax` | réel | Portée maximale en kilomètres. |
| `numTests` | entier | Nombre de simulations par altitude. |
| `eta` | réel | Seuil relatif de percolation. |
| `delta` | réel | Niveau d’erreur de Hoeffding. |
| `h_vals` | vecteur réel | Altitudes testées. |

### Sorties du script

| Variable | Description |
|---|---|
| `res` | Résultats de la simulation et des bornes. |
| `h_perc` | Altitude critique de percolation. |
| `h_conn` | Altitude critique de connexité. |
| `eps_H` | Demi-largeur de Hoeffding. |

### Figures produites

Une figure de probabilité de percolation en fonction de l’altitude.

---

## `Percolation/percolation_N.m`

### Objectif

Étudie la probabilité de percolation finie en fonction du nombre de satellites \(N\), pour \(h\) et \(d_{\max}\) fixés.

### Type

Script principal utilisant :

- `percolation_N_sweep`
- `plot_percolation_curve`

### Entrées du script

| Variable | Type | Description |
|---|---:|---|
| `R_E` | réel | Rayon terrestre. |
| `h` | réel | Altitude orbitale. |
| `dmax` | réel | Portée de lien. |
| `numTests` | entier | Nombre de simulations par valeur de \(N\). |
| `eta` | réel | Seuil relatif de percolation. |
| `delta` | réel | Niveau d’erreur de Hoeffding. |
| `N_vals` | vecteur d’entiers | Nombres de satellites testés. |
| `useEarthOccultation` | booléen | Paramètre déclaré pour une éventuelle contrainte d’occultation terrestre. |

### Sorties du script

| Variable | Description |
|---|---|
| `res` | Structure des résultats. |
| `N_perc` | Nombre critique de satellites pour la percolation. |
| `N_conn` | Nombre critique de satellites pour la connexité. |
| `p_link` | Probabilité de lien fixe. |
| `eps_H` | Demi-largeur de Hoeffding. |

---

## `Percolation/percolation_N_sweep.m`

### Objectif

Calcule la probabilité de percolation lorsque \(\alpha_{\max}\) est fixé et que le nombre de satellites varie.

### Signature

```matlab
res = percolation_N_sweep(N_vals, alpha_max, numTests, eta, delta)
```

### Entrées

| Argument | Type | Description |
|---|---:|---|
| `N_vals` | vecteur d’entiers | Nombres de satellites testés. |
| `alpha_max` | réel | Angle maximal de connexion. |
| `numTests` | entier | Nombre de simulations par valeur de \(N\). |
| `eta` | réel | Fraction minimale de satellites dans la plus grande composante. |
| `delta` | réel | Niveau d’erreur de Hoeffding. |

### Sortie

La structure `res` contient :

| Champ | Description |
|---|---|
| `P_perc_hat` | Probabilité empirique de percolation. |
| `Hoeffding_low` | Borne inférieure de Hoeffding. |
| `Hoeffding_up` | Borne supérieure de Hoeffding. |
| `Bound_upper_nonisolated` | Borne supérieure par satellites non isolés. |
| `Bound_upper_edges` | Borne supérieure par nombre d’arêtes. |
| `Bound_upper_tree` | Borne supérieure par arbres couvrants. |
| `Bound_upper_math` | Minimum des bornes supérieures. |
| `Bound_lower_star` | Borne inférieure par étoile. |
| `P_threshold_approx` | Approximation sigmoïde. |
| `threshold_percolation` | Seuil critique en nombre de satellites. |
| `threshold_connectivity` | Seuil de connexité en nombre de satellites. |

---

## `Percolation/plot_percolation_curve.m`

### Objectif

Trace de manière standardisée les résultats d’une étude de percolation.

### Signature

```matlab
plot_percolation_curve(x_vals, res, xlab, titre, eta)
```

### Entrées

| Argument | Type | Description |
|---|---:|---|
| `x_vals` | vecteur réel | Valeurs de l’axe horizontal. |
| `res` | structure | Résultats produits par une fonction de sweep. |
| `xlab` | chaîne | Étiquette de l’axe horizontal. |
| `titre` | chaîne | Format du titre MATLAB. |
| `eta` | réel | Seuil relatif de percolation, utilisé dans le titre. |

### Sorties

Cette fonction ne retourne aucune variable.

Elle produit une figure contenant :

- la zone de confiance de Hoeffding ;
- l’estimation Monte-Carlo ;
- l’approximation par seuils si elle existe ;
- trois bornes supérieures ;
- une borne inférieure.

---

---

# Sous-dossier `Routage`

## `Routage/routage_dmax.m`

### Objectif

Étudie la probabilité de routage multi-sauts en fonction de la distance maximale de connexion \(d_{\max}\).

Le script compare :

1. la probabilité obtenue par simulation Monte-Carlo ;
2. une approximation de composante géante de type Erdős-Rényi ;
3. une approximation géométrique corrigée par un seuil de percolation ;
4. la probabilité de lien direct.

La probabilité simulée est calculée à partir des tailles \(s_k\) des composantes connexes :

\[
P_{\mathrm{routing}}
=
\frac{\sum_k s_k(s_k-1)}
{N(N-1)}.
\]

Cette expression est la probabilité que deux satellites distincts tirés uniformément appartiennent à la même composante connexe.

### Type

Script principal avec une fonction locale :

- `giant_component_fraction_ER`

### Entrées du script

| Variable | Type | Description |
|---|---:|---|
| `dmax_values` | vecteur réel | Distances maximales testées, en kilomètres. |
| `nSim` | entier | Nombre de réalisations Monte-Carlo pour chaque valeur de \(d_{\max}\). |
| `R_earth` | réel | Rayon terrestre en kilomètres. |
| `h` | réel | Altitude orbitale en kilomètres. |
| `R` | réel | Rayon orbital. |
| `lambda` | réel | Intensité satellitaire en satellites/km². |
| `k_crit_geo` | réel | Seuil critique effectif du degré moyen pour l’approximation géométrique. |

### Sorties du script

| Variable | Description |
|---|---|
| `P_routing_mean` | Probabilité moyenne de routage obtenue par simulation. |
| `P_routing_std` | Écart-type entre les réalisations. |
| `N_mean` | Nombre moyen empirique de satellites. |
| `P_direct_theory` | Probabilité théorique de lien direct. |
| `P_routing_ER_theory` | Approximation de routage issue du modèle Erdős-Rényi. |
| `P_routing_geo_theory` | Approximation géométrique corrigée. |
| `k_mean_theory` | Degré moyen théorique. |
| `Nbar` | Nombre moyen théorique de satellites. |

---

## `Routage/routage_h.m`

### Objectif

Étudie la probabilité de routage multi-sauts en fonction de l’altitude orbitale \(h\).

La densité satellitaire \(\lambda\) et la portée `dmax` sont fixées. Lorsque l’altitude varie, le rayon orbital et la surface de la sphère changent, ce qui modifie :

- le nombre moyen de satellites ;
- la probabilité de lien direct ;
- le degré moyen ;
- la taille de la composante géante.

### Type

Script principal avec une fonction locale :

- `giant_component_fraction_ER`

### Entrées du script

| Variable | Type | Description |
|---|---:|---|
| `R_earth` | réel | Rayon terrestre en kilomètres. |
| `lambda` | réel | Intensité satellitaire en satellites/km². |
| `dmax` | réel | Distance maximale de connexion en kilomètres. |
| `h_values` | vecteur réel | Altitudes testées. |
| `nSim` | entier | Nombre de simulations par altitude. |
| `k_crit_geo` | réel | Seuil critique utilisé pour la correction géométrique. |

### Sorties du script

| Variable | Description |
|---|---|
| `P_routing_mean` | Probabilité moyenne de routage simulée. |
| `P_routing_std` | Écart-type de la probabilité simulée. |
| `N_mean` | Nombre moyen empirique de satellites. |
| `Nbar_theory` | Nombre moyen théorique de satellites. |
| `P_direct_theory` | Probabilité de lien direct. |
| `P_routing_ER_theory` | Approximation Erdős-Rényi. |
| `P_routing_geo_theory` | Approximation géométrique corrigée. |
| `k_mean_theory` | Degré moyen théorique. |

---

## `Routage/routage_lambda.m`

### Objectif

Étudie la probabilité de routage multi-sauts en fonction de la densité satellitaire \(\lambda\).

Le script compare la simulation à trois grandeurs :

1. la probabilité de lien direct ;
2. l’approximation Erdős-Rényi de la composante géante ;
3. une approximation géométrique corrigée par le seuil
   \[
   k_{\mathrm{crit}}=4.512.
   \]

### Type

Script principal avec une fonction locale :

- `giant_component_fraction_ER`

### Entrées du script

| Variable | Type | Description |
|---|---:|---|
| `R_earth` | réel | Rayon terrestre en kilomètres. |
| `h` | réel | Altitude orbitale en kilomètres. |
| `R` | réel | Rayon orbital. |
| `surface_sphere` | réel | Surface de la sphère orbitale. |
| `dmax` | réel | Distance maximale de connexion. |
| `lambda_values` | vecteur réel | Densités satellitaires testées. |
| `nSim` | entier | Nombre de simulations par densité. |
| `k_crit_geo` | réel | Seuil critique de l’approximation géométrique. |

### Sorties du script

| Variable | Description |
|---|---|
| `P_routing_mean` | Probabilité moyenne de routage simulée. |
| `P_routing_std` | Écart-type des simulations. |
| `N_mean` | Nombre moyen empirique de satellites. |
| `P_direct_theory` | Probabilité de lien direct. |
| `P_routing_ER_theory` | Approximation Erdős-Rényi. |
| `P_routing_geo_theory` | Approximation géométrique corrigée. |
| `k_mean_theory` | Degré moyen théorique. |
| `lambda_values` | Densités satellitaires étudiées. |

---

## `Routage/test_routage.m`

### Objectif

Teste l’existence d’un chemin entre deux satellites distincts choisis au hasard dans un graphe existant.

Le script retourne également le nombre de sauts du plus court chemin :

\[
H=\text{longueur du chemin}-1.
\]

### Type

Script principal sans fonction locale.

### Dépendances

Le script suppose que les variables suivantes sont déjà présentes :

| Variable | Description |
|---|---|
| `G` | Objet MATLAB `graph`. |
| `N` | Nombre de satellites du graphe. |

### Entrées du script

| Variable | Type | Description |
|---|---:|---|
| `G` | objet `graph` | Graphe de la constellation. |
| `N` | entier | Nombre de sommets. |

### Sorties du script

| Variable | Description |
|---|---|
| `pair` | Paire de satellites tirée. |
| `s` | Satellite source. |
| `t` | Satellite destination. |
| `path` | Plus court chemin entre les deux satellites. |
| `nb_hops` | Nombre de sauts du chemin. |

---

---

# Sous-dossier `Valeurs moyennes`

## `Valeurs moyennes/mean_degre.m`

### Objectif

Vérifie la formule théorique du **degré moyen** dans un graphe géométrique aléatoire construit à partir de satellites uniformément distribués sur la sphère unité.

Deux satellites sont reliés lorsque leur séparation angulaire \(\alpha\) vérifie

\[
\alpha \leq \alpha_{\max}.
\]

Pour deux points indépendants uniformes sur la sphère, la probabilité de lien vaut

\[
p_{\mathrm{link}}
=
\frac{1-\cos(\alpha_{\max})}{2}.
\]

Le degré théorique moyen est alors

\[
\mathbb E[\deg]
=
(N-1)p_{\mathrm{link}}.
\]

Le script compare cette expression à une estimation obtenue par simulation Monte-Carlo.

### Type

Script principal sans fonction locale.

### Entrées du script

| Variable | Type | Description |
|---|---:|---|
| `N` | entier | Nombre de satellites dans chaque réalisation. |
| `numTests` | entier | Nombre de réalisations Monte-Carlo effectuées pour chaque valeur de \(\alpha_{\max}\). |
| `alpha_vals` | vecteur réel | Valeurs testées de l’angle maximal de connexion, en radians. |

### Sorties du script

| Variable | Description |
|---|---|
| `E_deg_sim` | Estimation Monte-Carlo du degré moyen pour chaque valeur de `alpha_vals`. |
| `E_deg_theo` | Valeur théorique du degré moyen. |
| `alpha_vals` | Valeurs de \(\alpha_{\max}\) étudiées. |
| `deg_mean_tests` | Degrés moyens obtenus dans les réalisations pour la valeur courante de \(\alpha_{\max}\). |
| `positions` | Dernier ensemble de positions généré. |
| `A` | Dernière matrice d’adjacence construite. |
| `deg` | Degrés des satellites dans la dernière réalisation. |

---

## `Valeurs moyennes/mean_edges.m`

### Objectif

Vérifie la formule théorique du **nombre moyen d’arêtes** dans un graphe de satellites uniformément distribués sur la sphère unité.

Pour une paire de satellites,

\[
p_{\mathrm{link}}
=
\frac{1-\cos(\alpha_{\max})}{2}.
\]

Comme il existe

\[
\binom{N}{2}
\]

paires distinctes, le nombre moyen théorique d’arêtes est

\[
\mathbb E[|E|]
=
\binom{N}{2}p_{\mathrm{link}}.
\]

### Type

Script principal sans fonction locale.

### Entrées du script

| Variable | Type | Description |
|---|---:|---|
| `N` | entier | Nombre de satellites par réalisation. |
| `numTests` | entier | Nombre de réalisations Monte-Carlo pour chaque angle maximal. |
| `alpha_vals` | vecteur réel | Valeurs de \(\alpha_{\max}\) testées, en radians. |

### Sorties du script

| Variable | Description |
|---|---|
| `E_edges_sim` | Nombre moyen d’arêtes estimé par Monte-Carlo. |
| `E_edges_theo` | Nombre moyen théorique d’arêtes. |
| `alpha_vals` | Valeurs de \(\alpha_{\max}\) étudiées. |
| `edges_tests` | Nombre d’arêtes de chaque réalisation pour l’angle courant. |
| `positions` | Dernières positions générées. |
| `A` | Dernière matrice d’adjacence. |
| `nb_edges` | Nombre d’arêtes dans la dernière réalisation. |

---

## `Valeurs moyennes/mean_plink.m`

### Objectif

Vérifie la formule théorique de la **probabilité moyenne de lien** entre deux satellites uniformément distribués sur la sphère unité.

La probabilité théorique est

\[
p_{\mathrm{link}}
=
\mathbb P(\alpha\leq\alpha_{\max})
=
\frac{1-\cos(\alpha_{\max})}{2}.
\]

Cette expression correspond au rapport entre l’aire d’une calotte sphérique d’angle \(\alpha_{\max}\) et l’aire totale de la sphère.

### Type

Script principal sans fonction locale.

### Entrées du script

| Variable | Type | Description |
|---|---:|---|
| `numTests` | entier | Nombre de paires de satellites générées pour chaque valeur de \(\alpha_{\max}\). |
| `alpha_vals` | vecteur réel | Valeurs de l’angle maximal de connexion, en radians. |

### Sorties du script

| Variable | Description |
|---|---|
| `P_sim` | Probabilité de lien estimée par Monte-Carlo. |
| `P_theo` | Probabilité de lien théorique. |
| `alpha_vals` | Valeurs de \(\alpha_{\max}\) étudiées. |
| `alpha` | Angles calculés pour les paires de la dernière itération. |
| `dotProduct` | Produits scalaires des dernières paires générées. |
| `x1`, `y1`, `z1` | Coordonnées du premier ensemble de satellites. |
| `x2`, `y2`, `z2` | Coordonnées du second ensemble de satellites. |
