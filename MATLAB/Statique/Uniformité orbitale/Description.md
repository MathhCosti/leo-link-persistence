# Description des scripts et fonctions MATLAB

Ce document décrit les fichiers MATLAB, ainsi que leurs fonctions locales.

---

---

## `graphe_3D_delta.m`

### Objectif

Génère un **graphe statique tridimensionnel d’une constellation LEO Walker Delta stochastique** à l’instant $t$, construit les liens intersatellites respectant une portée maximale, puis affiche la constellation et ses arêtes dans une figure 3D.

Le nombre de satellites est tiré selon une loi de Poisson :

$$N\sim\mathcal P\!\left(\lambda 4\pi R^2\right),\qquad R=R_{\mathrm{Earth}}+h.$$

Chaque satellite possède un RAAN $\Omega$ et une phase initiale $u_0$ uniformes sur $[0,2\pi)$, puis

$$u(t)=u_0+\omega t\pmod{2\pi},\qquad \omega=\sqrt{\mu/R^3}.$$

Deux satellites sont reliés lorsque $\|\mathbf r_i-\mathbf r_j\|\le d_{\max}$.

### Type

Fonction MATLAB principale avec une fonction locale : `walker_delta_positions`.

### Signature

```matlab
[G, positions, A, params] = ...
    graphe_3D_delta(lambda, inc_deg, dmax, t, h, rng_seed)
```

Tous les arguments sont optionnels.

### Entrées

| Argument | Type | Unité | Valeur par défaut | Description |
|---|---:|---:|---:|---|
| `lambda` | réel positif | satellites/km² | `4e-7` | Densité satellitaire sur la sphère orbitale. |
| `inc_deg` | réel dans $[0,90]$ | degrés | `58` | Inclinaison orbitale commune. |
| `dmax` | réel positif | km | `2500` dans le code | Distance maximale d’un lien intersatellite. |
| `t` | réel positif ou nul | s | `0` | Instant auquel la topologie est construite. |
| `h` | réel positif | km | `550` | Altitude orbitale. |
| `rng_seed` | entier ou `[]` | — | `[]` | Graine aléatoire ; `[]` conserve l’état courant. |

### Sorties

| Sortie | Type ou dimensions | Description |
|---|---:|---|
| `G` | objet MATLAB `graph` | Graphe non orienté de la constellation. |
| `positions` | matrice `N x 3` | Positions cartésiennes des satellites en kilomètres. |
| `A` | matrice creuse `N x N` | Matrice d’adjacence symétrique. |
| `params` | structure MATLAB | Paramètres physiques, tirages orbitaux et informations sur le graphe. |

### Contenu de `params`

`R_earth`, `h`, `R`, `mu`, `omega`, `lambda`, `N_mean`, `N`, `inc_deg`, `inc`, `dmax`, `t`, `Omega`, `u0`, `u_t`, `num_edges`.

---

---

# Sous-dossier `Valeurs moyennes`

## `Valeurs moyennes/mean_degre.m`

### Objectif

Calcule et compare le **degré moyen** d’un graphe de satellites Walker Delta à uniformité orbitale en fonction de l’angle maximal de connexion \(\alpha_{\max}\).

Le script compare :

1. une estimation empirique par simulation de Monte-Carlo ;
2. le modèle semi-analytique Walker Delta obtenu par quadrature ;
3. la référence correspondant à des satellites uniformément distribués sur toute la sphère.

La formule théorique utilisée est

\[
\mathbb E[\deg]=(N-1)p_{\mathrm{link}}^\Delta.
\]

### Type

Script principal avec trois fonctions locales :

- `sample_walker_delta`
- `plink_delta_quadrature`
- `gauss_legendre_interval`

### Entrées du script

| Variable | Type | Description |
|---|---:|---|
| `N` | entier | Nombre de satellites générés dans chaque réalisation. |
| `numTests` | entier | Nombre de réalisations Monte-Carlo pour chaque valeur de \(\alpha_{\max}\). |
| `alpha_vals` | vecteur réel | Valeurs de l’angle géocentrique maximal de connexion, en radians. |
| `inc_deg` | réel | Inclinaison orbitale en degrés. |
| `inc` | réel | Inclinaison orbitale en radians. |
| `nQuad` | entier | Ordre de la quadrature de Gauss-Legendre dans chaque dimension. |

### Sorties du script

| Variable | Description |
|---|---|
| `E_deg_sim` | Degré moyen empirique estimé par Monte-Carlo. |
| `E_deg_theo` | Degré moyen théorique Walker Delta. |
| `E_deg_sphere` | Degré moyen sous l’hypothèse uniforme sur la sphère. |
| `alpha_vals` | Valeurs de \(\alpha_{\max}\) utilisées. |

Le script produit également deux figures et sauvegarde :

```matlab
prob_degre_walker_delta_results.mat
```

avec :

```matlab
N, alpha_vals, inc_deg, nQuad, ...
E_deg_sim, E_deg_theo, E_deg_sphere
```

---

## `Valeurs moyennes/mean_edges.m`

### Objectif

Calcule et compare le **nombre moyen total d’arêtes** du graphe Walker Delta en fonction de l’angle maximal de connexion.

La formule utilisée est

\[
\mathbb E[|E|]
=
\binom{N}{2}p_{\mathrm{link}}^\Delta.
\]

### Type

Script principal avec trois fonctions locales :

- `sample_walker_delta`
- `plink_delta_quadrature`
- `gauss_legendre_interval`

### Entrées du script

| Variable | Type | Description |
|---|---:|---|
| `N` | entier | Nombre de satellites. |
| `numTests` | entier | Nombre de simulations Monte-Carlo pour chaque valeur de \(\alpha_{\max}\). |
| `alpha_vals` | vecteur réel | Valeurs de l’angle maximal de connexion, en radians. |
| `inc_deg` | réel | Inclinaison en degrés. |
| `inc` | réel | Inclinaison en radians. |
| `nQuad` | entier | Ordre de quadrature par dimension. |

### Sorties du script

| Variable | Description |
|---|---|
| `E_edges_sim` | Nombre moyen d’arêtes obtenu par simulation. |
| `E_edges_theo` | Nombre moyen théorique d’arêtes Walker Delta. |
| `E_edges_sphere` | Nombre moyen d’arêtes sous l’hypothèse uniforme sur la sphère. |
| `alpha_vals` | Angles maximaux étudiés. |

Le script produit deux figures et sauvegarde :

```matlab
prob_edges_walker_delta_results.mat
```

avec :

```matlab
N, alpha_vals, inc_deg, nQuad, ...
E_edges_sim, E_edges_theo, E_edges_sphere
```

---

## `Valeurs moyennes/mean_plink.m`

### Objectif

Compare la **probabilité moyenne globale de lien** entre deux satellites Walker Delta obtenue par :

1. simulation Monte-Carlo ;
2. quadrature semi-analytique ;
3. modèle uniforme sur la sphère.

La grandeur calculée est

\[
p_{\mathrm{link}}^\Delta
=
\mathbb P(\alpha\leq\alpha_{\max}),
\]

où \(\alpha\) est l’angle géocentrique séparant deux satellites tirés indépendamment selon l’uniformité orbitale Walker Delta.

### Type

Script principal avec trois fonctions locales :

- `sample_walker_delta`
- `plink_delta_quadrature`
- `gauss_legendre_interval`

### Entrées du script

| Variable | Type | Description |
|---|---:|---|
| `numTests` | entier | Nombre de paires de satellites tirées pour chaque valeur de \(\alpha_{\max}\). |
| `alpha_vals` | vecteur réel | Valeurs de l’angle maximal de connexion. |
| `inc_deg` | réel | Inclinaison orbitale en degrés. |
| `inc` | réel | Inclinaison en radians. |
| `nQuad` | entier | Ordre de quadrature par dimension. |

### Sorties du script

| Variable | Description |
|---|---|
| `P_sim` | Probabilité de lien estimée par Monte-Carlo. |
| `P_theo` | Probabilité de lien Walker Delta calculée par quadrature. |
| `P_sphere` | Probabilité de lien pour une distribution uniforme sur la sphère. |
| `alpha_vals` | Angles maximaux étudiés. |

Le script produit deux figures et sauvegarde :

```matlab
prob_lien_walker_delta_results.mat
```

avec :

```matlab
alpha_vals, inc_deg, nQuad, P_sim, P_theo, P_sphere
```

---

## `Valeurs moyennes/plink_phi.m`

### Objectif

Calcule la **probabilité locale de lien conditionnée par la latitude** d’un premier satellite :

\[
p_{\mathrm{link}}(\phi)
=
\mathbb P(1\leftrightarrow2\mid\Phi_1=\phi).
\]

Le script calcule également la densité orbitale de latitude \(f_\Phi(\phi)\) et reconstruit la probabilité globale par

\[
p_{\mathrm{link}}^\Delta
=
\int_{-i}^{i}
p_{\mathrm{link}}(\phi)f_\Phi(\phi)\,d\phi.
\]

### Type

Script principal avec une fonction locale :

- `gauss_legendre_interval`

### Entrées du script

| Variable | Type | Description |
|---|---:|---|
| `R` | réel | Rayon orbital en kilomètres. |
| `inc_deg` | réel | Inclinaison en degrés. |
| `inc` | réel | Inclinaison en radians. |
| `dmax` | réel | Distance maximale de connexion, en kilomètres. |
| `n_phi` | entier | Nombre de latitudes auxquelles \(p_{\mathrm{link}}(\phi)\) est évaluée. |
| `nQuad` | entier | Ordre de quadrature sur la phase orbitale du second satellite. |
| `eps_phi` | réel | Décalage évitant les singularités aux latitudes limites \(\pm i\). |

### Sorties du script

| Variable | Description |
|---|---|
| `phi_vals` | Latitudes auxquelles la probabilité locale est évaluée. |
| `p_link_phi` | Probabilité locale conditionnelle \(p_{\mathrm{link}}(\phi)\). |
| `f_phi` | Densité de latitude orbitale normalisée. |
| `p_link_global_from_phi` | Probabilité globale reconstruite par intégration sur la latitude. |
| `p_link_sphere` | Référence uniforme sur la sphère. |
| `alpha_max` | Angle géocentrique maximal correspondant à `dmax`. |

Le script produit une figure, affiche les résultats principaux dans la console et sauvegarde :

```matlab
plink_phi_walker_delta_results.mat
```

avec :

```matlab
R, inc_deg, inc, dmax, alpha_max, phi_vals, ...
p_link_phi, f_phi, p_link_global_from_phi, ...
p_link_sphere, nQuad
```

---

## `Valeurs moyennes/plink_quadrature.m`

### Objectif

Calcule la probabilité globale de lien du modèle Walker Delta par une **réduction semi-analytique**.

L’intégrale sur

\[
\Delta\Omega=\Omega_2-\Omega_1
\]

est résolue analytiquement. Il reste une quadrature déterministe bidimensionnelle sur \((u_1,u_2)\).

Le script teste plusieurs ordres de quadrature afin de vérifier la convergence.

### Type

Script principal avec une fonction locale :

- `gauss_legendre_interval`

### Entrées du script

| Variable | Type | Description |
|---|---:|---|
| `R_earth` | réel | Rayon terrestre, en kilomètres. |
| `h` | réel | Altitude orbitale, en kilomètres. |
| `R` | réel | Rayon orbital \(R=R_{\mathrm{earth}}+h\). |
| `inc_deg` | réel | Inclinaison orbitale en degrés. |
| `inc` | réel | Inclinaison en radians. |
| `dmax` | réel | Distance maximale de connexion, en kilomètres. |
| `quad_orders` | vecteur d’entiers | Ordres de quadrature testés. |
| `N` | entier | Nombre de satellites utilisé dans l’exemple final. |

### Sorties du script

| Variable | Description |
|---|---|
| `p_quad` | Valeur de \(p_{\mathrm{link}}^\Delta\) pour chaque ordre de quadrature. |
| `p_link_delta` | Valeur finale retenue. |
| `p_link_sphere` | Probabilité uniforme sur la sphère. |
| `alpha_max` | Angle maximal de connexion. |
| `cmax` | Cosinus du seuil angulaire. |
| `E_delta` | Nombre moyen d’arêtes Walker Delta pour l’exemple. |
| `k_delta` | Degré moyen Walker Delta pour l’exemple. |
| `E_sphere` | Nombre moyen d’arêtes sous l’hypothèse sphérique. |
| `k_sphere` | Degré moyen sous l’hypothèse sphérique. |

Le script produit une figure de convergence, affiche les résultats dans la console et sauvegarde :

```matlab
plink_walker_delta_semi_analytique_results.mat
```

avec :

```matlab
R, h, inc_deg, inc, dmax, alpha_max, cmax, ...
quad_orders, p_quad, p_link_delta, p_link_sphere, ...
N, E_delta, k_delta, E_sphere, k_sphere
```

---

---

# Sous-dossier `Betti`

## `Betti/betti_quadrature.m`

### Objectif

Estime les nombres moyens de petites composantes connexes d’un graphe Walker Delta :

- \(N_1\) : nombre de satellites isolés ;
- \(N_2\) : nombre de dimères isolés ;
- \(N_3\) : nombre de trimères isolés.

Le script en déduit ensuite une approximation du nombre de Betti d’ordre zéro :

\[
\mathbb E[\beta_0]
\approx
C_{\mathrm{macro}}
+
\mathbb E[N_1]
+
\mathbb E[N_2]
+
\mathbb E[N_3],
\]

avec, dans le code,

\[
C_{\mathrm{macro}}=2.
\]

Cette constante représente les composantes connexes principales supposées.

Le script compare les estimations obtenues par quadrature déterministe à des graphes statiques générés par simulation Monte-Carlo.

### Type

Script principal avec trois fonctions locales :

- `walker_delta_position`
- `sample_walker_delta_unit`
- `gauss_legendre_interval`

### Entrées du script

#### Paramètres physiques

| Variable | Type | Description |
|---|---:|---|
| `R` | réel | Rayon orbital en kilomètres. |
| `inc_deg` | réel | Inclinaison orbitale en degrés. |
| `inc` | réel | Inclinaison orbitale en radians. |
| `dmax` | réel | Distance euclidienne maximale de connexion, en kilomètres. |
| `N` | entier | Nombre total de satellites. |

#### Paramètres numériques

| Variable | Type | Description |
|---|---:|---|
| `n_phi` | entier | Nombre de latitudes utilisées pour représenter \(p_{\mathrm{link}}(\phi)\). |
| `nQuad` | entier | Ordre de quadrature utilisé pour calculer la probabilité locale de lien. |
| `nQuad_N1` | entier | Ordre de quadrature pour l’intégration finale de \(N_1\) sur la phase \(u_1\). |
| `nQuad_N2_outer` | entier | Ordre de la quadrature extérieure à trois dimensions pour \(N_2\). |
| `nQuad_N2_probe` | entier | Ordre de la quadrature servant à évaluer la probabilité qu’un satellite extérieur soit relié au dimère. |
| `nQuad_N3_outer` | entier | Ordre de la quadrature extérieure à cinq dimensions pour \(N_3\). |
| `nQuad_N3_probe` | entier | Ordre de la quadrature servant à évaluer la probabilité qu’un satellite extérieur soit relié au trimère. |
| `n_graph_realizations` | entier | Nombre de graphes statiques simulés pour valider les prédictions. |

### Sorties du script

#### Résultats théoriques

| Variable | Description |
|---|---|
| `N1_th_local` | Nombre moyen théorique de satellites isolés, calculé avec la probabilité de lien locale. |
| `N2_th_geom` | Nombre moyen théorique de dimères isolés. |
| `N3_th_geom` | Nombre moyen théorique de trimères isolés. |
| `C_macro` | Nombre supposé de composantes macroscopiques, fixé à 1. |
| `beta0_th_123` | Approximation de \(\mathbb E[\beta_0]\) par \(C_{\mathrm{macro}}+N_1+N_2+N_3\). |
| `p_link_phi` | Probabilité locale de lien représentée en fonction de la latitude. |
| `phi_vals` | Latitudes associées à `p_link_phi`. |

#### Résultats empiriques

| Variable | Description |
|---|---|
| `N1_emp` | Nombre de composantes de taille 1 dans chaque graphe simulé. |
| `N2_emp` | Nombre de composantes de taille 2 dans chaque graphe simulé. |
| `N3_emp` | Nombre de composantes de taille 3 dans chaque graphe simulé. |
| `beta0_emp` | Nombre total de composantes connexes dans chaque graphe simulé. |

#### Figures

Le script produit trois figures :

1. comparaison théorie/simulation pour \(N_1\), \(N_2\) et \(N_3\) ;
2. histogramme empirique de \(\beta_0\) avec l’approximation théorique ;
3. probabilité locale de lien \(p_{\mathrm{link}}(\phi)\).

#### Fichier sauvegardé

```matlab
N1_N2_N3_walker_delta_results.mat
```

avec notamment :

```matlab
R, inc_deg, inc, dmax, N, alpha_max, ...
phi_vals, p_link_phi, N1_th_local, ...
N2_th_geom, N3_th_geom, C_macro, beta0_th_123, ...
N1_emp, N2_emp, N3_emp, beta0_emp, ...
nQuad_N2_outer, nQuad_N2_probe, ...
nQuad_N3_outer, nQuad_N3_probe, ...
n_graph_realizations
```

---

## `Betti/isolated_phi.m`

### Objectif

Étudie les satellites isolés en tenant explicitement compte de la variation de la probabilité de lien avec la latitude.

Le script compare :

1. le nombre empirique de satellites isolés ;
2. une prédiction locale utilisant \(p_{\mathrm{link}}(\phi)\) ;
3. une prédiction globale utilisant une unique probabilité moyenne de lien.

La prédiction locale est

\[
\mathbb E[N_1]_{\mathrm{local}}
=
N\sum_b
w_b
\left(1-p_{\mathrm{link},b}\right)^{N-1},
\]

où les bandes de latitude sont indexées par \(b\).

La prédiction globale est

\[
\mathbb E[N_1]_{\mathrm{global}}
=
N
\left(1-\overline{p}_{\mathrm{link}}\right)^{N-1}.
\]

Le script permet donc de mesurer l’erreur introduite lorsque la moyenne spatiale est effectuée avant l’application de la fonction non linéaire.

### Type

Script principal sans fonction locale définie dans le fichier.

### Dépendance externe

Le script appelle la fonction suivante, qui doit être accessible dans le chemin MATLAB :

```matlab
walker_delta_static_sample(N, R, inc_rad)
```

Cette fonction doit retourner au minimum :

| Sortie utilisée | Description |
|---|---|
| `positions` | Positions cartésiennes des `N` satellites. |
| `lat` | Latitude de chaque satellite. |

Les deuxième et troisième sorties sont ignorées dans ce script avec la syntaxe `~`.

### Entrées du script

| Variable | Type | Description |
|---|---:|---|
| `R` | réel | Rayon orbital en kilomètres. |
| `inc_deg` | réel | Inclinaison orbitale en degrés. |
| `inc_rad` | réel | Inclinaison orbitale en radians. |
| `N` | entier | Nombre de satellites par réalisation. |
| `dmax` | réel | Distance maximale de connexion en kilomètres. |
| `n_realizations` | entier | Nombre de graphes statiques simulés. |
| `nbins` | entier | Nombre de bandes de latitude utilisées pour l’estimation locale. |
| `lat_edges` | vecteur réel | Bornes des bandes de latitude. |
| `lat_centers` | vecteur réel | Centres des bandes de latitude. |

### Sorties du script

#### Variables locales par latitude

| Variable | Description |
|---|---|
| `nodes_per_bin` | Nombre total de satellites observés dans chaque bande. |
| `degree_sum_bin` | Somme des degrés observés dans chaque bande. |
| `isolated_sum_bin` | Nombre total de satellites isolés observés dans chaque bande. |
| `mean_degree_lat` | Degré moyen conditionnel à la bande de latitude. |
| `p_link_lat` | Probabilité locale de lien estimée à partir du degré moyen. |
| `p_isolated_lat_emp` | Probabilité empirique d’être isolé dans chaque bande. |
| `p_isolated_lat_th` | Probabilité théorique locale d’être isolé. |
| `weights_lat` | Poids empirique de chaque bande dans la distribution de latitude. |

#### Résultats globaux

| Variable | Description |
|---|---|
| `isolated_emp` | Nombre de satellites isolés dans chaque réalisation. |
| `N1_th_local` | Prédiction locale du nombre moyen de satellites isolés. |
| `p_link_global` | Probabilité globale de lien obtenue par moyenne pondérée des probabilités locales. |
| `N1_th_global` | Prédiction globale utilisant une seule probabilité moyenne de lien. |

#### Figures

Le script produit trois figures :

1. probabilité de lien et degré moyen en fonction de la latitude ;
2. comparaison entre probabilité empirique et prédiction locale d’isolement ;
3. histogramme du nombre empirique de satellites isolés, avec les prédictions locale et globale.

#### Fichier sauvegardé

```matlab
verify_local_plink_isolated_results.mat
```

L’appel à `save` ne fournit pas de liste de variables : le fichier contient donc toutes les variables présentes dans le workspace au moment de la sauvegarde.