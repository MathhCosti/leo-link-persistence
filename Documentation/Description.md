# Description des fichiers

## `0812.0197v1.pdf`

**Titre :** *Zigzag Persistence* — Gunnar Carlsson et Vin de Silva.

Article fondateur sur la persistance zigzag, une généralisation de l’homologie persistante classique. Alors que la persistance standard suppose une suite d’inclusions toutes orientées dans le même sens, la persistance zigzag autorise des flèches dans les deux directions entre espaces successifs. Le document développe les fondements algébriques, la décomposition en intervalles, les barcodes associés et des algorithmes de calcul.

**Intérêt pour le projet :** ce texte fournit la base théorique des barcodes utilisés pour suivre la naissance, la disparition et la réapparition de composantes connexes dans une topologie satellitaire variable au cours du temps.

---

## `1704.06400v3.pdf`

**Titre :** *Counting k-Hop Paths in the Random Connection Model* — Alexander P. Kartun-Giles et Sunwoo Kim.

L’article étudie le nombre de chemins de exactement \(k\) sauts entre deux nœuds d’un graphe géométrique aléatoire. Les auteurs dérivent une expression générale du nombre moyen de chemins, puis des résultats fermés pour un modèle de connexion avec évanouissement de Rayleigh. Ils calculent également la variance du nombre de chemins à trois sauts et relient les moments factoriels à la probabilité d’existence d’au moins un chemin.

**Intérêt pour le projet :** il fournit une base théorique pour relier la distance géométrique, la densité de satellites, le nombre de sauts et la probabilité d’existence d’un chemin de routage.

---

## `A_distributed_routing_algorithm_for_datagram_traffic_in_LEO_satellite_networks.pdf`

**Titre :** *A Distributed Routing Algorithm for Datagram Traffic in LEO Satellite Networks* — Eylem Ekici, Ian F. Akyildiz et Michael D. Bender.

Cet article propose un algorithme distribué de routage par datagrammes pour les réseaux satellitaires LEO. Les décisions sont prises localement pour chaque paquet, sans échange global d’informations de topologie, en s’appuyant sur des positions logiques fixes occupées successivement par les satellites. L’objectif principal est de construire des chemins minimisant le délai de propagation malgré les variations des liens interplans et les coupures dans les régions polaires.

**Intérêt pour le projet :** il constitue une référence classique sur le routage géométrique distribué, les liens intra- et interplans, les coutures orbitales et l’impact de la dynamique des satellites sur les chemins.

---

## `A_dynamic_routing_concept_for_ATM-based_satellite_personal_communication_networks.pdf`

**Titre :** *A Dynamic Routing Concept for ATM-Based Satellite Personal Communication Networks* — Markus Werner.

Le document introduit le schéma **DT-DVTR** (*Discrete-Time Dynamic Virtual Topology Routing*) pour des réseaux LEO à topologie périodiquement variable. La période orbitale est découpée en intervalles temporels, une topologie virtuelle est calculée pour chaque intervalle, puis des séquences de chemins sont optimisées hors ligne. L’article traite notamment du changement de chemin, du délai, de la gigue et de l’implémentation dans un réseau ATM.

**Intérêt pour le projet :** cette référence historique formalise l’approche par snapshots et montre comment exploiter la périodicité orbitale pour anticiper les changements de routage.

---

## `A_survey_of_Mobility_management_for_mobile_networks_supporting_LEO_satellite_access.pdf`

**Titre :** *A Survey of Mobility Management for Mobile Networks Supporting LEO Satellite Access* — Peirong Xie, Qingyang Wang et Jie Chen.

Cette revue présente les mécanismes de gestion de mobilité applicables aux réseaux mobiles intégrant un accès LEO. Elle compare différentes approches IP et 5G, notamment MIPv6, HMIPv6 et les architectures NTN transparentes ou régénératives. Le document analyse les difficultés liées aux changements fréquents de satellite, aux mises à jour de localisation et à la signalisation.

**Intérêt pour le projet :** il apporte le contexte réseau nécessaire pour comprendre les handovers utilisateur–satellite, la continuité de service et l’architecture d’accès 5G/6G NTN.

---

## `An_operational_and_performance_overview_of_the_IRIDIUM_low_earth_orbit_satellite_system.pdf`

**Titre :** *An Operational and Performance Overview of the IRIDIUM Low Earth Orbit Satellite System* — Stephen R. Pratt, Richard A. Raines, Carl E. Fossa Jr. et Michael A. Temple.

Article de synthèse sur la constellation Iridium historique : 66 satellites répartis sur six plans quasi polaires à environ 780 km d’altitude. Il décrit les empreintes de couverture, les liens intersatellites, les liens permanents intraplans, les liens interplans dynamiques et les difficultés de routage autour des pôles et de la couture entre plans contrarotatifs. Des résultats de simulation évaluent également les performances et la robustesse du système.

**Intérêt pour le projet :** il fournit un exemple réel de constellation Walker Star et des paramètres concrets pour valider les hypothèses de visibilité, de durée de lien et de topologie dynamique.

---

## `Connectivity_of_LEO_Satellite_Mega_Constellations_An_Application_of_Percolation_Theory_on_a_Sphere.pdf`

**Titre :** *Connectivity of LEO Satellite Mega Constellations: An Application of Percolation Theory on a Sphere* — Hao Lin, Mustafa A. Kishk et Mohamed-Slim Alouini.

L’article applique la théorie de la percolation aux zones de couverture de satellites distribués sur une sphère. Il définit la percolation sur la surface terrestre, utilise une projection stéréographique pour relier le problème sphérique à un problème plan et dérive des seuils critiques portant sur le nombre de satellites, l’altitude et la portée maximale. L’objectif est de déterminer les conditions d’apparition d’une grande zone de couverture continue.

**Intérêt pour le projet :** il est directement pertinent pour les seuils de connectivité sur une sphère, les modèles ponctuels aléatoires et la transition entre régimes sous-critique et supercritique.

---

## `Continuum Percolation Thresholds - 1209.4936v2.pdf`

**Titre :** *Continuum Percolation Thresholds in Two Dimensions* — Stephan Mertens et Cristopher Moore.

Cet article calcule avec grande précision les seuils de percolation continue en deux dimensions pour des disques, carrés et bâtonnets. Les auteurs adaptent l’algorithme union-find de Newman–Ziff au cas continu et utilisent des conditions périodiques pour détecter les composantes qui traversent ou entourent le domaine. Le texte fournit notamment des valeurs critiques de facteur de remplissage et discute les effets de taille finie.

**Intérêt pour le projet :** il sert de référence pour le seuil critique du modèle de disques de Gilbert et pour la comparaison entre résultats simulés et seuils théoriques de connectivité.

---

## `Dynamic_Topology_Optimization_of_Mega-Constellation_Satellite_Networks_Based_on_Morphological_Computing.pdf`

**Titre :** *Dynamic Topology Optimization of Mega-Constellation Satellite Networks Based on Morphological Computing* — Yufeng Zhang et al.

Le papier propose une méthode d’optimisation de topologie appelée **MCDTO**, fondée sur le calcul morphologique. La topologie continue est décomposée en unités spatio-temporelles relativement stables, puis les liens sont évalués à partir de plusieurs indicateurs : géométrie, variation de distance, charge et gigue. Cette décomposition vise à réduire la complexité de gestion des méga-constellations et à améliorer la stabilité des liens critiques.

**Intérêt pour le projet :** il offre une approche récente de segmentation temporelle et de quantification de la stabilité des liens, utile pour comparer avec un modèle basé sur la persistance ou les probabilités de rupture.

---

## `Exploiting_topology_awareness_for_routing_in_LEO_satellite_constellations.pdf`

**Titre :** *Exploiting Topology Awareness for Routing in LEO Satellite Constellations* — Jonas W. Rabjerg, Israel Leyva-Mayorga, Beatriz Soret et Petar Popovski.

L’article étudie le routage unipath entre stations au sol à travers une constellation Walker Star. Il propose une métrique de routage consciente de la topologie, fondée sur les pertes de propagation, afin de privilégier les liens intersatellites à haut débit. Les performances sont comparées à des métriques basées sur le nombre de sauts et la latence, en prenant en compte propagation, transmission et files d’attente.

**Intérêt pour le projet :** il permet de relier la géométrie des chemins à leurs performances réelles et montre qu’un chemin minimal en nombre de sauts n’est pas nécessairement le meilleur en débit ou en délai.

---

## `Research_of_Adaptive_Routing_Scheme_for_LEO_Network.pdf`

**Titre :** *Research of Adaptive Routing Scheme for LEO Network* — Yuancao Lv et al.

Ce papier propose une stratégie de routage hybride combinant un mode statique par snapshots et un mode dynamique inspiré d’OSPF. En fonctionnement normal, les satellites utilisent des tables précalculées ; lorsqu’une panne ou une modification inattendue de topologie est détectée, le réseau bascule temporairement vers un calcul dynamique. Les auteurs cherchent ainsi à conserver la faible complexité du routage statique tout en améliorant sa robustesse.

**Intérêt pour le projet :** il est utile pour comprendre les compromis entre stabilité des chemins, charge de signalisation, fréquence des snapshots et réaction aux ruptures de liens.

---

## `Time-Varying_Topology_Model_for_Dynamic_Routing_in_LEO_Satellite_Constellation_Networks.pdf`

**Titre :** *Time-Varying Topology Model for Dynamic Routing in LEO Satellite Constellation Networks* — Zhenzhen Han et al.

L’article propose un graphe d’évolution espace-temps pondéré pour représenter une topologie LEO continuellement variable, plutôt qu’une simple succession de snapshots binaires. Les poids des liens combinent plusieurs attributs, notamment le SNR, la durée de disponibilité et l’état des files d’attente. À partir de ce modèle, les auteurs construisent l’algorithme **IUDR**, qui sélectionne les chemins en fonction de l’utilité globale des liens.

**Intérêt pour le projet :** ce document est particulièrement pertinent pour la modélisation de la durée de vie des liens et pour l’intégration de la stabilité temporelle dans une métrique de routage.
