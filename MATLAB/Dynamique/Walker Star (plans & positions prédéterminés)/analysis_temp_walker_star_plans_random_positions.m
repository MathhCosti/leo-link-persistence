clear; clc; close all;

%% ============================================================
%  ÉTUDE TOPOLOGIQUE TEMPORELLE D'UN RÉSEAU LEO
%  Sans animation
%
%  Sorties :
%  - beta0(t) : nombre de composantes connexes
%  - beta1_graph(t) : nombre de cycles du graphe non rempli
%  - beta0 et beta1 sur la suite zigzag par unions
%
%  Zigzag construit :
%  G1 -> G1 union G2 <- G2 -> G2 union G3 <- G3 ...
%% ============================================================

%% Paramètres physiques
R_earth = 6371;      % km
h = 550;             % km
R = R_earth + h;     % rayon orbital

mu = 398600;              % km^3/s^2
omega = sqrt(mu / R^3);   % vitesse angulaire orbitale rad/s

%% Paramètres de densité et de constellation Walker Star
lambda = 4e-7;       % satellites / km^2
surface_sphere = 4*pi*R^2;

% Walker Star aléatoire par plans :
% - les P plans orbitaux polaires sont fixés à l'avance ;
% - les satellites sont ensuite placés aléatoirement dans chaque plan.
% Ainsi, P garde un sens géométrique, mais les positions intra-plan ne sont
% pas régulièrement espacées.
P_walker = 28;       % nombre de plans orbitaux
rng('shuffle');      % tirages différents à chaque exécution

% On choisit S à partir de la densité souhaitée, puis on impose N = P*S.
N_target = poissrnd(lambda * surface_sphere);
S_walker = max(1, round(N_target/P_walker));
N = P_walker * S_walker;

fprintf('Nombre cible Poisson : N_cible = %d\n', N_target);
fprintf('Walker Star aléatoire par plans : P = %d plans, S = %d satellites/plan, N = %d satellites\n', ...
    P_walker, S_walker, N);

%% Génération Walker Star aléatoire par plans
% Les plans sont fixés comme dans une Walker Star : les longitudes des
% noeuds ascendants sont réparties sur 180 degrés.
% En revanche, les positions initiales dans chaque plan sont aléatoires :
% u_{p,s}(0) ~ U(0, 2*pi).

Omega_planes = (0:P_walker-1)' * pi/P_walker;   % Walker Star : 180 degrés

Omega = zeros(N,1);
u0 = zeros(N,1);
idx_sat = 1;

for p_idx = 0:P_walker-1
    Omega_p = Omega_planes(p_idx+1);

    for s_idx = 0:S_walker-1
        Omega(idx_sat) = Omega_p;

        % Position initiale aléatoire dans le plan orbital
        u0(idx_sat) = 2*pi*rand;

        idx_sat = idx_sat + 1;
    end
end

% Coordonnées cartésiennes sur une orbite circulaire polaire :
% r = R [cos(Omega) cos(u), sin(Omega) cos(u), sin(u)].
x = R * cos(Omega) .* cos(u0);
y = R * sin(Omega) .* cos(u0);
z = R * sin(u0);

positions0 = [x y z];

% Variables utiles pour la formule angulaire.
% Ici phi désigne le plan orbital, et non une longitude uniforme aléatoire.
phi = Omega;

%% Paramètres des liens et du temps
dmax = 1500;     % km
dt = 60;         % pas temporel en secondes
Tmax = 12000;     % durée totale de simulation

time_values = 0:dt:Tmax;
Nt = length(time_values);

%% Stockage
Positions = cell(Nt,1);
Adjacency = cell(Nt,1);

num_edges = zeros(Nt,1);

% Nombre de liens théorique issu directement de la géométrie Walker Star
% L_W(t) = sum_{i<j} 1_{d_ij(t) <= dmax}
% calculé via l'angle central, sans utiliser pdist.
num_edges_walker_theory = zeros(Nt,1);

beta0 = zeros(Nt,1);
beta1_graph = zeros(Nt,1);
largest_component = zeros(Nt,1);

% Seuil angulaire équivalent à la distance euclidienne dmax
% On tient aussi compte de la limite de visibilité géométrique LOS.
d_LOS = 2*sqrt(R^2 - R_earth^2);
d_eff = min(dmax, d_LOS);

alpha_max = 2*asin(d_eff/(2*R));
cos_alpha_max = cos(alpha_max);

% Minimum théorique approximatif : configuration uniforme sur la sphère
% L_min_th ≈ C(N,2) p_link, avec p_link = d_eff^2/(4R^2).
% Dans le Walker Star, cette valeur sert d'approximation du creux
% lorsque les satellites sont les moins concentrés.
p_link_uniform = d_eff^2/(4*R^2);
L_min_theory_uniform = N*(N-1)/2 * p_link_uniform;

% Maximum théorique approximatif : concentration polaire avec positions
% aléatoires dans les plans.
%
% Dans ce modèle, les plans sont fixés, mais les satellites ne sont pas
% régulièrement espacés dans chaque plan. L'ancienne formule avec
% Delta_u = 2*pi/S n'est donc plus une borne déterministe pertinente.
%
% Pour un pôle donné, un satellite tiré uniformément dans son plan se trouve
% dans une fenêtre angulaire de demi-largeur beta autour du pôle avec
% probabilité p_beta = beta/pi.
%
% Le nombre moyen de satellites dans une calotte polaire d'ouverture beta
% autour d'un pôle est donc M_beta = N*beta/pi.
%
% On utilise ensuite :
% - beta = alpha_max/2 pour une clique garantie en moyenne ;
% - beta = alpha_max pour une estimation haute plus large.

beta_low = alpha_max/2;
beta_up  = alpha_max;

M_pole_low_mean = N * beta_low/pi;
M_pole_up_mean  = N * beta_up/pi;

L_complete = N*(N-1)/2;

% Approximation avec des effectifs moyens réels, pas forcément entiers.
L_max_theory_low = min(L_complete, 2 * M_pole_low_mean*(M_pole_low_mean-1)/2);
L_max_theory_up  = min(L_complete, 2 * M_pole_up_mean *(M_pole_up_mean -1)/2);

% Par sécurité, si l'espérance dans la calotte est < 1, on évite une valeur négative.
L_max_theory_low = max(0, L_max_theory_low);
L_max_theory_up  = max(0, L_max_theory_up);

% Période théorique du nombre de liens.
% La période orbitale complète vaut T_orb = 2*pi/omega.
% Pour le nombre total de liens dans un Walker Star polaire, la configuration
% de connectivité se répète approximativement toutes les demi-orbites.
T_orb = 2*pi/omega;
T_links_theory = T_orb/2;     % = pi/omega

%% ============================================================
%  1. CONSTRUCTION DES GRAPHES TEMPORELS G(t)
%% ============================================================

for k = 1:Nt

    t = time_values(k);

    %% Mouvement orbital Walker Star
    % Chaque satellite reste dans son plan orbital Omega = phi.
    % Sa position angulaire dans le plan évolue à vitesse omega.
    u_t = u0 + omega*t;

    x_t = R * cos(phi) .* cos(u_t);
    y_t = R * sin(phi) .* cos(u_t);
    z_t = R * sin(u_t);

    positions_t = [x_t y_t z_t];

    %% Graphe de lien
    D = squareform(pdist(positions_t));
    A = (D <= dmax) & (D > 0);
    A = sparse(A);

    %% Nombre de liens issu de la formule angulaire Walker Star
    % Angle central gamma_ij(t) entre deux satellites :
    % cos(gamma_ij) = cos(u_i(t)) cos(u_j(t)) cos(Omega_i - Omega_j)
    %               + sin(u_i(t)) sin(u_j(t))
    % Le lien existe si gamma_ij <= alpha_max,
    % soit cos(gamma_ij) >= cos(alpha_max).
    cos_gamma = (cos(u_t)*cos(u_t)') .* cos(phi - phi') + ...
                sin(u_t)*sin(u_t)';
    A_theory = (cos_gamma >= cos_alpha_max) & ~eye(N);
    num_edges_walker_theory(k) = nnz(triu(A_theory,1));

    %% Stockage
    Positions{k} = positions_t;
    Adjacency{k} = A;

    %% Mesures topologiques sur le graphe
    G = graph(A);

    comp = conncomp(G);
    beta0(k) = max(comp);

    comp_sizes = accumarray(comp', 1);
    largest_component(k) = max(comp_sizes);

    E = nnz(triu(A,1));
    num_edges(k) = E;

    % Nombre cyclomatique du graphe :
    % beta1 = E - V + C
    beta1_graph(k) = E - N + beta0(k);
end

% Moyennes temporelles du nombre de liens
% - moyenne théorique Walker Star : calcul géométrique déterministe
% - moyenne empirique : moyenne des liens réellement mesurés dans la simulation
L_walker_mean_theory = mean(num_edges_walker_theory);
L_empirical_mean = mean(num_edges);

% Sinusoïde théorique du nombre de liens.
% On impose :
%   - la moyenne théorique Walker Star,
%   - le minimum théorique uniforme,
%   - la période théorique T_links_theory.
%
% Comme le minimum est supposé atteint à t = 0 :
% L_sin(0) = L_min_theory_uniform.
A_links_theory = L_walker_mean_theory - L_min_theory_uniform;
L_links_sinus_theory = L_walker_mean_theory ...
    - A_links_theory * cos(2*pi*time_values(:)/T_links_theory);

fprintf('Nombre moyen théorique de liens Walker Star : %.2f\n', L_walker_mean_theory);
fprintf('Nombre moyen empirique de liens par simulation : %.2f\n', L_empirical_mean);
fprintf('Minimum théorique uniforme approximatif : %.2f\n', L_min_theory_uniform);
fprintf('Paramètres Walker : P = %d plans, S = %d satellites/plan effectif\n', P_walker, S_walker);
fprintf('Approximation basse du maximum polaire aléatoire : %.2f\n', L_max_theory_low);
fprintf('Approximation haute du maximum polaire aléatoire : %.2f\n', L_max_theory_up);
fprintf('Période théorique des liens : %.2f s\n', T_links_theory);
fprintf('Amplitude de la sinusoïde théorique : %.2f\n', A_links_theory);

%% ============================================================
%  2. GRAPHES TEMPORELS CLASSIQUES
%% ============================================================

figure;
plot(time_values, beta0, 'LineWidth', 2);
grid on;
xlabel('Temps (s)');
ylabel('\beta_0');
title('\beta_0(t) : nombre de composantes connexes');

figure;
plot(time_values, beta1_graph, 'LineWidth', 2);
grid on;
xlabel('Temps (s)');
ylabel('\beta_1 graphe');
title('\beta_1(t) du graphe non rempli');

figure;
plot(time_values, largest_component / N, 'LineWidth', 2);
grid on;
xlabel('Temps (s)');
ylabel('|C_{max}| / N');
title('Fraction de satellites dans la plus grande composante');

figure;
plot(time_values, num_edges, 'LineWidth', 2); hold on;
plot(time_values, L_links_sinus_theory, 'm--', 'LineWidth', 2);
yline(L_walker_mean_theory, 'r--', ...
    sprintf('Moyenne th. Walker = %.1f', L_walker_mean_theory), ...
    'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
yline(L_empirical_mean, 'k:', ...
    sprintf('Moyenne empirique = %.1f', L_empirical_mean), ...
    'LineWidth', 1.8, 'LabelHorizontalAlignment', 'left');
yline(L_min_theory_uniform, 'g-.', ...
    sprintf('Minimum th. uniforme = %.1f', L_min_theory_uniform), ...
    'LineWidth', 1.8, 'LabelHorizontalAlignment', 'left');
yline(L_max_theory_low, 'c--', ...
    sprintf('Approx. basse max polaire = %.1f', L_max_theory_low), ...
    'LineWidth', 1.6, 'LabelHorizontalAlignment', 'left');
yline(L_max_theory_up, 'c:', ...
    sprintf('Approx. haute max polaire = %.1f', L_max_theory_up), ...
    'LineWidth', 1.8, 'LabelHorizontalAlignment', 'left');
grid on;
xlabel('Temps (s)');
ylabel('Nombre de liens');
title('Nombre de liens inter-satellites');
legend('Simulation', ...
       'Sinusoïde théorique', ...
       'Moyenne théorique Walker Star', ...
       'Moyenne empirique simulation', ...
       'Minimum théorique uniforme', ...
       'Approximation basse maximum polaire', ...
       'Approximation haute maximum polaire', ...
       'Location', 'best');


%% ============================================================
%  3. CONSTRUCTION DU ZIGZAG PAR UNIONS
%
%  G1 -> G1 U G2 <- G2 -> G2 U G3 <- G3 ...
%% ============================================================

Nz = 2*Nt - 1;

ZigzagAdjacency = cell(Nz,1);
ZigzagLabels = zeros(Nz,1);

idx = 1;

for k = 1:Nt

    % Graphe réel G_k
    ZigzagAdjacency{idx} = Adjacency{k};
    ZigzagLabels(idx) = k;
    idx = idx + 1;

    % Graphe union G_k U G_{k+1}
    if k < Nt
        ZigzagAdjacency{idx} = Adjacency{k} | Adjacency{k+1};
        ZigzagLabels(idx) = k + 0.5;
        idx = idx + 1;
    end
end

%% ============================================================
%  4. BETTI SUR LA SUITE ZIGZAG
%% ============================================================

beta0_zigzag = zeros(Nz,1);
beta1_zigzag_graph = zeros(Nz,1);
num_edges_zigzag = zeros(Nz,1);
largest_component_zigzag = zeros(Nz,1);

for k = 1:Nz

    A = ZigzagAdjacency{k};
    G = graph(A);

    comp = conncomp(G);
    beta0_zigzag(k) = max(comp);

    comp_sizes = accumarray(comp', 1);
    largest_component_zigzag(k) = max(comp_sizes);

    E = nnz(triu(A,1));
    num_edges_zigzag(k) = E;

    beta1_zigzag_graph(k) = E - N + beta0_zigzag(k);
end

%% ============================================================
%  5. GRAPHES SUR LA SUITE ZIGZAG
%% ============================================================

figure;
plot(ZigzagLabels, beta0_zigzag, '-o', 'LineWidth', 1.5);
grid on;
xlabel('Indice temporel / demi-indice');
ylabel('\beta_0');
title('\beta_0 sur le zigzag par unions');

figure;
plot(ZigzagLabels, beta1_zigzag_graph, '-o', 'LineWidth', 1.5);
grid on;
xlabel('Indice temporel / demi-indice');
ylabel('\beta_1 graphe');
title('\beta_1 du graphe sur le zigzag par unions');

figure;
plot(ZigzagLabels, largest_component_zigzag / N, '-o', 'LineWidth', 1.5);
grid on;
xlabel('Indice temporel / demi-indice');
ylabel('|C_{max}| / N');
title('Composante géante sur le zigzag par unions');

% Moyenne empirique du nombre de liens sur le zigzag par unions
L_empirical_mean_zigzag = mean(num_edges_zigzag);
fprintf('Nombre moyen empirique de liens sur le zigzag : %.2f\n', L_empirical_mean_zigzag);

figure;
plot(ZigzagLabels, num_edges_zigzag, '-o', 'LineWidth', 1.5); hold on;
yline(L_empirical_mean_zigzag, 'k:', ...
    sprintf('Moyenne empirique zigzag = %.1f', L_empirical_mean_zigzag), ...
    'LineWidth', 1.8, 'LabelHorizontalAlignment', 'left');
grid on;
xlabel('Indice temporel / demi-indice');
ylabel('Nombre de liens');
title('Nombre de liens sur le zigzag par unions');
legend('Zigzag par unions', 'Moyenne empirique zigzag', 'Location', 'best');

%% ============================================================
%  6. SAUVEGARDE DES DONNÉES
%% ============================================================

save('leo_zigzag_analysis_results.mat', ...
    'N', 'R', 'R_earth', 'h', 'mu', 'omega', 'lambda', 'dmax', 'dt', 'Tmax', ...
    'time_values', ...
    'Positions', 'Adjacency', ...
    'beta0', 'beta1_graph', 'largest_component', 'num_edges', ...
    'num_edges_walker_theory', 'L_walker_mean_theory', ...
    'L_empirical_mean', 'L_min_theory_uniform', ...
    'L_links_sinus_theory', 'A_links_theory', ...
    'T_orb', 'T_links_theory', ...
    'p_link_uniform', 'd_LOS', 'd_eff', 'alpha_max', ...
    'P_walker', 'S_walker', 'u0', 'phi', ...
    'L_max_theory_low', 'L_max_theory_up', ...
    'beta_low', 'beta_up', 'M_pole_low_mean', 'M_pole_up_mean', ...
    'ZigzagAdjacency', 'ZigzagLabels', ...
    'beta0_zigzag', 'beta1_zigzag_graph', ...
    'largest_component_zigzag', 'num_edges_zigzag', ...
    'L_empirical_mean_zigzag');

fprintf('\nAnalyse terminée.\n');
fprintf('Résultats sauvegardés dans leo_zigzag_analysis_results.mat\n');
