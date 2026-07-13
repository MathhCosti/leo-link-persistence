clear; clc; close all;

%% ============================================================
%  ETUDE TOPOLOGIQUE TEMPORELLE D'UN RESEAU LEO
%  Version sans calculs/graphes du nombre de liens
%  et sans calculs/graphes de p_disp
%
%  Sorties :
%  - beta0(t) : nombre de composantes connexes
%  - beta1_graph(t) : nombre de cycles du graphe non rempli
%  - beta0 et beta1 sur la suite zigzag par unions
%
%  Zigzag construit :
%  G1 -> G1 union G2 <- G2 -> G2 union G3 <- G3 ...
%% ============================================================

%% Parametres physiques
R_earth = 6371;      % km
h = 550;             % km
R = R_earth + h;     % rayon orbital

mu = 398600;              % km^3/s^2
omega = sqrt(mu / R^3);   % vitesse angulaire orbitale rad/s

%% Parametres du processus de Poisson
lambda = 4e-7;       % satellites / km^2
surface_sphere = 4*pi*R^2;

N = poissrnd(lambda * surface_sphere);

fprintf('Nombre de satellites generes : N = %d\n', N);

%% Generation uniforme sur l'orbite choisie (Walker Star)
% Omega est tire sur [0, 2*pi) afin de conserver une orientation dirigee
% du plan orbital autour de l'axe des poles. Les positions restent
% uniformes en longueur d'arc sur chaque orbite grace au tirage uniforme
% de u0 sur [0, 2*pi).
Omega = 2*pi * rand(N,1);
u0 = 2*pi * rand(N,1);

x = R * cos(u0) .* cos(Omega);
y = R * cos(u0) .* sin(Omega);
z = R * sin(u0);

positions0 = [x y z]; %#ok<NASGU>

%% Deux parties separees par le plan polaire y = 0
% Le plan y = 0 contient l'axe z et donc les poles Nord et Sud.
% Les satellites sont repartis en deux demi-espaces selon leur position
% initiale :
%   y0 >= 0  -> sens orbital +1 ;
%   y0 <  0  -> sens orbital -1.
% Le signe est fixe une seule fois a t = 0 et reste constant pendant toute
% la simulation, meme lorsque le satellite traverse ensuite le plan y = 0.
y0 = y;
rotation_sign = ones(N,1);
rotation_sign(y0 < 0) = -1;

group_plus = (rotation_sign == 1);
group_minus = (rotation_sign == -1);

fprintf('Partie y0 >= 0 / sens + : %d satellites | Partie y0 < 0 / sens - : %d satellites\n', ...
    nnz(group_plus), nnz(group_minus));

%% Parametres des liens et du temps
dmax = 1500;      % km
dt = 60;          % pas temporel en secondes
Tmax = 12000;     % duree totale de simulation

%% Approximations théoriques de beta0 pour le lambda choisi
% Approximation de la probabilité que deux satellites soient liés.
% ATTENTION : la formule ci-dessous est exacte pour deux points uniformes
% en surface sur la sphère. Avec u0 uniforme sur les orbites Walker Star,
% la densité spatiale n'est plus uniforme en latitude ; cette valeur sert
% donc uniquement de référence homogène.
alpha_max = 2 * asin(min(dmax/(2*R), 1));
p_link = (1 - cos(alpha_max)) / 2;

p = p_link;

% Nombre moyen de liens
E_theory = N * (N - 1) / 2 * p;

% Modèle "connectés" / forêt : chaque lien réduit beta0 d'environ 1,
% jusqu'à l'apparition d'une composante principale.
beta0_theory_connected = max(N - E_theory, 1);

%% ------------------------------------------------------------
%  Modèle analytique corrigé géométriquement jusqu'aux trimères
%
%  On utilise :
%     E[beta0] ~= 1 + N1 + N2 + N3
%
%  N1 : satellites isolés
%  N2 : composantes de taille 2, avec correction d'aire d'union c2
%  N3 : composantes de taille 3, avec correction de connexité c3_conn
%       et correction d'aire d'union c3_union
%
%  Les constantes c2_union et c3_conn viennent d'une approximation
%  géométrique locale plane du graphe de disque aléatoire.
%  c3_union est un coefficient effectif pour l'aire moyenne de l'union
%  de trois disques conditionnellement à la connexité du triplet.
%% ------------------------------------------------------------

% Coefficient moyen d'aire d'union de deux disques connectés :
% A_union,2 ~= c2_union * pi*r^2
c2_union = 1.4135;

% Probabilité de connexité interne d'un triplet géométrique :
% P(G3 connexe) ~= c3_conn * p^2
c3_conn = 1.827;

% Coefficient effectif d'aire d'union de trois disques connectés :
% A_union,3 ~= c3_union * pi*r^2
c3_union = 1.80;

% Sécurités numériques : les bases des puissances doivent rester positives.
q1_ext = max(1 - p, 0);
q2_ext = max(1 - c2_union*p, 0);
q3_ext = max(1 - c3_union*p, 0);

% Composantes de taille 1 : satellites isolés
N1_theory = N * q1_ext^(N - 1);

% Composantes de taille 2 : dimères isolés
if N >= 2
    N2_theory = nchoosek(N, 2) * p * q2_ext^(N - 2);
else
    N2_theory = 0;
end

% Composantes de taille 3 : trimères isolés
if N >= 3
    p_conn_3_geom = c3_conn * p^2;
    p_conn_3_geom = min(max(p_conn_3_geom, 0), 1);  % sécurité numérique
    N3_theory = nchoosek(N, 3) * p_conn_3_geom * q3_ext^(N - 3);
else
    p_conn_3_geom = 0;
    N3_theory = 0;
end

% Pour comparaison : ancienne version Erdos-Renyi indépendante
if N >= 2
    N2_theory_ER = nchoosek(N, 2) * p * (1 - p)^(2*(N - 2));
else
    N2_theory_ER = 0;
end

if N >= 3
    p_conn_3_ER = 3*p^2 - 2*p^3;
    N3_theory_ER = nchoosek(N, 3) * p_conn_3_ER * (1 - p)^(3*(N - 3));
else
    p_conn_3_ER = 0;
    N3_theory_ER = 0;
end

% Modèles successifs de beta0
% Le +2 représente les composantes principales restantes.
beta0_theory_isolated = 2 + N1_theory;
beta0_theory_isolated_dimers = 2 + N1_theory + N2_theory;
beta0_theory_isolated_dimers_trimers = 2 + N1_theory + N2_theory + N3_theory;

fprintf('p_link théorique : %.6f\n', p_link);
fprintf('|E| théorique moyen : %.3f\n', E_theory);
fprintf('c2_union : %.4f\n', c2_union);
fprintf('c3_conn : %.4f\n', c3_conn);
fprintf('c3_union : %.4f\n', c3_union);
fprintf('N1 théorie, isolés : %.3f\n', N1_theory);
fprintf('N2 théorie géométrique, dimères : %.3f\n', N2_theory);
fprintf('N3 théorie géométrique, trimères : %.3f\n', N3_theory);
fprintf('N2 théorie ER, dimères : %.3f\n', N2_theory_ER);
fprintf('N3 théorie ER, trimères : %.3f\n', N3_theory_ER);
fprintf('beta0 théorie connectés : %.3f\n', beta0_theory_connected);
fprintf('beta0 théorie isolés : %.3f\n', beta0_theory_isolated);
fprintf('beta0 théorie isolés + dimères géométriques : %.3f\n', beta0_theory_isolated_dimers);
fprintf('beta0 théorie isolés + dimères + trimères géométriques : %.3f\n', beta0_theory_isolated_dimers_trimers);

time_values = 0:dt:Tmax;
Nt = length(time_values);

%% Stockage
Positions = cell(Nt,1);
Adjacency = cell(Nt,1);

beta0 = zeros(Nt,1);
beta1_graph = zeros(Nt,1);
largest_component = zeros(Nt,1);

%% ============================================================
%  1. CONSTRUCTION DES GRAPHES TEMPORELS G(t)
%% ============================================================

for k = 1:Nt

    t = time_values(k);

    %% Mouvement orbital
    % Le plan orbital Omega reste constant. Seule la phase orbitale evolue.
    u_t = u0 + rotation_sign * omega * t;

    x_t = R * cos(u_t) .* cos(Omega);
    y_t = R * cos(u_t) .* sin(Omega);
    z_t = R * sin(u_t);

    positions_t = [x_t y_t z_t];

    %% Graphe de lien
    D = squareform(pdist(positions_t));
    A = (D <= dmax) & (D > 0);
    A = sparse(A);

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

    % Nombre cyclomatique du graphe : beta1 = E - V + C
    beta1_graph(k) = E - N + beta0(k);
end

%% ============================================================
%  2. GRAPHES TEMPORELS CLASSIQUES
%% ============================================================

figure;
plot(time_values, beta0, 'LineWidth', 2); hold on;
yline(beta0_theory_connected, '--', 'Théorie connectés', 'LineWidth', 1.5);
yline(beta0_theory_isolated, ':', 'Théorie isolés', 'LineWidth', 1.5);
yline(beta0_theory_isolated_dimers, '-.', 'Théorie isolés + dimères geom.', 'LineWidth', 1.5);
yline(beta0_theory_isolated_dimers_trimers, '-', 'Théorie isolés + dimères + trimères geom.', 'LineWidth', 1.0);
grid on;
xlabel('Temps (s)');
ylabel('\beta_0');
title('\beta_0(t) : nombre de composantes connexes');
legend('Simulation', ...
       'Théorie sparse', ...
       'Théorie isolés', ...
       'Théorie isolés + dimères',...
       'Théorie isolés + trimères', ...
       'Location', 'best');

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

    % Graphe reel G_k
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
largest_component_zigzag = zeros(Nz,1);

for k = 1:Nz

    A = ZigzagAdjacency{k};
    G = graph(A);

    comp = conncomp(G);
    beta0_zigzag(k) = max(comp);

    comp_sizes = accumarray(comp', 1);
    largest_component_zigzag(k) = max(comp_sizes);

    E = nnz(triu(A,1));

    beta1_zigzag_graph(k) = E - N + beta0_zigzag(k);
end

%% ============================================================
%  5. GRAPHES SUR LA SUITE ZIGZAG
%% ============================================================

figure;
plot(ZigzagLabels, beta0_zigzag, '-o', 'LineWidth', 1.5); hold on;
yline(beta0_theory_connected, '--', 'Théorie connectés', 'LineWidth', 1.5);
yline(beta0_theory_isolated, ':', 'Théorie isolés', 'LineWidth', 1.5);
yline(beta0_theory_isolated_dimers, '-.', 'Théorie isolés + dimères geom.', 'LineWidth', 1.5);
yline(beta0_theory_isolated_dimers_trimers, '-', 'Théorie isolés + dimères + trimères geom.', 'LineWidth', 1.0);
grid on;
xlabel('Indice temporel / demi-indice');
ylabel('\beta_0');
title('\beta_0 sur le zigzag par unions');
legend('Simulation zigzag', ...
       'Théorie sparse', ...
       'Théorie isolés', ...
       'Théorie isolés + dimères',...
       'Théorie isolés + trimères', ...
       'Location', 'best');

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
title('Composante geante sur le zigzag par unions');

%% ============================================================
%  6. SAUVEGARDE DES DONNEES
%% ============================================================

save('leo_zigzag_analysis_results.mat', ...
    'N', 'R', 'h', 'lambda', 'dmax', 'dt', 'Tmax', 'Omega', 'u0', ...
    'y0', 'rotation_sign', 'group_plus', 'group_minus', ...
    'p_link', 'E_theory', ...
    'N1_theory', 'N2_theory', 'N3_theory', ...
    'N2_theory_ER', 'N3_theory_ER', ...
    'c2_union', 'c3_conn', 'c3_union', ...
    'beta0_theory_connected', 'beta0_theory_isolated', ...
    'beta0_theory_isolated_dimers', 'beta0_theory_isolated_dimers_trimers', ...
    'time_values', ...
    'Positions', 'Adjacency', ...
    'beta0', 'beta1_graph', 'largest_component', ...
    'ZigzagAdjacency', 'ZigzagLabels', ...
    'beta0_zigzag', 'beta1_zigzag_graph', ...
    'largest_component_zigzag');

fprintf('\nAnalyse terminee.\n');
fprintf('Resultats sauvegardes dans leo_zigzag_analysis_results.mat\n');
