clear; clc; close all;

%% ============================================================
%  ÉTUDE TOPOLOGIQUE TEMPORELLE D'UN RÉSEAU LEO
%  Mouvement avec vecteurs aléatoires tangentiels
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
omega = sqrt(mu / R^3);   % vitesse angulaire rad/s

%% Paramètres du processus de Poisson
lambda = 5e-7;       % satellites / km^2
surface_sphere = 4*pi*R^2;

N = poissrnd(lambda * surface_sphere);

fprintf('Nombre de satellites générés : N = %d\n', N);

%% Génération uniforme des positions initiales sur la sphère
u = rand(N,1);
phi = 2*pi*rand(N,1);
theta = acos(1 - 2*u);

x = R * sin(theta) .* cos(phi);
y = R * sin(theta) .* sin(phi);
z = R * cos(theta);

positions0 = [x y z];

%% ============================================================
%  Directions aléatoires tangentes à la sphère
%% ============================================================

% Vecteurs radiaux unitaires
r0 = positions0 / R;              % N x 3

% Vecteurs aléatoires 3D quelconques
a = randn(N,3);

% Projection dans le plan tangent à la sphère :
% v = a - (a.r0) r0
v = a - sum(a .* r0, 2) .* r0;

% Normalisation : chaque v est unitaire et tangent à la sphère
v = v ./ vecnorm(v, 2, 2);

% Optionnel : sens aléatoire +1 / -1
sens = sign(rand(N,1) - 0.5);
sens(sens == 0) = 1;
v = sens .* v;

%% Paramètres des liens et du temps
dmax = 1500;     % km
dt = 60;         % pas temporel en secondes
Tmax = 12000;    % durée totale de simulation

time_values = 0:dt:Tmax;
Nt = length(time_values);

%% ============================================================
%  Valeurs théoriques approximatives
%% ============================================================

% Angle central maximal correspondant à une distance corde dmax
alpha_max = 2 * asin(min(1, dmax/(2*R)));

% Probabilité théorique qu'une paire de satellites soit liée
% pour deux points uniformes indépendants sur la sphère
p_link_theory = (1 - cos(alpha_max)) / 2;

% Forme équivalente pour une distance corde : p = dmax^2/(4R^2)
% tant que dmax <= 2R
p_link_theory_corde = dmax^2 / (4*R^2);

num_pairs = nchoosek(N,2);
num_edges_theory = num_pairs * p_link_theory;

% Degré moyen théorique
mean_degree_theory = (N-1) * p_link_theory;

% Approximation simple de beta0 par le nombre de sommets isolés.
% Attention : c'est une approximation, surtout valable quand les petites
% composantes non isolées sont rares.
beta0_theory_iso = N * (1 - p_link_theory)^(N-1);

% Approximation de beta1 par la formule cyclomatique E - V + C,
% en remplaçant E et C par leurs approximations théoriques.
beta1_theory_graph = max(0, num_edges_theory - N + beta0_theory_iso);

% Pour les graphes par unions, approximation naïve si G(t) et G(t+dt)
% étaient indépendants. En réalité ils sont corrélés, donc cette valeur
% sert surtout de repère.
p_link_union_theory_ind = 1 - (1 - p_link_theory)^2;
num_edges_union_theory_ind = num_pairs * p_link_union_theory_ind;
beta0_union_theory_iso_ind = N * (1 - p_link_union_theory_ind)^(N-1);
beta1_union_theory_graph_ind = max(0, num_edges_union_theory_ind - N + beta0_union_theory_iso_ind);

fprintf('p_link théorique = %.6f, soit %.3f %%\n', p_link_theory, 100*p_link_theory);
fprintf('Nombre moyen théorique de liens = %.2f\n', num_edges_theory);
fprintf('Degré moyen théorique = %.2f\n', mean_degree_theory);
fprintf('Approx beta0 théorique isolés = %.2f\n', beta0_theory_iso);
fprintf('Approx beta1 théorique graphe = %.2f\n', beta1_theory_graph);

%% Stockage
Positions = cell(Nt,1);
Adjacency = cell(Nt,1);

num_edges = zeros(Nt,1);
beta0 = zeros(Nt,1);
beta1_graph = zeros(Nt,1);
largest_component = zeros(Nt,1);

%% ============================================================
%  1. CONSTRUCTION DES GRAPHES TEMPORELS G(t)
%% ============================================================

for k = 1:Nt

    t = time_values(k);

    %% Mouvement sur grand cercle avec direction tangentielle aléatoire
    % Formule géodésique sur la sphère :
    % r(t) = R [ r0 cos(omega t) + v sin(omega t) ]
    positions_t = R * (r0 * cos(omega*t) + v * sin(omega*t));

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
    num_edges(k) = E;

    % Nombre cyclomatique du graphe :
    % beta1 = E - V + C
    beta1_graph(k) = E - N + beta0(k);
end

%% ============================================================
%  2. GRAPHES TEMPORELS CLASSIQUES
%% ============================================================

figure;
plot(time_values, beta0, 'LineWidth', 2);
hold on;
yline(beta0_theory_iso, '--', 'Approx. théorique isolés', 'LineWidth', 1.5);
grid on;
xlabel('Temps (s)');
ylabel('\beta_0');
title('\beta_0(t) : nombre de composantes connexes');
legend('Simulation', 'Approx. théorique', 'Location', 'best');

figure;
plot(time_values, beta1_graph, 'LineWidth', 2);
hold on;
yline(beta1_theory_graph, '--', 'Approx. théorique', 'LineWidth', 1.5);
grid on;
xlabel('Temps (s)');
ylabel('\beta_1 graphe');
title('\beta_1(t) du graphe non rempli');
legend('Simulation', 'Approx. théorique', 'Location', 'best');

figure;
plot(time_values, largest_component / N, 'LineWidth', 2);
grid on;
xlabel('Temps (s)');
ylabel('|C_{max}| / N');
title('Fraction de satellites dans la plus grande composante');

figure;
plot(time_values, num_edges, 'LineWidth', 2);
hold on;
yline(num_edges_theory, '--', 'E[liens] théorique', 'LineWidth', 1.5);
grid on;
xlabel('Temps (s)');
ylabel('Nombre de liens');
title('Nombre de liens inter-satellites');
legend('Simulation', 'Valeur théorique', 'Location', 'best');

figure;
p_link_emp = num_edges / num_pairs;
plot(time_values, p_link_emp, 'LineWidth', 2);
hold on;
yline(p_link_theory, '--', 'p_{link} théorique', 'LineWidth', 1.5);
grid on;
xlabel('Temps (s)');
ylabel('Probabilité de lien');
title('Probabilité empirique de lien p_{link}(t)');
legend('Simulation', 'Valeur théorique', 'Location', 'best');

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
hold on;
yline(beta0_theory_iso, '--', 'Instantané théorique', 'LineWidth', 1.5);
yline(beta0_union_theory_iso_ind, ':', 'Union théorique indép.', 'LineWidth', 1.5);
grid on;
xlabel('Indice temporel / demi-indice');
ylabel('\beta_0');
title('\beta_0 sur le zigzag par unions');
legend('Simulation zigzag', 'Approx. instantanée', 'Approx. union indépendante', 'Location', 'best');

figure;
plot(ZigzagLabels, beta1_zigzag_graph, '-o', 'LineWidth', 1.5);
hold on;
yline(beta1_theory_graph, '--', 'Instantané théorique', 'LineWidth', 1.5);
yline(beta1_union_theory_graph_ind, ':', 'Union théorique indép.', 'LineWidth', 1.5);
grid on;
xlabel('Indice temporel / demi-indice');
ylabel('\beta_1 graphe');
title('\beta_1 du graphe sur le zigzag par unions');
legend('Simulation zigzag', 'Approx. instantanée', 'Approx. union indépendante', 'Location', 'best');

figure;
plot(ZigzagLabels, largest_component_zigzag / N, '-o', 'LineWidth', 1.5);
grid on;
xlabel('Indice temporel / demi-indice');
ylabel('|C_{max}| / N');
title('Composante géante sur le zigzag par unions');

figure;
plot(ZigzagLabels, num_edges_zigzag, '-o', 'LineWidth', 1.5);
hold on;
yline(num_edges_theory, '--', 'Instantané théorique', 'LineWidth', 1.5);
yline(num_edges_union_theory_ind, ':', 'Union théorique indép.', 'LineWidth', 1.5);
grid on;
xlabel('Indice temporel / demi-indice');
ylabel('Nombre de liens');
title('Nombre de liens sur le zigzag par unions');
legend('Simulation zigzag', 'E[liens] instantané', 'E[liens] union indépendante', 'Location', 'best');

figure;
p_link_zigzag_emp = num_edges_zigzag / num_pairs;
plot(ZigzagLabels, p_link_zigzag_emp, '-o', 'LineWidth', 1.5);
hold on;
yline(p_link_theory, '--', 'p_{link} instantané', 'LineWidth', 1.5);
yline(p_link_union_theory_ind, ':', 'p_{link} union indép.', 'LineWidth', 1.5);
grid on;
xlabel('Indice temporel / demi-indice');
ylabel('Probabilité de lien');
title('Probabilité de lien sur le zigzag par unions');
legend('Simulation zigzag', 'p_{link} instantané', 'p_{link} union indépendante', 'Location', 'best');

%% ============================================================
%  6. SAUVEGARDE DES DONNÉES
%% ============================================================

save('leo_zigzag_analysis_random_vectors_results.mat', ...
    'N', 'R', 'h', 'lambda', 'dmax', 'dt', 'Tmax', ...
    'alpha_max', 'p_link_theory', 'p_link_theory_corde', ...
    'num_pairs', 'num_edges_theory', 'mean_degree_theory', ...
    'beta0_theory_iso', 'beta1_theory_graph', ...
    'p_link_union_theory_ind', 'num_edges_union_theory_ind', ...
    'beta0_union_theory_iso_ind', 'beta1_union_theory_graph_ind', ...
    'time_values', ...
    'positions0', 'r0', 'v', ...
    'Positions', 'Adjacency', ...
    'beta0', 'beta1_graph', 'largest_component', 'num_edges', ...
    'ZigzagAdjacency', 'ZigzagLabels', ...
    'beta0_zigzag', 'beta1_zigzag_graph', ...
    'largest_component_zigzag', 'num_edges_zigzag');

fprintf('\nAnalyse terminée.\n');
fprintf('Résultats sauvegardés dans leo_zigzag_analysis_random_vectors_results.mat\n');
