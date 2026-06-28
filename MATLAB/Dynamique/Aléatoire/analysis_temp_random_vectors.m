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
dt = 5;         % pas temporel en secondes
Tmax = 2000;    % durée totale de simulation

time_values = 0:dt:Tmax;
Nt = length(time_values);

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
plot(time_values, num_edges, 'LineWidth', 2);
grid on;
xlabel('Temps (s)');
ylabel('Nombre de liens');
title('Nombre de liens inter-satellites');

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

figure;
plot(ZigzagLabels, num_edges_zigzag, '-o', 'LineWidth', 1.5);
grid on;
xlabel('Indice temporel / demi-indice');
ylabel('Nombre de liens');
title('Nombre de liens sur le zigzag par unions');

%% ============================================================
%  6. SAUVEGARDE DES DONNÉES
%% ============================================================

save('leo_zigzag_analysis_random_vectors_results.mat', ...
    'N', 'R', 'h', 'lambda', 'dmax', 'dt', 'Tmax', ...
    'time_values', ...
    'positions0', 'r0', 'v', ...
    'Positions', 'Adjacency', ...
    'beta0', 'beta1_graph', 'largest_component', 'num_edges', ...
    'ZigzagAdjacency', 'ZigzagLabels', ...
    'beta0_zigzag', 'beta1_zigzag_graph', ...
    'largest_component_zigzag', 'num_edges_zigzag');

fprintf('\nAnalyse terminée.\n');
fprintf('Résultats sauvegardés dans leo_zigzag_analysis_random_vectors_results.mat\n');
