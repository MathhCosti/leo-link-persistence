clear; clc; close all;

%% ============================================================
%  ETUDE TOPOLOGIQUE TEMPORELLE D'UN RESEAU LEO
%  Version Walker-Delta
%
%  Difference avec la version Walker-Star :
%  - les plans orbitaux ne passent pas tous par les poles ;
%  - ils ont une inclinaison commune inc ;
%  - leurs noeuds ascendants sont repartis sur 360 degres ;
%  - les satellites avancent dans leur plan orbital via l'argument de latitude u.
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
lambda = 2e-7;       % satellites / km^2
surface_sphere = 4*pi*R^2;

N = poissrnd(lambda * surface_sphere);
fprintf('Nombre de satellites generes : N = %d\n', N);

%% Parametres Walker-Delta
inc_deg = 53;                  % inclinaison commune des plans orbitaux, en degres
inc = deg2rad(inc_deg);        % radians
P = max(1, round(sqrt(N)));     % nombre de plans orbitaux, choix simple modifiable

% En Walker-Delta, les RAAN sont repartis sur 360 degres.
Omega_planes = 2*pi*(0:P-1)'/P;

% Attribution des satellites aux plans.
plane_id = mod((0:N-1)', P) + 1;
Omega = Omega_planes(plane_id);

% Position initiale dans chaque plan orbital.
% Pour un Walker deterministe strict, remplacer cette ligne par un phasage regulier.
u0 = 2*pi*rand(N,1);           % argument de latitude initial

% Option de phasage regulier de type Walker-Delta, decommenter si besoin :
% f = 1;
% sat_rank = floor((0:N-1)'/P);
% u0 = 2*pi*sat_rank/max(1,ceil(N/P)) + 2*pi*f*(plane_id-1)/N;
% u0 = mod(u0, 2*pi);

%% Positions initiales Walker-Delta
positions0 = walker_delta_positions(R, inc, Omega, u0);

%% Parametres des liens et du temps
dmax = 1500;     % km
dt = 60;         % pas temporel en secondes
Tmax = 12000;    % duree totale de simulation

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

    %% Mouvement orbital Walker-Delta
    u_t = u0 + omega*t;
    positions_t = walker_delta_positions(R, inc, Omega, u_t);

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

    % Nombre cyclomatique du graphe : beta1 = E - V + C
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
title('\beta_0(t) : nombre de composantes connexes - Walker-Delta');

figure;
plot(time_values, beta1_graph, 'LineWidth', 2);
grid on;
xlabel('Temps (s)');
ylabel('\beta_1 graphe');
title('\beta_1(t) du graphe non rempli - Walker-Delta');

figure;
plot(time_values, largest_component / N, 'LineWidth', 2);
grid on;
xlabel('Temps (s)');
ylabel('|C_{max}| / N');
title('Fraction de satellites dans la plus grande composante - Walker-Delta');

figure;
plot(time_values, num_edges, 'LineWidth', 2);
grid on;
xlabel('Temps (s)');
ylabel('Nombre de liens');
title('Nombre de liens inter-satellites - Walker-Delta');

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
title('\beta_0 sur le zigzag par unions - Walker-Delta');

figure;
plot(ZigzagLabels, beta1_zigzag_graph, '-o', 'LineWidth', 1.5);
grid on;
xlabel('Indice temporel / demi-indice');
ylabel('\beta_1 graphe');
title('\beta_1 du graphe sur le zigzag par unions - Walker-Delta');

figure;
plot(ZigzagLabels, largest_component_zigzag / N, '-o', 'LineWidth', 1.5);
grid on;
xlabel('Indice temporel / demi-indice');
ylabel('|C_{max}| / N');
title('Composante geante sur le zigzag par unions - Walker-Delta');

figure;
plot(ZigzagLabels, num_edges_zigzag, '-o', 'LineWidth', 1.5);
grid on;
xlabel('Indice temporel / demi-indice');
ylabel('Nombre de liens');
title('Nombre de liens sur le zigzag par unions - Walker-Delta');

%% ============================================================
%  6. SAUVEGARDE DES DONNEES
%% ============================================================

save('leo_zigzag_analysis_results_delta.mat', ...
    'N', 'R', 'h', 'lambda', 'dmax', 'dt', 'Tmax', 'mu', 'omega', ...
    'inc_deg', 'inc', 'P', 'Omega', 'Omega_planes', 'plane_id', 'u0', ...
    'time_values', ...
    'Positions', 'Adjacency', ...
    'beta0', 'beta1_graph', 'largest_component', 'num_edges', ...
    'ZigzagAdjacency', 'ZigzagLabels', ...
    'beta0_zigzag', 'beta1_zigzag_graph', ...
    'largest_component_zigzag', 'num_edges_zigzag');

fprintf('\nAnalyse terminee.\n');
fprintf('Modele orbital : Walker-Delta, i = %.1f deg, P = %d plans.\n', inc_deg, P);
fprintf('Resultats sauvegardes dans leo_zigzag_analysis_results_delta.mat\n');

%% ============================================================
%  FONCTION LOCALE WALKER-DELTA
%% ============================================================

function positions = walker_delta_positions(R, inc, Omega, u)
    % Positions cartesiennes pour des orbites circulaires Walker-Delta.
    %
    % R     : rayon orbital
    % inc   : inclinaison commune
    % Omega : RAAN de chaque satellite
    % u     : argument de latitude de chaque satellite
    %
    % Formule orbitale circulaire :
    % r = R3(Omega) * R1(inc) * [R cos(u); R sin(u); 0]

    x = R * (cos(Omega).*cos(u) - sin(Omega).*sin(u).*cos(inc));
    y = R * (sin(Omega).*cos(u) + cos(Omega).*sin(u).*cos(inc));
    z = R * (sin(u).*sin(inc));

    positions = [x y z];
end
