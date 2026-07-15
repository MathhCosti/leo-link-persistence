clear; clc; close all;

%% ============================================================
%  ETUDE TOPOLOGIQUE TEMPORELLE D'UN RESEAU LEO
%  Version orbites aleatoires a inclinaison fixe
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
lambda = 4e-7;       % satellites / km^2
surface_sphere = 4*pi*R^2;

N = poissrnd(lambda * surface_sphere);
fprintf('Nombre de satellites generes : N = %d\n', N);

%% Parametres orbitaux aleatoires a inclinaison deterministe
inc_deg = 53;                  % inclinaison commune imposee, en degres
inc = deg2rad(inc_deg);        % radians

% Chaque satellite recoit une orientation de plan aleatoire :
% le RAAN Omega est tire uniformement sur [0, 2*pi[.
% Les plans ne sont donc plus espaces regulierement.
Omega = 2*pi*rand(N,1);

% La phase initiale dans le plan orbital est egalement aleatoire.
u0 = 2*pi*rand(N,1);

% Pour compatibilite avec les sauvegardes et affichages precedents,
% on considere ici un plan individuel par satellite.
P = N;
plane_id = (1:N)';
Omega_planes = Omega;

%% Positions initiales orbites aleatoires a inclinaison fixe
positions0 = walker_delta_positions(R, inc, Omega, u0);

%% Parametres des liens et du temps
dmax = 1500;     % km
dt = 20;         % pas temporel en secondes
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

% Distribution empirique des tailles de composantes :
% component_size_counts(k,s) = nombre de composantes de taille s au temps k.
component_size_counts = zeros(Nt, N);

%% ============================================================
%  1. CONSTRUCTION DES GRAPHES TEMPORELS G(t)
%% ============================================================

for k = 1:Nt

    t = time_values(k);

    %% Mouvement orbital orbites aleatoires a inclinaison fixe
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

    % Comptage du nombre de composantes pour chaque taille s.
    size_hist = accumarray(comp_sizes, 1, [N, 1]);
    component_size_counts(k,:) = size_hist.';

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
title('\beta_0(t) : nombre de composantes connexes - orbites aleatoires a inclinaison fixe');

figure;
plot(time_values, beta1_graph, 'LineWidth', 2);
grid on;
xlabel('Temps (s)');
ylabel('\beta_1 graphe');
title('\beta_1(t) du graphe non rempli - orbites aleatoires a inclinaison fixe');

figure;
plot(time_values, largest_component / N, 'LineWidth', 2);
grid on;
xlabel('Temps (s)');
ylabel('|C_{max}| / N');
title('Fraction de satellites dans la plus grande composante - orbites aleatoires a inclinaison fixe');

figure;
plot(time_values, num_edges, 'LineWidth', 2);
grid on;
xlabel('Temps (s)');
ylabel('Nombre de liens');
title('Nombre de liens inter-satellites - orbites aleatoires a inclinaison fixe');


%% ============================================================
%  2.b DISTRIBUTION EMPIRIQUE MOYENNE DES TAILLES DE COMPOSANTES
%% ============================================================

% Nombre moyen temporel de composantes de taille s.
mean_component_count_by_size = mean(component_size_counts, 1);

% Taille maximale effectivement observee.
last_nonzero_size = find(mean_component_count_by_size > 0, 1, 'last');
if isempty(last_nonzero_size)
    last_nonzero_size = 1;
end

component_sizes_axis = 1:last_nonzero_size;
mean_component_count_plot = mean_component_count_by_size(component_sizes_axis);

% Fraction moyenne des composantes appartenant a chaque classe de taille.
total_mean_components = sum(mean_component_count_plot);
if total_mean_components > 0
    component_size_fraction = mean_component_count_plot / total_mean_components;
else
    component_size_fraction = zeros(size(mean_component_count_plot));
end

figure;
bar(component_sizes_axis, mean_component_count_plot);
grid on;
xlabel('Taille s de la composante');
ylabel('Nombre moyen temporel de composantes de taille s');
title('Distribution empirique moyenne des tailles de composantes');

figure;
bar(component_sizes_axis, component_size_fraction);
grid on;
xlabel('Taille s de la composante');
ylabel('Fraction moyenne des composantes');
title('Distribution normalisee des tailles de composantes');


%% ============================================================
%  2.c COMPARAISON EMPIRIQUE / THEORIQUE POUR N1, N2 ET N3
%% ============================================================

if dmax >= 2*R
    alpha_max = pi;
    p_link = 1;
else
    alpha_max = 2*asin(dmax/(2*R));
    p_link = 2*mean(num_edges)/(N*(N-1));
end

c2 = 1 + 3*sqrt(3)/(4*pi);
c3_conn = 1 + 3*sqrt(3)/(2*pi);
c3_union = 1.8;

N1_theory = N*(1-p_link)^(N-1);

if N >= 2
    N2_theory = nchoosek(N,2)*p_link ...
        * max(0,1-c2*p_link)^(N-2);
else
    N2_theory = 0;
end

if N >= 3
    N3_theory = nchoosek(N,3)*c3_conn*p_link^2 ...
        * max(0,1-c3_union*p_link)^(N-3);
else
    N3_theory = 0;
end

N1_emp_time = component_size_counts(:,1);
N2_emp_time = zeros(Nt,1);
N3_emp_time = zeros(Nt,1);

if N >= 2
    N2_emp_time = component_size_counts(:,2);
end
if N >= 3
    N3_emp_time = component_size_counts(:,3);
end

N1_emp_mean = mean(N1_emp_time);
N2_emp_mean = mean(N2_emp_time);
N3_emp_mean = mean(N3_emp_time);

figure;
hold on;
grid on;
plot(time_values, N1_emp_time, 'LineWidth', 1.3, ...
    'DisplayName', 'N_1 empirique');
plot(time_values, N2_emp_time, 'LineWidth', 1.3, ...
    'DisplayName', 'N_2 empirique');
plot(time_values, N3_emp_time, 'LineWidth', 1.3, ...
    'DisplayName', 'N_3 empirique');

yline(N1_theory, '--', 'LineWidth', 1.8, ...
    'DisplayName', sprintf('N_1 theorique = %.3f', N1_theory));
yline(N2_theory, '--', 'LineWidth', 1.8, ...
    'DisplayName', sprintf('N_2 theorique = %.3f', N2_theory));
yline(N3_theory, '--', 'LineWidth', 1.8, ...
    'DisplayName', sprintf('N_3 theorique = %.3f', N3_theory));

xlabel('Temps (s)');
ylabel('Nombre de composantes');
title('Comparaison temporelle de N_1, N_2 et N_3');
legend('Location', 'best');
hold off;

empirical_N123 = [N1_emp_mean, N2_emp_mean, N3_emp_mean];
theoretical_N123 = [N1_theory, N2_theory, N3_theory];

figure;
bar([empirical_N123; theoretical_N123].');
grid on;
xticks(1:3);
xticklabels({'N_1 : isoles', 'N_2 : dimeres', 'N_3 : trimeres'});
ylabel('Nombre moyen de composantes');
title('N_1, N_2, N_3 : empirique vs theorie geometrique');
legend('Empirique', 'Theorie', 'Location', 'best');

relative_error_N123 = abs(theoretical_N123-empirical_N123) ...
    ./ max(empirical_N123, eps);

figure;
bar(100*relative_error_N123);
grid on;
xticks(1:3);
xticklabels({'N_1', 'N_2', 'N_3'});
ylabel('Erreur relative (%)');
title('Erreur relative des approximations de N_1, N_2 et N_3');

fprintf('\n--- Comparaison empirique / theorique N1, N2, N3 ---\n');
fprintf('N1 empirique = %.4f | theorie = %.4f | erreur = %.2f %%\n', ...
    N1_emp_mean, N1_theory, 100*relative_error_N123(1));
fprintf('N2 empirique = %.4f | theorie = %.4f | erreur = %.2f %%\n', ...
    N2_emp_mean, N2_theory, 100*relative_error_N123(2));
fprintf('N3 empirique = %.4f | theorie = %.4f | erreur = %.2f %%\n', ...
    N3_emp_mean, N3_theory, 100*relative_error_N123(3));

% Taille moyenne d'une composante, moyennee dans le temps.
mean_component_size_time = N ./ beta0;
mean_component_size = mean(mean_component_size_time);

fprintf('\n--- Distribution empirique des tailles de composantes ---\n');
fprintf('Taille moyenne temporelle d''une composante : %.4f satellites\n', ...
    mean_component_size);
fprintf('Taille maximale observee                    : %d satellites\n', ...
    last_nonzero_size);

% Affichage des classes de taille les plus frequentes.
[sorted_counts, sorted_sizes] = sort(mean_component_count_by_size, 'descend');
n_display = min(10, nnz(sorted_counts > 0));

fprintf('Tailles les plus representees en moyenne :\n');
for ii = 1:n_display
    fprintf('  taille %d : %.4f composantes en moyenne\n', ...
        sorted_sizes(ii), sorted_counts(ii));
end

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
title('\beta_0 sur le zigzag par unions - orbites aleatoires a inclinaison fixe');

figure;
plot(ZigzagLabels, beta1_zigzag_graph, '-o', 'LineWidth', 1.5);
grid on;
xlabel('Indice temporel / demi-indice');
ylabel('\beta_1 graphe');
title('\beta_1 du graphe sur le zigzag par unions - orbites aleatoires a inclinaison fixe');

figure;
plot(ZigzagLabels, largest_component_zigzag / N, '-o', 'LineWidth', 1.5);
grid on;
xlabel('Indice temporel / demi-indice');
ylabel('|C_{max}| / N');
title('Composante geante sur le zigzag par unions - orbites aleatoires a inclinaison fixe');

figure;
plot(ZigzagLabels, num_edges_zigzag, '-o', 'LineWidth', 1.5);
grid on;
xlabel('Indice temporel / demi-indice');
ylabel('Nombre de liens');
title('Nombre de liens sur le zigzag par unions - orbites aleatoires a inclinaison fixe');

%% ============================================================
%  6. SAUVEGARDE DES DONNEES
%% ============================================================

save('leo_zigzag_analysis_results_delta.mat', ...
    'N', 'R', 'h', 'lambda', 'dmax', 'dt', 'Tmax', 'mu', 'omega', ...
    'inc_deg', 'inc', 'P', 'Omega', 'Omega_planes', 'plane_id', 'u0', ...
    'time_values', ...
    'Positions', 'Adjacency', ...
    'beta0', 'beta1_graph', 'largest_component', 'num_edges', ...
    'component_size_counts', 'mean_component_count_by_size', ...
    'component_size_fraction', 'mean_component_size_time', ...
    'mean_component_size', ...
    'N1_emp_time', 'N2_emp_time', 'N3_emp_time', ...
    'N1_emp_mean', 'N2_emp_mean', 'N3_emp_mean', ...
    'N1_theory', 'N2_theory', 'N3_theory', ...
    'p_link', 'alpha_max', 'c2', 'c3_conn', 'c3_union', ...
    'relative_error_N123', ...
    'ZigzagAdjacency', 'ZigzagLabels', ...
    'beta0_zigzag', 'beta1_zigzag_graph', ...
    'largest_component_zigzag', 'num_edges_zigzag');

fprintf('\nAnalyse terminee.\n');
fprintf('Modele orbital : RAAN et phases initiaux aleatoires, inclinaison fixe i = %.1f deg, P = %d plans.\n', inc_deg, P);
fprintf('Resultats sauvegardes dans leo_zigzag_analysis_results_delta.mat\n');

%% ============================================================
%  FONCTION LOCALE WALKER-DELTA
%% ============================================================

function positions = walker_delta_positions(R, inc, Omega, u)
    % Positions cartesiennes pour des orbites circulaires orbites aleatoires a inclinaison fixe.
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
