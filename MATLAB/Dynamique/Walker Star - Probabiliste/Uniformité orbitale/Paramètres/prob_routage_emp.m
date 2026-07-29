clear; clc; close all;

%% ============================================================
%  PROBABILITE DE ROUTAGE TEMPORELLE
%  Version orbites aleatoires a inclinaison fixe
%% ============================================================

%% Parametres physiques
R_earth = 6371;      % km
h = 550;             % km
R = R_earth + h;     % rayon orbital

mu = 398600;              % km^3/s^2
omega = sqrt(mu / R^3);   % vitesse angulaire orbitale rad/s

%% Parametres reseau
lambda = 2e-7;       % satellites / km^2
surface_sphere = 4*pi*R^2;
N = poissrnd(lambda * surface_sphere);

dmax = 1500;         % km

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

%% Temps de simulation
dt = 10;              % pas de temps en secondes
Tmax = 6000;          % duree totale en secondes
time_values = 0:dt:Tmax;

P_routing_time = zeros(length(time_values), 1);
num_edges_time = zeros(length(time_values), 1);
largest_component_time = zeros(length(time_values), 1);

for k = 1:length(time_values)

    t = time_values(k);

    %% Mouvement orbital orbites aleatoires a inclinaison fixe
    u_t = u0 + omega*t;
    positions_t = walker_delta_positions(R, inc, Omega, u_t);

    %% Construction du graphe a l'instant t
    D = squareform(pdist(positions_t));
    A = (D <= dmax) & (D > 0);
    G = graph(A);

    %% Probabilite de routage multi-sauts
    comp = conncomp(G);
    component_sizes = accumarray(comp', 1);

    if N > 1
        P_routing_time(k) = sum(component_sizes .* (component_sizes - 1)) / (N * (N - 1));
    else
        P_routing_time(k) = 0;
    end

    num_edges_time(k) = nnz(triu(A,1));
    largest_component_time(k) = max(component_sizes) / N;
end

%% Affichage
figure;
plot(time_values, P_routing_time, 'LineWidth', 1.5);
grid on;
xlabel('Temps (s)');
ylabel('P(routage multi-sauts)');
title(sprintf('Evolution temporelle de la probabilite de routage - orbites aleatoires a inclinaison fixe, i = %.1f deg', inc_deg));

figure;
plot(time_values, largest_component_time, 'LineWidth', 1.5);
grid on;
xlabel('Temps (s)');
ylabel('|C_{max}| / N');
title('Composante geante temporelle - orbites aleatoires a inclinaison fixe');

figure;
plot(time_values, num_edges_time, 'LineWidth', 1.5);
grid on;
xlabel('Temps (s)');
ylabel('Nombre de liens');
title('Nombre de liens temporel - orbites aleatoires a inclinaison fixe');

fprintf('Modele orbital : orbites aleatoires a inclinaison fixe, i = %.1f deg, P = %d plans, N = %d satellites.\n', inc_deg, P, N);

%% Fonction locale orbites aleatoires a inclinaison fixe
function positions = walker_delta_positions(R, inc, Omega, u)
    x = R * (cos(Omega).*cos(u) - sin(Omega).*sin(u).*cos(inc));
    y = R * (sin(Omega).*cos(u) + cos(Omega).*sin(u).*cos(inc));
    z = R * (sin(u).*sin(inc));
    positions = [x y z];
end
