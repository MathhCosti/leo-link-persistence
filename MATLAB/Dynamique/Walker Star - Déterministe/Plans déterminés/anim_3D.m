clear; clc; close all;

%% Paramètres physiques
R_earth = 6371;      % km
h = 550;             % km
R = R_earth + h;     % rayon orbital

mu = 398600;              % km^3/s^2
omega = sqrt(mu / R^3);   % vitesse angulaire orbitale rad/s

%% Paramètres de densité et de constellation Walker Star
lambda = 5e-7;       % satellites / km^2
surface_sphere = 4*pi*R^2;

% Walker Star : P plans orbitaux polaires répartis sur 180 degrés.
N_target = poissrnd(lambda * surface_sphere);
P_walker = max(1, round(sqrt(N_target)));       % nombre de plans orbitaux
F_walker = 1;        % facteur de phasage Walker

S_walker = max(1, round(N_target/P_walker));
N = P_walker * S_walker;

fprintf('Nombre cible Poisson : N_cible = %d\n', N_target);
fprintf('Walker Star : P = %d plans, S = %d satellites/plan, N = %d satellites\n', ...
    P_walker, S_walker, N);

%% Génération Walker Star des positions initiales
Omega_planes = (0:P_walker-1)' * pi/P_walker;   % Walker Star : 180 degrés

Omega = zeros(N,1);
u0 = zeros(N,1);
idx_sat = 1;

for p_idx = 0:P_walker-1
    Omega_p = Omega_planes(p_idx+1);
    phase_p = 2*pi*F_walker*p_idx/N;

    for s_idx = 0:S_walker-1
        Omega(idx_sat) = Omega_p;
        u0(idx_sat) = 2*pi*s_idx/S_walker + phase_p;
        idx_sat = idx_sat + 1;
    end
end

x = R * cos(Omega) .* cos(u0);
y = R * sin(Omega) .* cos(u0);
z = R * sin(u0);

positions0 = [x y z];

%% Paramètres des liens et de l'animation
dmax = 1500;     % km
dt = 30;         % pas temporel en secondes
Tmax = 6000;     % durée totale de simulation
time_values = 0:dt:Tmax;

%% Création de la figure
figure;
hold on;

axis equal;
grid on;
xlabel('x (km)');
ylabel('y (km)');
zlabel('z (km)');
view(3);
rotate3d on;

title_handle = title('');

%% Initialisation des objets graphiques

% Satellites
sat_handle = scatter3(positions0(:,1), positions0(:,2), positions0(:,3), ...
    25, 'filled');

% Liens : un seul objet graphique optimisé
link_handle = plot3(NaN, NaN, NaN, 'k-', 'LineWidth', 0.5);

%% Boucle d'animation
for k = 1:length(time_values)

    t = time_values(k);

    %% Mouvement orbital Walker Star
    % Chaque satellite reste dans son plan orbital Omega.
    % Sa position angulaire dans le plan évolue à vitesse omega.
    u_t = u0 + omega*t;

    x_t = R * cos(Omega) .* cos(u_t);
    y_t = R * sin(Omega) .* cos(u_t);
    z_t = R * sin(u_t);

    positions_t = [x_t y_t z_t];

    %% Construction du graphe à l'instant t
    D = squareform(pdist(positions_t));
    A = (D <= dmax) & (D > 0);

    % Liste des liens
    [row, col] = find(triu(A, 1));
    E = length(row);

    % Construction optimisée des segments
    Xlinks = NaN(3*E, 1);
    Ylinks = NaN(3*E, 1);
    Zlinks = NaN(3*E, 1);

    Xlinks(1:3:end) = x_t(row);
    Xlinks(2:3:end) = x_t(col);

    Ylinks(1:3:end) = y_t(row);
    Ylinks(2:3:end) = y_t(col);

    Zlinks(1:3:end) = z_t(row);
    Zlinks(2:3:end) = z_t(col);

    %% Mise à jour graphique
    set(sat_handle, ...
        'XData', x_t, ...
        'YData', y_t, ...
        'ZData', z_t);

    set(link_handle, ...
        'XData', Xlinks, ...
        'YData', Ylinks, ...
        'ZData', Zlinks);

    set(title_handle, 'String', ...
        sprintf('Walker Star | P = %d | S = %d | t = %.0f s | N = %d | E = %d', ...
        P_walker, S_walker, t, N, E));

    drawnow;
    pause(0.5);
end

hold off;