clear; clc; close all;

%% Paramètres physiques
R_earth = 6371;      % km
h = 550;             % km
R = R_earth + h;     % rayon orbital

mu = 398600;              % km^3/s^2
omega = sqrt(mu / R^3);   % vitesse angulaire orbitale rad/s

%% Paramètres du processus de Poisson
lambda = 5e-7;       % satellites / km^2
surface_sphere = 4*pi*R^2;

N = poissrnd(lambda * surface_sphere);

%% Génération uniforme des positions initiales sur la sphère
u = rand(N,1);
phi = 2*pi*rand(N,1);        % longitude / plan orbital
theta = acos(1 - 2*u);       % colatitude, tirage uniforme correct sur sphère

x = R * sin(theta) .* cos(phi);
y = R * sin(theta) .* sin(phi);
z = R * cos(theta);

positions0 = [x y z];

%% Sens de rotation défini par un plan séparateur passant par les pôles
% +1 : sens croissant de theta
% -1 : sens décroissant de theta
% Plan séparateur choisi : y = 0
% Il contient l'axe z, donc les pôles Nord/Sud.

rotation_sign = ones(N, 1);

% Une moitié de la sphère tourne dans un sens
rotation_sign(y >= 0) = 1;

% L'autre moitié tourne dans le sens opposé
rotation_sign(y < 0) = -1;

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

    %% Mouvement orbital avec pôles Nord/Sud communs
    % Le satellite reste dans son plan orbital : phi constant
    phi_t = phi;

    % Le satellite avance sur son grand cercle méridien
    theta_t = theta + rotation_sign * omega * t;

    % Coordonnées sphériques -> cartésiennes
    x_t = R * sin(theta_t) .* cos(phi_t);
    y_t = R * sin(theta_t) .* sin(phi_t);
    z_t = R * cos(theta_t);

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
        sprintf('Graphe LEO dynamique à pôles communs | t = %.0f s | N = %d | E = %d', ...
        t, N, E));

    drawnow;
    pause(0.5);
end

hold off;