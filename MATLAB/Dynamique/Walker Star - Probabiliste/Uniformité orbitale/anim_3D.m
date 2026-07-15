clear; clc; close all;

%% Paramètres physiques
R_earth = 6371;      % km
h = 550;             % km
R = R_earth + h;     % rayon orbital

mu = 398600;              % km^3/s^2
omega = sqrt(mu / R^3);   % vitesse angulaire orbitale rad/s

%% Paramètres du processus de Poisson
lambda = 4e-7;       % satellites / km^2
surface_sphere = 4*pi*R^2;

N = poissrnd(lambda * surface_sphere);

%% Génération uniforme sur l'orbite choisie (Walker Star)
% Omega choisit le plan orbital méridien. Comme les plans Omega et
% Omega + pi décrivent le même grand cercle, Omega est tiré sur [0, pi).
Omega = 2*pi * rand(N,1);

% u0 est l'argument orbital initial. Un tirage uniforme sur [0, 2*pi)
% rend chaque satellite uniforme en longueur d'arc sur son orbite.
u0 = 2*pi * rand(N,1);

% Paramétrisation du grand cercle contenu dans le plan défini par Omega
x = R * cos(u0) .* cos(Omega);
y = R * cos(u0) .* sin(Omega);
z = R * sin(u0);

positions0 = [x y z];

%% Sens de rotation défini par deux demi-espaces séparés par y = 0
% Le plan y = 0 contient l'axe des pôles et sépare réellement l'espace
% en deux moitiés :
%   y0 >= 0  -> sens orbital +1 ;
%   y0 <  0  -> sens orbital -1.
%
% IMPORTANT : Omega est tiré sur [0, 2*pi). Ainsi, contrairement au cas
% Omega dans [0, pi), le signe de y0 n'est pas artificiellement corrélé
% au signe de cos(u0). Les deux moitiés ne sont donc pas toutes dirigées
% vers le même pôle au temps initial.
%
% Le signe est fixé à t = 0 et reste constant pendant toute la simulation.

y0 = y;
rotation_sign = ones(N,1);
rotation_sign(y0 < 0) = -1;

fprintf('Sens + : %d satellites | Sens - : %d satellites\n', ...
    nnz(rotation_sign == 1), nnz(rotation_sign == -1));

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
colors = zeros(N,3);
colors(rotation_sign == 1,:) = repmat([0 0 1], nnz(rotation_sign == 1), 1);
colors(rotation_sign == -1,:) = repmat([1 0 0], nnz(rotation_sign == -1), 1);

sat_handle = scatter3(positions0(:,1), positions0(:,2), positions0(:,3), ...
    25, colors, 'filled');

% Liens : un seul objet graphique optimisé
link_handle = plot3(NaN, NaN, NaN, 'k-', 'LineWidth', 0.5);

%% Boucle d'animation
for k = 1:length(time_values)

    t = time_values(k);

    %% Mouvement orbital avec pôles Nord/Sud communs
    % Le plan orbital Omega reste constant et la phase évolue
    % déterministement à la vitesse angulaire orbitale omega.
    u_t = u0 + rotation_sign * omega * t;

    % Position sur le grand cercle orbital choisi
    x_t = R * cos(u_t) .* cos(Omega);
    y_t = R * cos(u_t) .* sin(Omega);
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
        sprintf('Graphe LEO dynamique à pôles communs | t = %.0f s | N = %d | E = %d', ...
        t, N, E));

    drawnow;
    pause(0.5);
end

hold off;