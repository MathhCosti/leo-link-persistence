clear; clc; close all;

%% ============================================================
%  ANIMATION 3D D'UN RESEAU LEO
%  Mouvement avec vecteurs aléatoires tangentiels
%
%  Idée :
%  - chaque satellite part d'une position uniforme sur la sphère ;
%  - on lui associe une direction aléatoire tangente à la sphère ;
%  - il se déplace ensuite sur le grand cercle défini par cette direction.
%
%  Avantage : les satellites restent exactement à distance R du centre.
%% ============================================================

%% Paramètres physiques
R_earth = 6371;      % km
h = 550;             % km
R = R_earth + h;     % rayon orbital

mu = 398600;              % km^3/s^2
omega = sqrt(mu / R^3);   % vitesse angulaire rad/s

%% Paramètres du processus de Poisson
lambda = 10e-7;       % satellites / km^2
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

% Projection sur le plan tangent :
% v = a - (a.r0) r0
v = a - sum(a .* r0, 2) .* r0;

% Normalisation des directions tangentes
v = v ./ vecnorm(v, 2, 2);

% Optionnel : sens aléatoire +1 / -1
sens = sign(rand(N,1) - 0.5);
sens(sens == 0) = 1;
v = sens .* v;

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
sat_handle = scatter3(positions0(:,1), positions0(:,2), positions0(:,3), ...
    25, 'filled');

link_handle = plot3(NaN, NaN, NaN, 'k-', 'LineWidth', 0.5);

%% Boucle d'animation
for k = 1:length(time_values)

    t = time_values(k);

    %% Mouvement sur grand cercle avec direction tangentielle aléatoire
    % Formule géodésique sur la sphère :
    % r(t) = R [ r0 cos(omega t) + v sin(omega t) ]
    positions_t = R * (r0 * cos(omega*t) + v * sin(omega*t));

    x_t = positions_t(:,1);
    y_t = positions_t(:,2);
    z_t = positions_t(:,3);

    %% Construction du graphe à l'instant t
    D = squareform(pdist(positions_t));
    A = (D <= dmax) & (D > 0);

    [row, col] = find(triu(A, 1));
    E = length(row);

    %% Construction optimisée des segments de liens
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
        sprintf('Graphe LEO dynamique avec vecteurs aléatoires | t = %.0f s | N = %d | E = %d', ...
        t, N, E));

    drawnow;
    pause(0.5);
end

hold off;
