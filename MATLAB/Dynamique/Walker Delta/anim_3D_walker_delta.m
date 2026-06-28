clear; clc; close all;

%% ============================================================
%  ANIMATION 3D D'UN RESEAU LEO
%  Version Walker-Delta
%% ============================================================

%% Parametres physiques
R_earth = 6371;      % km
h = 550;             % km
R = R_earth + h;     % rayon orbital

mu = 398600;              % km^3/s^2
omega = sqrt(mu / R^3);   % vitesse angulaire orbitale rad/s

%% Parametres du processus de Poisson
lambda = 5e-7;       % satellites / km^2
surface_sphere = 4*pi*R^2;

N = poissrnd(lambda * surface_sphere);
fprintf('Nombre de satellites generes : N = %d\n', N);

%% Parametres Walker-Delta
inc_deg = 53;                  % inclinaison commune, en degres
inc = deg2rad(inc_deg);        % radians
P = max(1, round(sqrt(N)));     % nombre de plans orbitaux, modifiable

% Walker-Delta : RAAN repartis sur 360 degres.
Omega_planes = 2*pi*(0:P-1)'/P;
plane_id = mod((0:N-1)', P) + 1;
Omega = Omega_planes(plane_id);

% Position initiale dans chaque plan.
u0 = 2*pi*rand(N,1);

positions0 = walker_delta_positions(R, inc, Omega, u0);

%% Parametres des liens et de l'animation
dmax = 1500;     % km
dt = 30;         % pas temporel en secondes
Tmax = 6000;     % duree totale de simulation
time_values = 0:dt:Tmax;

%% Creation de la figure
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

%% Option : affichage des plans orbitaux
for p = 1:P
    uu = linspace(0, 2*pi, 300)';
    OO = Omega_planes(p) * ones(size(uu));
    plane_curve = walker_delta_positions(R, inc, OO, uu);
    plot3(plane_curve(:,1), plane_curve(:,2), plane_curve(:,3), ':', ...
        'LineWidth', 0.5);
end

%% Boucle d'animation
for k = 1:length(time_values)

    t = time_values(k);

    %% Mouvement orbital Walker-Delta
    u_t = u0 + omega*t;
    positions_t = walker_delta_positions(R, inc, Omega, u_t);

    x_t = positions_t(:,1);
    y_t = positions_t(:,2);
    z_t = positions_t(:,3);

    %% Construction du graphe a l'instant t
    D = squareform(pdist(positions_t));
    A = (D <= dmax) & (D > 0);

    [row, col] = find(triu(A, 1));
    E = length(row);

    Xlinks = NaN(3*E, 1);
    Ylinks = NaN(3*E, 1);
    Zlinks = NaN(3*E, 1);

    Xlinks(1:3:end) = x_t(row);
    Xlinks(2:3:end) = x_t(col);

    Ylinks(1:3:end) = y_t(row);
    Ylinks(2:3:end) = y_t(col);

    Zlinks(1:3:end) = z_t(row);
    Zlinks(2:3:end) = z_t(col);

    %% Mise a jour graphique
    set(sat_handle, ...
        'XData', x_t, ...
        'YData', y_t, ...
        'ZData', z_t);

    set(link_handle, ...
        'XData', Xlinks, ...
        'YData', Ylinks, ...
        'ZData', Zlinks);

    set(title_handle, 'String', ...
        sprintf('Graphe LEO dynamique Walker-Delta | i = %.1f deg | t = %.0f s | N = %d | E = %d', ...
        inc_deg, t, N, E));

    drawnow;
    pause(0.5);
end

hold off;

%% Fonction locale Walker-Delta
function positions = walker_delta_positions(R, inc, Omega, u)
    x = R * (cos(Omega).*cos(u) - sin(Omega).*sin(u).*cos(inc));
    y = R * (sin(Omega).*cos(u) + cos(Omega).*sin(u).*cos(inc));
    z = R * (sin(u).*sin(inc));
    positions = [x y z];
end
