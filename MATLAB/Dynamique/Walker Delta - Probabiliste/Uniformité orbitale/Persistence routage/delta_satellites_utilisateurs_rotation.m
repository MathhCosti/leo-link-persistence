clear; clc; close all;

%% ============================================================
%  MODELE DELTA STOCHASTIQUE AVEC UTILISATEURS AU SOL
%
%  - Satellites : orbites circulaires, inclinaison commune,
%    RAAN et phases initiales aleatoires.
%  - Utilisateurs : fixes dans le repere terrestre, donc en rotation
%    dans le repere inertiel utilise pour les satellites.
%  - Ce code genere et anime uniquement les positions.
%% ============================================================

rng(1); % Reproductibilite ; commenter pour changer chaque execution

%% ============================================================
%  1. PARAMETRES PHYSIQUES
%% ============================================================

R_earth = 6371;          % Rayon terrestre [km]
h = 550;                 % Altitude orbitale [km]
R_orbit = R_earth + h;   % Rayon orbital [km]

mu = 398600;                         % Parametre gravitationnel [km^3/s^2]
omega_sat = sqrt(mu / R_orbit^3);    % Vitesse angulaire satellite [rad/s]

T_sidereal = 86164;                  % Jour sideral [s]
omega_earth = 2*pi / T_sidereal;     % Rotation terrestre [rad/s]

%% ============================================================
%  2. GENERATION DES SATELLITES
%% ============================================================

lambda_sat = 4e-7;                   % Intensite [satellites/km^2]
surface_orbital_sphere = 4*pi*R_orbit^2;

N_sat = poissrnd(lambda_sat * surface_orbital_sphere);

inc_deg = 58;                        % Inclinaison commune [degres]
inc = deg2rad(inc_deg);

% Un RAAN et une phase orbitale initiale independants par satellite.
Omega = 2*pi*rand(N_sat,1);
u0 = 2*pi*rand(N_sat,1);

fprintf('Nombre de satellites : %d\n', N_sat);
fprintf('Inclinaison commune  : %.1f deg\n', inc_deg);

%% ============================================================
%  3. GENERATION DES UTILISATEURS AU SOL
%% ============================================================

N_users = 2;

% Deux modes possibles :
%   'random' : utilisateurs uniformes sur la surface terrestre ;
%   'fixed'  : latitudes et longitudes imposees ci-dessous.
user_generation_mode = 'random';

switch lower(user_generation_mode)
    case 'random'
        % Uniformite surfacique sur la sphere terrestre :
        % longitude uniforme et sin(latitude) uniforme.
        user_lon0 = -pi + 2*pi*rand(N_users,1);
        user_lat = asin((2*rand(N_users,1)-1) * sin(inc));

    case 'fixed'
        % Exemple : Paris, Montreal, Tokyo et Sydney.
        user_lat_deg = [48.8566; 45.5019; 35.6762; -33.8688];
        user_lon_deg = [ 2.3522; -73.5674; 139.6503; 151.2093];

        user_lat = deg2rad(user_lat_deg);
        user_lon0 = deg2rad(user_lon_deg);
        N_users = numel(user_lat);

    otherwise
        error('Mode de generation utilisateur inconnu : %s', ...
            user_generation_mode);
end

fprintf('Nombre d''utilisateurs : %d\n', N_users);

%% ============================================================
%  4. PARAMETRES TEMPORELS
%% ============================================================

dt = 20;                  % Pas temporel [s]
Tmax = 6000;              % Duree totale [s]
time_values = 0:dt:Tmax;
Nt = numel(time_values);

% Stockage facultatif des trajectoires.
SatellitePositions = cell(Nt,1);
UserPositions = cell(Nt,1);

%% ============================================================
%  5. FIGURE 3D
%% ============================================================

figure('Color','w');
hold on;
axis equal;
grid on;
view(3);
rotate3d on;

xlabel('x [km]');
ylabel('y [km]');
zlabel('z [km]');

axis_limit = 1.08 * R_orbit;
xlim([-axis_limit axis_limit]);
ylim([-axis_limit axis_limit]);
zlim([-axis_limit axis_limit]);

% Terre.
[xe, ye, ze] = sphere(100);
earth_handle = surf(R_earth*xe, R_earth*ye, R_earth*ze, ...
    'FaceAlpha', 0.25, ...
    'EdgeColor', 'none', ...
    'FaceColor', [0.3 0.6 0.9]);

% Positions initiales.
sat_pos0 = walker_delta_positions(R_orbit, inc, Omega, u0);
user_pos0 = ground_user_positions(R_earth, user_lat, user_lon0);

sat_handle = scatter3(sat_pos0(:,1), sat_pos0(:,2), sat_pos0(:,3), ...
    24, 'filled', 'MarkerFaceColor', [0.85 0.2 0.1]);

user_handle = scatter3(user_pos0(:,1), user_pos0(:,2), user_pos0(:,3), ...
    45, 'filled', 'MarkerFaceColor', [0.1 0.65 0.2], ...
    'MarkerEdgeColor', 'k');

legend([earth_handle, sat_handle, user_handle], ...
    {'Terre', 'Satellites', 'Utilisateurs'}, ...
    'Location', 'bestoutside');

title_handle = title('');

%% ============================================================
%  6. EVOLUTION TEMPORELLE
%% ============================================================

for k = 1:Nt

    t = time_values(k);

    % Satellites : mouvement orbital dans le repere inertiel.
    u_t = mod(u0 + omega_sat*t, 2*pi);
    sat_positions_t = walker_delta_positions( ...
        R_orbit, inc, Omega, u_t);

    % Utilisateurs : longitude inertielle = longitude terrestre initiale
    % + rotation de la Terre.
    user_lon_t = mod(user_lon0 + omega_earth*t + pi, 2*pi) - pi;
    user_positions_t = ground_user_positions( ...
        R_earth, user_lat, user_lon_t);

    % Stockage.
    SatellitePositions{k} = sat_positions_t;
    UserPositions{k} = user_positions_t;

    % Mise a jour graphique.
    set(sat_handle, ...
        'XData', sat_positions_t(:,1), ...
        'YData', sat_positions_t(:,2), ...
        'ZData', sat_positions_t(:,3));

    set(user_handle, ...
        'XData', user_positions_t(:,1), ...
        'YData', user_positions_t(:,2), ...
        'ZData', user_positions_t(:,3));

    set(title_handle, 'String', sprintf([ ...
        'Delta stochastique + utilisateurs au sol | ' ...
        't = %.0f s | N_{sat} = %d | N_{users} = %d'], ...
        t, N_sat, N_users));

    drawnow;
    pause(0.03);
end

%% ============================================================
%  7. SAUVEGARDE
%% ============================================================

save('delta_satellites_utilisateurs_rotation.mat', ...
    'R_earth', 'h', 'R_orbit', 'mu', ...
    'omega_sat', 'omega_earth', ...
    'lambda_sat', 'N_sat', 'inc_deg', 'inc', ...
    'Omega', 'u0', ...
    'N_users', 'user_lat', 'user_lon0', ...
    'dt', 'Tmax', 'time_values', ...
    'SatellitePositions', 'UserPositions');

fprintf(['Donnees sauvegardees dans ' ...
    'delta_satellites_utilisateurs_rotation.mat\n']);

%% ============================================================
%  FONCTIONS LOCALES
%% ============================================================

function positions = walker_delta_positions(R, inc, Omega, u)
    % Positions ECI de satellites sur des orbites circulaires ayant
    % une inclinaison commune inc.

    x = R .* (cos(Omega).*cos(u) ...
        - sin(Omega).*sin(u).*cos(inc));

    y = R .* (sin(Omega).*cos(u) ...
        + cos(Omega).*sin(u).*cos(inc));

    z = R .* sin(u).*sin(inc);

    positions = [x, y, z];
end

function positions = ground_user_positions(R_earth, latitude, longitude)
    % Positions ECI des utilisateurs a partir de leur latitude fixe
    % et de leur longitude inertielle courante.

    x = R_earth .* cos(latitude).*cos(longitude);
    y = R_earth .* cos(latitude).*sin(longitude);
    z = R_earth .* sin(latitude);

    positions = [x, y, z];
end
