clear; clc; close all;

%% ============================================================
% PROBABILITE D'OUTAGE
% Comparaison empirique / theorie exacte / approximation locale
%% ============================================================

rng(1);

%% Parametres physiques et orbitaux
R_earth = 6371;                  % [km]
h = 550;                         % [km]
R = R_earth + h;                 % [km]
mu = 398600;                     % [km^3/s^2]

omega_sat = sqrt(mu/R^3);        % [rad/s]
omega_earth = 2*pi/86164;        % [rad/s]

lambda = 4e-7;                   % [satellites/km^2]
N_mean = lambda*4*pi*R^2;

inc_deg = 58;
inc = deg2rad(inc_deg);

elevation_min_deg = 20;
elevation_min = deg2rad(elevation_min_deg);

%% Utilisateurs
user_lat_deg = [-50 -35 -20 0 20 35 50].';
user_lat = deg2rad(user_lat_deg);
N_users = numel(user_lat);

user_lon0 = 2*pi*rand(N_users,1)-pi;

%% Monte-Carlo
N_realizations = 80;
dt = 30;                         % [s]
Tmax = 12000;                    % [s]

time_values = (0:dt:Tmax-dt).';
Nt = numel(time_values);

sum_outage = zeros(N_users,1);
num_samples = N_realizations*Nt;

for r = 1:N_realizations

    N_sat = poissrnd(N_mean);

    Omega = 2*pi*rand(N_sat,1);
    u0 = 2*pi*rand(N_sat,1);

    for k = 1:Nt

        t = time_values(k);

        %% Positions satellites Walker Delta
        u_t = mod(u0 + omega_sat*t,2*pi);

        x_sat = R*(cos(Omega).*cos(u_t) ...
              - sin(Omega).*sin(u_t).*cos(inc));

        y_sat = R*(sin(Omega).*cos(u_t) ...
              + cos(Omega).*sin(u_t).*cos(inc));

        z_sat = R*(sin(u_t).*sin(inc));

        sat_pos = [x_sat y_sat z_sat];

        %% Positions utilisateurs
        lon_t = mod(user_lon0 + omega_earth*t + pi,2*pi)-pi;

        x_user = R_earth*cos(user_lat).*cos(lon_t);
        y_user = R_earth*cos(user_lat).*sin(lon_t);
        z_user = R_earth*sin(user_lat);

        user_pos = [x_user y_user z_user];

        %% Outage empirique
        for q = 1:N_users

            rho = sat_pos - user_pos(q,:);
            dist = sqrt(sum(rho.^2,2));

            zenith = user_pos(q,:)/R_earth;

            sin_el = (rho*zenith.') ./ dist;
            el = asin(max(-1,min(1,sin_el)));

            nvis = nnz(el >= elevation_min);

            sum_outage(q) = ...
                sum_outage(q) + (nvis == 0);
        end
    end
end

Outage_emp = sum_outage/num_samples;

%% ============================================================
% Theorie
%% ============================================================

psi_max = acos((R_earth/R)*cos(elevation_min)) - elevation_min;

MeanVisible_exact = zeros(N_users,1);
MeanVisible_local = zeros(N_users,1);

Outage_exact = zeros(N_users,1);
Outage_local = zeros(N_users,1);

Nu = 20000;
u_grid = linspace(0,2*pi,Nu);

for q = 1:N_users

    phi_u = user_lat(q);

    %% Latitude satellite
    phi_s = asin(sin(inc).*sin(u_grid));

    denominator = cos(phi_u).*cos(phi_s);
    numerator = cos(psi_max)-sin(phi_u).*sin(phi_s);

    qq = numerator ./ denominator;

    fraction = zeros(size(u_grid));

    full = qq <= -1;
    none = qq >= 1;
    partial = ~(full | none);

    fraction(full) = 1;
    fraction(none) = 0;
    fraction(partial) = acos(qq(partial))/pi;

    singular = abs(denominator) < 1e-14;

    if any(singular)
        angular_distance = abs(phi_s(singular)-phi_u);
        fraction(singular) = double(angular_distance <= psi_max);
    end

    %% Nombre moyen visible exact
    p_vis_exact = trapz(u_grid,fraction)/(2*pi);
    MeanVisible_exact(q) = N_mean*p_vis_exact;

    %% Nombre moyen visible local
    denom_local = sqrt(max( ...
        sin(inc)^2 - sin(phi_u)^2, eps));

    MeanVisible_local(q) = ...
        N_mean*(1-cos(psi_max))/(pi*denom_local);

    %% Probabilite d'outage
    Outage_exact(q) = exp(-MeanVisible_exact(q));
    Outage_local(q) = exp(-MeanVisible_local(q));
end

%% Tableau
Results = table( ...
    user_lat_deg, ...
    Outage_emp, ...
    Outage_exact, ...
    Outage_local, ...
    'VariableNames', { ...
    'Latitude_deg', ...
    'Outage_emp', ...
    'Outage_exact', ...
    'Outage_local'});

disp(Results);

fprintf('\nN moyen = %.2f\n',N_mean);
fprintf('Inclinaison = %.1f deg\n',inc_deg);
fprintf('Elevation minimale = %.1f deg\n',elevation_min_deg);
fprintf('psi_max = %.2f deg\n',rad2deg(psi_max));

%% Erreur relative
err_outage_exact = ...
    abs(Outage_emp-Outage_exact) ...
    ./ max(Outage_emp,eps);

fprintf('\nErreur relative moyenne outage : %.2f %%\n', ...
    100*mean(err_outage_exact));

%% Graphique
figure;

plot(user_lat_deg,Outage_emp, ...
    'o-','LineWidth',1.5, ...
    'DisplayName','Empirique');

hold on;

plot(user_lat_deg,Outage_exact, ...
    's--','LineWidth',1.5, ...
    'DisplayName','Theorie exacte');

plot(user_lat_deg,Outage_local, ...
    'd:','LineWidth',1.5, ...
    'DisplayName','Approximation locale');

grid on;
xlabel('Latitude utilisateur (deg)');
ylabel('Probabilite d''outage');
title('Probabilite d''outage');
legend('Location','best');

%% Sauvegarde
save('p_outage_results.mat');
