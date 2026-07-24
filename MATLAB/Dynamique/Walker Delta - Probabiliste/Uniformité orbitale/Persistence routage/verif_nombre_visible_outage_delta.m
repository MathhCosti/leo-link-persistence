clear; clc; close all;

%% ============================================================
% VERIFICATION EMPIRIQUE / THEORIQUE : NOMBRE DE SATELLITES VISIBLES
% ET PROBABILITE D'OUTAGE
%
% Modele :
% - orbites circulaires d'inclinaison commune i ;
% - Omega et u0 uniformes ;
% - nombre de satellites Poisson de moyenne lambda*4*pi*R^2 ;
% - utilisateurs fixes sur la Terre tournante ;
% - visibilite si elevation >= elevation_min.
%
% Comparaisons :
% 1) moyenne empirique du nombre visible ;
% 2) formule exacte obtenue par integration sur u et la longitude ;
% 3) approximation locale par densite surfacique ;
% 4) outage empirique et approximation exp(-E[N_vis]).
%% ============================================================

rng(1);

%% Parametres physiques et orbitaux
R_earth = 6371;                  % [km]
h = 550;                         % [km]
R = R_earth + h;                 % [km]
mu = 398600;                     % [km^3/s^2]
omega_sat = sqrt(mu/R^3);        % [rad/s]
omega_earth = 2*pi/86164;        % [rad/s]

lambda = 4e-7;                   % [satellites/km^2 sur sphere orbitale]
N_mean = lambda*4*pi*R^2;

inc_deg = 58;
inc = deg2rad(inc_deg);

elevation_min_deg = 20;
elevation_min = deg2rad(elevation_min_deg);

%% Utilisateurs : latitudes imposees pour voir l'effet spatial
user_lat_deg = [-50 -35 -20 0 20 35 50].';
user_lat = deg2rad(user_lat_deg);
N_users = numel(user_lat);
user_lon0 = 2*pi*rand(N_users,1)-pi;

%% Monte-Carlo temporel
N_realizations = 80;
dt = 30;                         % [s]
Tmax = 12000;                    % [s]
time_values = (0:dt:Tmax-dt).';
Nt = numel(time_values);

sum_visible = zeros(N_users,1);
sum_outage = zeros(N_users,1);
num_samples = N_realizations*Nt;

for r = 1:N_realizations
    N_sat = poissrnd(N_mean);
    Omega = 2*pi*rand(N_sat,1);
    u0 = 2*pi*rand(N_sat,1);

    for k = 1:Nt
        t = time_values(k);

        u_t = mod(u0 + omega_sat*t,2*pi);
        sat_pos = walker_delta_positions(R,inc,Omega,u_t);

        lon_t = mod(user_lon0 + omega_earth*t + pi,2*pi)-pi;
        user_pos = ground_user_positions(R_earth,user_lat,lon_t);

        for q = 1:N_users
            rho = sat_pos-user_pos(q,:);
            dist = sqrt(sum(rho.^2,2));
            zenith = user_pos(q,:)/R_earth;
            sin_el = (rho*zenith.') ./ dist;
            el = asin(max(-1,min(1,sin_el)));

            nvis = nnz(el >= elevation_min);
            sum_visible(q) = sum_visible(q)+nvis;
            sum_outage(q) = sum_outage(q)+(nvis == 0);
        end
    end
end

MeanVisible_emp = sum_visible/num_samples;
Outage_emp = sum_outage/num_samples;

%% Theorie geometrique
psi_max = acos((R_earth/R)*cos(elevation_min))-elevation_min;

p_vis_exact = zeros(N_users,1);
MeanVisible_exact = zeros(N_users,1);
MeanVisible_local = zeros(N_users,1);
Outage_exact = zeros(N_users,1);
Outage_local = zeros(N_users,1);

for q = 1:N_users
    phi_u = user_lat(q);

    % Integration directe sur l'argument de latitude u. Cette ecriture
    % evite la singularite de la densite de latitude aux bords +/-i.
    integrand = @(u) visible_longitude_fraction(u,phi_u,inc,psi_max);
    p_vis_exact(q) = integral(integrand,0,2*pi, ...
        'RelTol',1e-9,'AbsTol',1e-11)/(2*pi);

    MeanVisible_exact(q) = N_mean*p_vis_exact(q);

    % Approximation locale, valable loin des bords de bande.
    denom = sqrt(max(sin(inc)^2-sin(phi_u)^2,eps));
    MeanVisible_local(q) = ...
        N_mean*(1-cos(psi_max))/(pi*denom);

    Outage_exact(q) = exp(-MeanVisible_exact(q));
    Outage_local(q) = exp(-MeanVisible_local(q));
end

%% Tableau de comparaison
Results = table(user_lat_deg,MeanVisible_emp,MeanVisible_exact, ...
    MeanVisible_local,Outage_emp,Outage_exact,Outage_local, ...
    'VariableNames',{'Latitude_deg','MeanVisible_emp', ...
    'MeanVisible_exact','MeanVisible_local','Outage_emp', ...
    'Outage_exact','Outage_local'});

disp(Results);

fprintf('\nParametres : N moyen = %.2f, i = %.1f deg, elevation min = %.1f deg\n', ...
    N_mean,inc_deg,elevation_min_deg);
fprintf('psi_max = %.4f rad = %.2f deg\n',psi_max,rad2deg(psi_max));

%% Erreurs relatives
err_visible_exact = abs(MeanVisible_emp-MeanVisible_exact) ...
    ./ max(MeanVisible_emp,eps);
err_outage_exact = abs(Outage_emp-Outage_exact) ...
    ./ max(Outage_emp,eps);

fprintf('\nErreur relative moyenne E[N_vis], theorie exacte : %.2f %%\n', ...
    100*mean(err_visible_exact));
fprintf('Erreur relative moyenne outage, exp(-E[N_vis]) : %.2f %%\n', ...
    100*mean(err_outage_exact));

%% Graphiques
figure;
plot(user_lat_deg,MeanVisible_emp,'o-','LineWidth',1.5, ...
    'DisplayName','Empirique');
hold on;
plot(user_lat_deg,MeanVisible_exact,'s--','LineWidth',1.5, ...
    'DisplayName','Theorie exacte');
plot(user_lat_deg,MeanVisible_local,'d:','LineWidth',1.5, ...
    'DisplayName','Approximation locale');
grid on;
xlabel('Latitude utilisateur (deg)');
ylabel('Nombre moyen de satellites visibles');
title('Nombre moyen visible : empirique et theorie');
legend('Location','best');

figure;
plot(user_lat_deg,Outage_emp,'o-','LineWidth',1.5, ...
    'DisplayName','Empirique');
hold on;
plot(user_lat_deg,Outage_exact,'s--','LineWidth',1.5, ...
    'DisplayName','exp(-m_{vis exact})');
plot(user_lat_deg,Outage_local,'d:','LineWidth',1.5, ...
    'DisplayName','exp(-m_{vis local})');
grid on;
xlabel('Latitude utilisateur (deg)');
ylabel('Fraction d''outage');
title('Outage : empirique et approximation de Poisson');
legend('Location','best');

writetable(Results,'verification_nombre_visible_outage.csv');
save('verification_nombre_visible_outage.mat');

%% Fonctions locales
function positions = walker_delta_positions(R,inc,Omega,u)
    x = R*(cos(Omega).*cos(u)-sin(Omega).*sin(u).*cos(inc));
    y = R*(sin(Omega).*cos(u)+cos(Omega).*sin(u).*cos(inc));
    z = R*(sin(u).*sin(inc));
    positions = [x y z];
end

function positions = ground_user_positions(R,lat,lon)
    positions = [R*cos(lat).*cos(lon), ...
                 R*cos(lat).*sin(lon), ...
                 R*sin(lat)];
end

function fraction = visible_longitude_fraction(u,phi_u,inc,psi_max)
    % Pour un u donne, la latitude du satellite est fixee tandis que sa
    % longitude est uniforme grace au RAAN uniforme.
    phi_s = asin(sin(inc).*sin(u));

    denominator = cos(phi_u).*cos(phi_s);
    numerator = cos(psi_max)-sin(phi_u).*sin(phi_s);

    q = numerator./denominator;
    fraction = zeros(size(u));

    full = q <= -1;
    none = q >= 1;
    partial = ~(full | none);

    fraction(full) = 1;
    fraction(none) = 0;
    fraction(partial) = acos(q(partial))/pi;

    % Traitement robuste des poles eventuels.
    singular = abs(denominator) < 1e-14;
    if any(singular)
        angular_distance = abs(phi_s(singular)-phi_u);
        fraction(singular) = double(angular_distance <= psi_max);
    end
end
