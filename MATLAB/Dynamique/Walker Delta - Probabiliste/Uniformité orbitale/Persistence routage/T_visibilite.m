clear; clc; close all;

%% ============================================================
% VERIFICATION EMPIRIQUE / THEORIQUE : DUREE DE VISIBILITE GSL
%
% On suit chaque paire utilisateur-satellite et on extrait chaque episode
% continu durant lequel elevation >= elevation_min.
%
% Theorie testee :
%   T_vis ~= pi*psi_max/(2*nu_rel)
%
% Deux vitesses sont comparees :
% - approximation simple nu_rel = sqrt(omega_sat^2 + omega_earth^2*cos(phi_u)^2) ;
% - vitesse relative angulaire moyenne mesuree numeriquement aux instants
%   visibles, pour isoler l'erreur due a l'approximation de corde.
%% ============================================================

rng(2);

R_earth = 6371;
h = 550;
R = R_earth+h;
mu = 398600;
omega_sat = sqrt(mu/R^3);
omega_earth = 2*pi/86164;

lambda = 4e-7;
N_mean = lambda*4*pi*R^2;

inc_deg = 58;
inc = deg2rad(inc_deg);

elevation_min_deg = 20;
elevation_min = deg2rad(elevation_min_deg);
psi_max = acos((R_earth/R)*cos(elevation_min))-elevation_min;

user_lat_deg = [-45 -25 0 25 45].';
user_lat = deg2rad(user_lat_deg);
N_users = numel(user_lat);
user_lon0 = 2*pi*rand(N_users,1)-pi;

N_realizations = 40;
dt = 2;                         % fin pour mesurer les episodes [s]
Tmax = 12000;
time_values = (0:dt:Tmax).';
Nt = numel(time_values);

all_durations = cell(N_users,1);
sum_rel_speed = zeros(N_users,1);
count_rel_speed = zeros(N_users,1);

for r = 1:N_realizations
    N_sat = poissrnd(N_mean);
    Omega = 2*pi*rand(N_sat,1);
    u0 = 2*pi*rand(N_sat,1);

    visible_prev = false(N_sat,N_users);
    episode_start = zeros(N_sat,N_users);

    for k = 1:Nt
        t = time_values(k);
        u_t = mod(u0+omega_sat*t,2*pi);
        sat_pos = walker_delta_positions(R,inc,Omega,u_t);

        lon_t = mod(user_lon0+omega_earth*t+pi,2*pi)-pi;
        user_pos = ground_user_positions(R_earth,user_lat,lon_t);

        visible_now = false(N_sat,N_users);

        for q = 1:N_users
            rho = sat_pos-user_pos(q,:);
            dist = sqrt(sum(rho.^2,2));
            zenith = user_pos(q,:)/R_earth;
            sin_el = (rho*zenith.') ./ dist;
            el = asin(max(-1,min(1,sin_el)));
            visible_now(:,q) = el >= elevation_min;

            % Vitesse relative angulaire instantanee des directions radiales.
            sat_vel = walker_delta_velocities(R,inc,Omega,u_t,omega_sat);
            user_vel = omega_earth*[-user_pos(q,2),user_pos(q,1),0];

            sat_unit = sat_pos/R;
            user_unit = user_pos(q,:)/R_earth;
            sat_unit_dot = sat_vel/R;
            user_unit_dot = user_vel/R_earth;
            relative_tangent = sat_unit_dot-user_unit_dot;
            nu_inst = sqrt(sum(relative_tangent.^2,2));

            idx_vis = visible_now(:,q);
            sum_rel_speed(q) = sum_rel_speed(q)+sum(nu_inst(idx_vis));
            count_rel_speed(q) = count_rel_speed(q)+nnz(idx_vis);
        end

        if k == 1
            episode_start(visible_now) = t;
        else
            entering = visible_now & ~visible_prev;
            leaving = ~visible_now & visible_prev;
            episode_start(entering) = t;

            [sat_ids,user_ids] = find(leaving);
            for j = 1:numel(sat_ids)
                q = user_ids(j);
                duration = t-episode_start(sat_ids(j),q);
                if duration > 0
                    all_durations{q}(end+1,1) = duration; %#ok<SAGROW>
                end
            end
        end

        visible_prev = visible_now;
    end

    % Episodes encore ouverts a la fin : exclus, car tronques par Tmax.
end

MeanDuration_emp = NaN(N_users,1);
MinDuration_emp = NaN(N_users,1);
MaxDuration_emp = NaN(N_users,1);
NumberEpisodes = zeros(N_users,1);

for q = 1:N_users
    d = all_durations{q};
    NumberEpisodes(q) = numel(d);
    if ~isempty(d)
        MeanDuration_emp(q) = mean(d);
        MinDuration_emp(q) = min(d);
        MaxDuration_emp(q) = max(d);
    end
end

nu_simple = sqrt(omega_sat^2 + ...
    (omega_earth*cos(user_lat)).^2);
nu_measured = sum_rel_speed ./ max(count_rel_speed,1);

MeanDuration_theory_simple = pi*psi_max ./ (2*nu_simple);
MeanDuration_theory_measuredSpeed = pi*psi_max ./ (2*nu_measured);

Results = table(user_lat_deg,NumberEpisodes,MeanDuration_emp, ...
    MinDuration_emp,MaxDuration_emp,nu_simple,nu_measured, ...
    MeanDuration_theory_simple,MeanDuration_theory_measuredSpeed, ...
    'VariableNames',{'Latitude_deg','NumberEpisodes','MeanDuration_emp_s', ...
    'MinDuration_emp_s','MaxDuration_emp_s','NuSimple_rad_s', ...
    'NuMeasured_rad_s','TheorySimple_s','TheoryMeasuredSpeed_s'});

disp(Results);

figure;
plot(user_lat_deg,MeanDuration_emp,'o-','LineWidth',1.5, ...
    'DisplayName','Empirique');
hold on;
plot(user_lat_deg,MeanDuration_theory_simple,'s--','LineWidth',1.5, ...
    'DisplayName','Théorie');
plot(user_lat_deg,MeanDuration_theory_measuredSpeed,'d:','LineWidth',1.5, ...
    'DisplayName','Corde + vitesse empirique');
grid on;
xlabel('Latitude utilisateur (deg)');
ylabel('Duree moyenne de visibilite (s)');
title('Duree GSL visible : empirique et approximation');
legend('Location','best');

% writetable(Results,'verification_duree_visibilite.csv');
save('T_visibilite_results.mat');

%% Fonctions
function positions = walker_delta_positions(R,inc,Omega,u)
    x = R*(cos(Omega).*cos(u)-sin(Omega).*sin(u).*cos(inc));
    y = R*(sin(Omega).*cos(u)+cos(Omega).*sin(u).*cos(inc));
    z = R*(sin(u).*sin(inc));
    positions = [x y z];
end

function velocities = walker_delta_velocities(R,inc,Omega,u,omega)
    dxdu = R*(-cos(Omega).*sin(u)-sin(Omega).*cos(u).*cos(inc));
    dydu = R*(-sin(Omega).*sin(u)+cos(Omega).*cos(u).*cos(inc));
    dzdu = R*(cos(u).*sin(inc));
    velocities = omega*[dxdu dydu dzdu];
end

function positions = ground_user_positions(R,lat,lon)
    positions = [R*cos(lat).*cos(lon), ...
                 R*cos(lat).*sin(lon), ...
                 R*sin(lat)];
end
