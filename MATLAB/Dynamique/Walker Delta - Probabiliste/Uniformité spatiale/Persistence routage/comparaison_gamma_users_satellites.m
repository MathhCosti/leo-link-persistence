clear; clc; close all;

%% ============================================================
% COMPARAISON DE L'ECART ANGULAIRE
% WALKER DELTA - UNIFORMITE SPATIALE
%
% Objectif :
% verifier si l'on peut approximer
%
%   gamma_satellites_assignes ~= gamma_utilisateurs
%
% dans le modele spatial, et verifier si la qualite de cette
% approximation depend de l'instant t.
%% ============================================================

rng(5);

%% Parametres physiques et orbitaux
R_earth = 6371;
h = 550;
R = R_earth+h;

mu = 398600;
omega_sat = sqrt(mu/R^3);
T_orbit = 2*pi/omega_sat;
T_spatial = T_orbit/2;

lambda = 4e-7;
N_mean = lambda*4*pi*R^2;

inc_deg = 58;
inc = deg2rad(inc_deg);

elevation_min_deg = 20;
elevation_min = deg2rad(elevation_min_deg);

%% Instants testes dans une demi-periode fondamentale
time_samples = linspace(0,T_spatial,5).';
N_times = numel(time_samples);

%% Parametres Monte-Carlo
N_realizations = 50;
N_pairs_per_realization = 100;

N_pairs_total = ...
    N_times*N_realizations*N_pairs_per_realization;

%% Stockage global
time_sample_s = NaN(N_pairs_total,1);
gamma_users_deg = NaN(N_pairs_total,1);
gamma_sats_deg = NaN(N_pairs_total,1);
angular_error_deg = NaN(N_pairs_total,1);
absolute_error_deg = NaN(N_pairs_total,1);
same_satellite = false(N_pairs_total,1);
both_covered = false(N_pairs_total,1);

idx_global = 0;

%% Monte-Carlo
for it = 1:N_times

    t = time_samples(it);

    for r = 1:N_realizations

        %% Constellation Walker Delta spatiale
        N_sat = poissrnd(N_mean);

        Omega = 2*pi*rand(N_sat,1);

        % Uniformite spatiale a t=0.
        u0 = sample_u_spatial(N_sat);

        % Evolution jusqu'a l'instant teste.
        u_t = mod(u0+omega_sat*t,2*pi);

        sat_pos = walker_delta_positions(R,inc,Omega,u_t);

        %% Paires d'utilisateurs uniformes en surface dans [-i,+i]
        sin_lat_A = ...
            -sin(inc)+2*sin(inc)*rand(N_pairs_per_realization,1);

        sin_lat_B = ...
            -sin(inc)+2*sin(inc)*rand(N_pairs_per_realization,1);

        lat_A = asin(sin_lat_A);
        lat_B = asin(sin_lat_B);

        lon_A = 2*pi*rand(N_pairs_per_realization,1)-pi;
        lon_B = 2*pi*rand(N_pairs_per_realization,1)-pi;

        userA_pos = [ ...
            R_earth*cos(lat_A).*cos(lon_A), ...
            R_earth*cos(lat_A).*sin(lon_A), ...
            R_earth*sin(lat_A)];

        userB_pos = [ ...
            R_earth*cos(lat_B).*cos(lon_B), ...
            R_earth*cos(lat_B).*sin(lon_B), ...
            R_earth*sin(lat_B)];

        %% Boucle sur les paires
        for p = 1:N_pairs_per_realization

            idx_global = idx_global+1;
            time_sample_s(idx_global) = t;

            uA = userA_pos(p,:);
            uB = userB_pos(p,:);

            %% Separation utilisateurs
            cos_gamma_users = dot(uA,uB)/R_earth^2;
            cos_gamma_users = max(-1,min(1,cos_gamma_users));

            gamma_users = acos(cos_gamma_users);
            gamma_users_deg(idx_global) = rad2deg(gamma_users);

            %% Assignation utilisateur A
            satA = closest_visible_satellite( ...
                sat_pos,uA,R_earth,elevation_min);

            %% Assignation utilisateur B
            satB = closest_visible_satellite( ...
                sat_pos,uB,R_earth,elevation_min);

            %% Comparaison
            if isfinite(satA) && isfinite(satB)

                both_covered(idx_global) = true;
                same_satellite(idx_global) = (satA==satB);

                rA = sat_pos(satA,:);
                rB = sat_pos(satB,:);

                cos_gamma_sats = dot(rA,rB)/R^2;
                cos_gamma_sats = max(-1,min(1,cos_gamma_sats));

                gamma_sats = acos(cos_gamma_sats);

                gamma_sats_deg(idx_global) = ...
                    rad2deg(gamma_sats);

                angular_error_deg(idx_global) = ...
                    gamma_sats_deg(idx_global) ...
                    - gamma_users_deg(idx_global);

                absolute_error_deg(idx_global) = ...
                    abs(angular_error_deg(idx_global));
            end
        end
    end

    fprintf('Instant %d/%d termine : t = %.1f s\n', ...
        it,N_times,t);
end

%% Nettoyage
valid = both_covered;

time_valid = time_sample_s(valid);
gamma_user = gamma_users_deg(valid);
gamma_sat = gamma_sats_deg(valid);
error_signed = angular_error_deg(valid);
error_abs = absolute_error_deg(valid);
same_sat_valid = same_satellite(valid);

%% Statistiques globales
N_valid = nnz(valid);

coverage_pair_probability = ...
    N_valid/N_pairs_total;

mean_gamma_user = mean(gamma_user);
mean_gamma_sat = mean(gamma_sat);

mean_error = mean(error_signed);
mean_abs_error = mean(error_abs);
median_abs_error = median(error_abs);
rmse_error = sqrt(mean(error_signed.^2));
max_abs_error = max(error_abs);

correlation_gamma = corr(gamma_user,gamma_sat);
same_satellite_probability = mean(same_sat_valid);

%% Statistiques par instant
MeanAbsError_by_time = NaN(N_times,1);
RMSE_by_time = NaN(N_times,1);
Correlation_by_time = NaN(N_times,1);
Coverage_by_time = NaN(N_times,1);
PairCount_by_time = zeros(N_times,1);

for it = 1:N_times

    mask_all = abs(time_sample_s-time_samples(it)) < 1e-9;
    mask = valid & mask_all;

    PairCount_by_time(it) = nnz(mask);

    Coverage_by_time(it) = ...
        nnz(mask)/(N_realizations*N_pairs_per_realization);

    if nnz(mask) >= 2

        gu = gamma_users_deg(mask);
        gs = gamma_sats_deg(mask);
        er = angular_error_deg(mask);

        MeanAbsError_by_time(it) = mean(abs(er));
        RMSE_by_time(it) = sqrt(mean(er.^2));
        Correlation_by_time(it) = corr(gu,gs);
    end
end

TimeStatistics = table( ...
    time_samples,PairCount_by_time,Coverage_by_time, ...
    MeanAbsError_by_time,RMSE_by_time,Correlation_by_time, ...
    'VariableNames',{ ...
    'Time_s','CoveredPairs','DoubleCoverageProbability', ...
    'MeanAbsoluteError_deg','RMSE_deg','Correlation'});

fprintf('\n============================================================\n');
fprintf('COMPARAISON GAMMA - DELTA UNIFORMITE SPATIALE\n');
fprintf('============================================================\n');
fprintf('Nombre total de paires                  : %d\n',N_pairs_total);
fprintf('Paires avec double couverture           : %d\n',N_valid);
fprintf('Probabilite de double couverture        : %.4f\n',coverage_pair_probability);
fprintf('Erreur moyenne gamma_sat-gamma_user     : %.4f deg\n',mean_error);
fprintf('Erreur absolue moyenne                  : %.4f deg\n',mean_abs_error);
fprintf('Erreur absolue mediane                  : %.4f deg\n',median_abs_error);
fprintf('RMSE                                    : %.4f deg\n',rmse_error);
fprintf('Correlation globale                     : %.6f\n',correlation_gamma);
fprintf('Probabilite meme satellite assigne      : %.4f\n',same_satellite_probability);
fprintf('============================================================\n');

disp(TimeStatistics);

%% Table detaillee
Results = table( ...
    time_valid,gamma_user,gamma_sat,error_signed,error_abs,same_sat_valid, ...
    'VariableNames',{ ...
    'Time_s','GammaUsers_deg','GammaAssignedSatellites_deg', ...
    'AngularError_deg','AbsoluteAngularError_deg', ...
    'SameAssignedSatellite'});

%% Figure 1 : comparaison directe
figure;
scatter(gamma_user,gamma_sat,15,time_valid,'filled');
hold on;

gamma_lim = [0 max([gamma_user;gamma_sat])];
plot(gamma_lim,gamma_lim,'k--','LineWidth',1.5);

grid on;
axis equal;
xlim(gamma_lim);
ylim(gamma_lim);
colorbar;

xlabel('\gamma_{users} (deg)');
ylabel('\gamma_{sat} (deg)');
title(sprintf( ...
    'Delta spatial : r = %.4f, MAE = %.2f deg', ...
    correlation_gamma,mean_abs_error));

%% Figure 2 : erreur selon le temps
figure;
plot(time_samples,MeanAbsError_by_time,'o-','LineWidth',1.5);
grid on;
xlabel('Temps (s)');
ylabel('Erreur absolue moyenne (deg)');
title('Erreur de \gamma selon la phase spatiale');

%% Figure 3 : correlation selon le temps
figure;
plot(time_samples,Correlation_by_time,'o-','LineWidth',1.5);
grid on;
xlabel('Temps (s)');
ylabel('Correlation');
title('Correlation \gamma_{users} / \gamma_{sat} selon le temps');

%% Sauvegarde
save('comparaison_gamma_users_satellites_results.mat', ...
    'Results','TimeStatistics', ...
    'time_samples','T_orbit','T_spatial', ...
    'gamma_users_deg','gamma_sats_deg', ...
    'angular_error_deg','absolute_error_deg', ...
    'same_satellite','both_covered', ...
    'mean_abs_error','rmse_error','correlation_gamma', ...
    'coverage_pair_probability');

%% ============================================================
% FONCTIONS LOCALES
%% ============================================================

function u = sample_u_spatial(N)
    s = 2*rand(N,1)-1;
    a = asin(s);
    branch = rand(N,1)<0.5;

    u = a;
    u(~branch) = pi-a(~branch);
    u = mod(u,2*pi);
end

function positions = walker_delta_positions(R,inc,Omega,u)
    x = R*(cos(Omega).*cos(u)-sin(Omega).*sin(u).*cos(inc));
    y = R*(sin(Omega).*cos(u)+cos(Omega).*sin(u).*cos(inc));
    z = R*(sin(u).*sin(inc));
    positions = [x y z];
end

function sat_id = closest_visible_satellite( ...
    sat_pos,user_pos,R_earth,elevation_min)

    rho = sat_pos-user_pos;
    dist = sqrt(sum(rho.^2,2));

    zenith = user_pos/R_earth;
    sin_el = (rho*zenith.') ./ dist;
    el = asin(max(-1,min(1,sin_el)));

    visible = el >= elevation_min;

    if any(visible)
        ids = find(visible);
        [~,j] = min(dist(ids));
        sat_id = ids(j);
    else
        sat_id = NaN;
    end
end
