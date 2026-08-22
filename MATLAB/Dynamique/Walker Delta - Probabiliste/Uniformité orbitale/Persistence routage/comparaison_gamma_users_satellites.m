clear; clc; close all;

%% ============================================================
% COMPARAISON DE L'ECART ANGULAIRE ENTRE DEUX UTILISATEURS
% ET ENTRE LES SATELLITES QUI LEUR SONT ASSIGNES
%
% Objectif :
% verifier si l'on peut approximer
%
%       gamma_satellites_assignes ~= gamma_utilisateurs
%
% Assignation :
% chaque utilisateur choisit le satellite visible minimisant la
% distance sol-satellite, avec un critere d'elevation minimale.
%% ============================================================

rng(5);

%% Parametres physiques et orbitaux
R_earth = 6371;                  % [km]
h = 550;                         % [km]
R = R_earth + h;                 % [km]

lambda = 4e-7;                   % [satellites/km^2]
N_mean = lambda*4*pi*R^2;

inc_deg = 58;
inc = deg2rad(inc_deg);

elevation_min_deg = 20;
elevation_min = deg2rad(elevation_min_deg);

%% Parametres Monte-Carlo
N_realizations = 100;
N_pairs_per_realization = 200;
N_pairs_total = N_realizations*N_pairs_per_realization;

%% Stockage
gamma_users_deg = NaN(N_pairs_total,1);
gamma_sats_deg = NaN(N_pairs_total,1);
angular_error_deg = NaN(N_pairs_total,1);
absolute_error_deg = NaN(N_pairs_total,1);
same_satellite = false(N_pairs_total,1);
both_covered = false(N_pairs_total,1);

idx_global = 0;

%% Monte-Carlo
for r = 1:N_realizations

    %% Constellation Walker Delta stochastique
    N_sat = poissrnd(N_mean);

    Omega = 2*pi*rand(N_sat,1);
    u = 2*pi*rand(N_sat,1);

    x_sat = R*(cos(Omega).*cos(u) ...
        - sin(Omega).*sin(u).*cos(inc));

    y_sat = R*(sin(Omega).*cos(u) ...
        + cos(Omega).*sin(u).*cos(inc));

    z_sat = R*(sin(u).*sin(inc));

    sat_pos = [x_sat y_sat z_sat];

    %% Paires d'utilisateurs uniformes en surface dans la bande [-i,+i]
    sin_lat_A = -sin(inc) + 2*sin(inc)*rand(N_pairs_per_realization,1);
    sin_lat_B = -sin(inc) + 2*sin(inc)*rand(N_pairs_per_realization,1);

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

        idx_global = idx_global + 1;

        uA = userA_pos(p,:);
        uB = userB_pos(p,:);

        %% Ecart angulaire entre utilisateurs
        cos_gamma_users = dot(uA,uB)/R_earth^2;
        cos_gamma_users = max(-1,min(1,cos_gamma_users));
        gamma_users = acos(cos_gamma_users);

        gamma_users_deg(idx_global) = rad2deg(gamma_users);

        %% Assignation utilisateur A
        rho_A = sat_pos-uA;
        dist_A = sqrt(sum(rho_A.^2,2));
        zenith_A = uA/R_earth;

        sin_el_A = (rho_A*zenith_A.') ./ dist_A;
        el_A = asin(max(-1,min(1,sin_el_A)));

        visible_A = el_A >= elevation_min;

        if any(visible_A)
            ids_A = find(visible_A);
            [~,jA] = min(dist_A(ids_A));
            satA = ids_A(jA);
        else
            satA = NaN;
        end

        %% Assignation utilisateur B
        rho_B = sat_pos-uB;
        dist_B = sqrt(sum(rho_B.^2,2));
        zenith_B = uB/R_earth;

        sin_el_B = (rho_B*zenith_B.') ./ dist_B;
        el_B = asin(max(-1,min(1,sin_el_B)));

        visible_B = el_B >= elevation_min;

        if any(visible_B)
            ids_B = find(visible_B);
            [~,jB] = min(dist_B(ids_B));
            satB = ids_B(jB);
        else
            satB = NaN;
        end

        %% Comparaison
        if isfinite(satA) && isfinite(satB)

            both_covered(idx_global) = true;
            same_satellite(idx_global) = (satA == satB);

            rA = sat_pos(satA,:);
            rB = sat_pos(satB,:);

            cos_gamma_sats = dot(rA,rB)/R^2;
            cos_gamma_sats = max(-1,min(1,cos_gamma_sats));
            gamma_sats = acos(cos_gamma_sats);

            gamma_sats_deg(idx_global) = rad2deg(gamma_sats);

            angular_error_deg(idx_global) = ...
                gamma_sats_deg(idx_global)-gamma_users_deg(idx_global);

            absolute_error_deg(idx_global) = ...
                abs(angular_error_deg(idx_global));
        end
    end

    fprintf('Realisation %3d/%3d terminee\n',r,N_realizations);
end

%% Nettoyage
valid = both_covered;

gamma_user = gamma_users_deg(valid);
gamma_sat = gamma_sats_deg(valid);

error_signed = angular_error_deg(valid);
error_abs = absolute_error_deg(valid);

same_sat_valid = same_satellite(valid);

%% Statistiques
N_valid = nnz(valid);
coverage_pair_probability = N_valid/N_pairs_total;

mean_gamma_user = mean(gamma_user);
mean_gamma_sat = mean(gamma_sat);

mean_error = mean(error_signed);
mean_abs_error = mean(error_abs);
median_abs_error = median(error_abs);
rmse_error = sqrt(mean(error_signed.^2));
max_abs_error = max(error_abs);

correlation_gamma = corr(gamma_user,gamma_sat);
same_satellite_probability = mean(same_sat_valid);

fprintf('\n============================================================\n');
fprintf('COMPARAISON DES ECARTS ANGULAIRES\n');
fprintf('============================================================\n');
fprintf('Nombre total de paires                  : %d\n',N_pairs_total);
fprintf('Paires avec double couverture           : %d\n',N_valid);
fprintf('Probabilite de double couverture        : %.4f\n',coverage_pair_probability);
fprintf('Ecart angulaire moyen utilisateurs      : %.4f deg\n',mean_gamma_user);
fprintf('Ecart angulaire moyen satellites        : %.4f deg\n',mean_gamma_sat);
fprintf('Erreur moyenne gamma_sat-gamma_user     : %.4f deg\n',mean_error);
fprintf('Erreur absolue moyenne                  : %.4f deg\n',mean_abs_error);
fprintf('Erreur absolue mediane                  : %.4f deg\n',median_abs_error);
fprintf('RMSE                                    : %.4f deg\n',rmse_error);
fprintf('Erreur absolue maximale                 : %.4f deg\n',max_abs_error);
fprintf('Correlation gamma_user / gamma_sat      : %.6f\n',correlation_gamma);
fprintf('Probabilite meme satellite assigne      : %.4f\n',same_satellite_probability);
fprintf('============================================================\n');

%% Table de resultats
Results = table( ...
    gamma_user,gamma_sat,error_signed,error_abs,same_sat_valid, ...
    'VariableNames',{ ...
    'GammaUsers_deg', ...
    'GammaAssignedSatellites_deg', ...
    'AngularError_deg', ...
    'AbsoluteAngularError_deg', ...
    'SameAssignedSatellite'});

%% Graphique 1 : comparaison directe
figure;
scatter(gamma_user,gamma_sat,18,'filled');
hold on;

gamma_lim = [0 max([gamma_user;gamma_sat])];
plot(gamma_lim,gamma_lim,'--','LineWidth',1.5);

grid on;
axis equal;
xlim(gamma_lim);
ylim(gamma_lim);

xlabel('Ecart angulaire utilisateurs \gamma_{users} (deg)');
ylabel('Ecart angulaire satellites assignes \gamma_{sat} (deg)');
title(sprintf('Comparaison des ecarts angulaires | r = %.4f, MAE = %.2f deg', ...
    correlation_gamma,mean_abs_error));
legend('Paires','y=x','Location','best');

%% Graphique 2 : erreur signee
figure;
scatter(gamma_user,error_signed,18,'filled');
hold on;
yline(0,'--','LineWidth',1.5);
grid on;

xlabel('Ecart angulaire utilisateurs \gamma_{users} (deg)');
ylabel('\gamma_{sat} - \gamma_{users} (deg)');
title('Erreur sur la separation angulaire apres assignation');

%% Graphique 3 : distribution de l'erreur absolue
figure;
histogram(error_abs,40,'Normalization','probability');
grid on;

xlabel('|\gamma_{sat}-\gamma_{users}| (deg)');
ylabel('Probabilite');
title(sprintf('Erreur absolue | moyenne = %.2f deg, mediane = %.2f deg', ...
    mean_abs_error,median_abs_error));

%% Graphique 4 : erreur absolue moyenne selon gamma_users
gamma_edges = 0:15:180;
gamma_centers = (gamma_edges(1:end-1)+gamma_edges(2:end))/2;

N_bins = numel(gamma_centers);
mean_error_by_gamma = NaN(N_bins,1);
median_error_by_gamma = NaN(N_bins,1);
pair_count_by_gamma = zeros(N_bins,1);

for b = 1:N_bins

    in_bin = gamma_user >= gamma_edges(b) & ...
             gamma_user < gamma_edges(b+1);

    pair_count_by_gamma(b) = nnz(in_bin);

    if any(in_bin)
        mean_error_by_gamma(b) = mean(error_abs(in_bin));
        median_error_by_gamma(b) = median(error_abs(in_bin));
    end
end

figure;
plot(gamma_centers,mean_error_by_gamma,'o-','LineWidth',1.5, ...
    'DisplayName','Erreur absolue moyenne');
hold on;

plot(gamma_centers,median_error_by_gamma,'s--','LineWidth',1.5, ...
    'DisplayName','Erreur absolue mediane');

grid on;
xlabel('Ecart angulaire utilisateurs \gamma_{users} (deg)');
ylabel('Erreur angulaire absolue (deg)');
title('Erreur selon la separation des utilisateurs');
legend('Location','best');