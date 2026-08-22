clear; clc; close all;

%% ============================================================
% MATRICE EMPIRIQUE T_route(phi_A,phi_B)
% WALKER DELTA STOCHASTIQUE EN UNIFORMITE ORBITALE
%
% Lignes   : phi_A
% Colonnes : phi_B
%
% Les longitudes lambda_A et lambda_B sont imposees et identiques
% pour toutes les realisations. Leur ecart reste constant pendant
% la rotation terrestre.
%
% Pour chaque couple de latitudes, on mesure les episodes de route fixe :
% - assignation au satellite visible le plus proche ;
% - plus court chemin ISL memorise ;
% - fin si assignation A/B change ou si un ISL du chemin disparait.
%% ============================================================

rng(8);

%% Parametres physiques
R_earth = 6371;
h = 550;
R = R_earth+h;
mu = 398600;
omega_sat = sqrt(mu/R^3);
omega_earth = 2*pi/86164;

%% Constellation
lambda = 4e-7;
inc_deg = 58;
inc = deg2rad(inc_deg);

%% Liens
dmax = 1500;
elevation_min_deg = 20;
elevation_min = deg2rad(elevation_min_deg);

%% Latitudes testees
user_lat_deg = (-55:5:55).';
user_lat = deg2rad(user_lat_deg);
N_lat = numel(user_lat_deg);

%% Longitudes imposees
% Meme choix que dans le code theorique.
user_lon_A_deg = 0;
user_lon_B_deg = 90;

user_lon0_A_fixed = deg2rad(user_lon_A_deg);
user_lon0_B_fixed = deg2rad(user_lon_B_deg);

Delta_lon_deg = abs(user_lon_A_deg-user_lon_B_deg);

%% Simulation
N_realizations = 50;
dt = 10;
Tmax = 15000;
time_values = (0:dt:Tmax-dt).';
Nt = numel(time_values);

min_episode_samples = 1;

%% Stockage
route_durations = cell(N_lat,N_lat);
route_hops = cell(N_lat,N_lat);
route_end_cause = cell(N_lat,N_lat);

samples_both_assigned = zeros(N_lat,N_lat);
samples_route_exists = zeros(N_lat,N_lat);
route_starts = zeros(N_lat,N_lat);

total_samples = N_realizations*Nt;

%% ============================================================
% MONTE-CARLO
%% ============================================================

for r = 1:N_realizations

    N_mean = lambda*4*pi*R^2;
    N_sat = poissrnd(N_mean);

    if N_sat < 2
        continue;
    end

    Omega = 2*pi*rand(N_sat,1);
    u0 = 2*pi*rand(N_sat,1);

    % Longitudes imposees pour toutes les latitudes.
    % Tous les utilisateurs A partagent la longitude user_lon_A_deg,
    % et tous les utilisateurs B partagent la longitude user_lon_B_deg.
    %
    % Les deux familles tournent ensuite avec la Terre, donc
    % Delta_lambda reste constant dans le temps.
    user_lon0_A = user_lon0_A_fixed*ones(N_lat,1);
    user_lon0_B = user_lon0_B_fixed*ones(N_lat,1);

    % Etat d'une route pour chaque couple.
    route_active = false(N_lat,N_lat);
    route_start_time = NaN(N_lat,N_lat);
    access_A_initial = zeros(N_lat,N_lat);
    access_B_initial = zeros(N_lat,N_lat);
    route_hops_current = NaN(N_lat,N_lat);
    route_path = cell(N_lat,N_lat);

    for k = 1:Nt

        t = time_values(k);

        %% Positions satellites
        u_t = mod(u0+omega_sat*t,2*pi);

        x_sat = R*(cos(Omega).*cos(u_t) ...
            - sin(Omega).*sin(u_t).*cos(inc));
        y_sat = R*(sin(Omega).*cos(u_t) ...
            + cos(Omega).*sin(u_t).*cos(inc));
        z_sat = R*(sin(u_t).*sin(inc));

        sat_pos = [x_sat y_sat z_sat];

        %% Graphe ISL
        D = squareform(pdist(sat_pos));
        Adj = (D <= dmax) & (D > 0);
        G = graph(Adj);

        %% Positions utilisateurs A
        lon_A = mod(user_lon0_A+omega_earth*t+pi,2*pi)-pi;

        user_pos_A = [ ...
            R_earth*cos(user_lat).*cos(lon_A), ...
            R_earth*cos(user_lat).*sin(lon_A), ...
            R_earth*sin(user_lat)];

        %% Positions utilisateurs B
        lon_B = mod(user_lon0_B+omega_earth*t+pi,2*pi)-pi;

        user_pos_B = [ ...
            R_earth*cos(user_lat).*cos(lon_B), ...
            R_earth*cos(user_lat).*sin(lon_B), ...
            R_earth*sin(user_lat)];

        %% Assignations A
        assignment_A = zeros(N_lat,1);

        for qA = 1:N_lat
            rho = sat_pos-user_pos_A(qA,:);
            dist = sqrt(sum(rho.^2,2));
            zenith = user_pos_A(qA,:)/R_earth;
            sin_el = (rho*zenith.') ./ dist;
            el = asin(max(-1,min(1,sin_el)));
            visible = el >= elevation_min;

            if any(visible)
                ids = find(visible);
                [~,j] = min(dist(ids));
                assignment_A(qA) = ids(j);
            end
        end

        %% Assignations B
        assignment_B = zeros(N_lat,1);

        for qB = 1:N_lat
            rho = sat_pos-user_pos_B(qB,:);
            dist = sqrt(sum(rho.^2,2));
            zenith = user_pos_B(qB,:)/R_earth;
            sin_el = (rho*zenith.') ./ dist;
            el = asin(max(-1,min(1,sin_el)));
            visible = el >= elevation_min;

            if any(visible)
                ids = find(visible);
                [~,j] = min(dist(ids));
                assignment_B(qB) = ids(j);
            end
        end

        %% Tous les couples
        for qA = 1:N_lat
            for qB = 1:N_lat

                a_now = assignment_A(qA);
                b_now = assignment_B(qB);

                both_assigned = (a_now > 0) && (b_now > 0);
                samples_both_assigned(qA,qB) = ...
                    samples_both_assigned(qA,qB) + both_assigned;

                %% Persistance d'une route active
                if route_active(qA,qB)

                    cause = "";

                    if a_now ~= access_A_initial(qA,qB)
                        cause = "assignment_A";

                    elseif b_now ~= access_B_initial(qA,qB)
                        cause = "assignment_B";

                    else
                        path_current = route_path{qA,qB};

                        if numel(path_current) <= 1
                            alive = true;
                        else
                            idx = sub2ind(size(Adj), ...
                                path_current(1:end-1), ...
                                path_current(2:end));
                            alive = all(Adj(idx));
                        end

                        if ~alive
                            cause = "ISL_break";
                        end
                    end

                    if strlength(cause) > 0

                        duration = t-route_start_time(qA,qB);

                        if duration >= min_episode_samples*dt
                            route_durations{qA,qB}(end+1,1) = duration;
                            route_hops{qA,qB}(end+1,1) = ...
                                route_hops_current(qA,qB);
                            route_end_cause{qA,qB}(end+1,1) = cause;
                        end

                        route_active(qA,qB) = false;
                        route_path{qA,qB} = [];
                    end
                end

                %% Nouvelle route
                if ~route_active(qA,qB) && both_assigned

                    if a_now == b_now
                        path_now = a_now;
                    else
                        path_now = shortestpath(G,a_now,b_now, ...
                            'Method','unweighted');
                    end

                    if ~isempty(path_now)
                        route_active(qA,qB) = true;
                        route_start_time(qA,qB) = t;
                        route_path{qA,qB} = path_now;
                        access_A_initial(qA,qB) = a_now;
                        access_B_initial(qA,qB) = b_now;
                        route_hops_current(qA,qB) = ...
                            max(0,numel(path_now)-1);
                        route_starts(qA,qB) = route_starts(qA,qB)+1;
                    end
                end

                if route_active(qA,qB)
                    samples_route_exists(qA,qB) = ...
                        samples_route_exists(qA,qB)+1;
                end
            end
        end
    end

    fprintf('Realisation %2d/%2d terminee, N=%d\n', ...
        r,N_realizations,N_sat);
end

%% ============================================================
% MATRICES DE RESULTATS
%% ============================================================

Troute_emp_matrix = NaN(N_lat,N_lat);
Troute_median_matrix = NaN(N_lat,N_lat);
Troute_std_matrix = NaN(N_lat,N_lat);
EpisodeCount_matrix = zeros(N_lat,N_lat);
MeanHopCount_matrix = NaN(N_lat,N_lat);

for qA = 1:N_lat
    for qB = 1:N_lat

        d = route_durations{qA,qB};
        EpisodeCount_matrix(qA,qB) = numel(d);

        if ~isempty(d)
            Troute_emp_matrix(qA,qB) = mean(d);
            Troute_median_matrix(qA,qB) = median(d);
            Troute_std_matrix(qA,qB) = std(d);
        end

        hh = route_hops{qA,qB};

        if ~isempty(hh)
            MeanHopCount_matrix(qA,qB) = mean(hh);
        end
    end
end

P_both_assigned_matrix = samples_both_assigned/total_samples;
P_route_available_matrix = samples_route_exists/total_samples;

%% Table
lat_labels = matlab.lang.makeValidName( ...
    compose('phiB_%gdeg',user_lat_deg));

Troute_emp_table = array2table( ...
    Troute_emp_matrix, ...
    'VariableNames',lat_labels);

Troute_emp_table = addvars( ...
    Troute_emp_table,user_lat_deg, ...
    'Before',1, ...
    'NewVariableNames','phiA_deg');

disp(Troute_emp_table);

%% ============================================================
% CARTE T_route
%% ============================================================

figure;
imagesc(user_lat_deg,user_lat_deg,Troute_emp_matrix);
set(gca,'YDir','normal');
colorbar;

xlabel('\phi_B (deg)');
ylabel('\phi_A (deg)');
title(sprintf( ...
    'T_{route}^{emp}(\\phi_A,\\phi_B), \\lambda_A=%.0f deg, \\lambda_B=%.0f deg', ...
    user_lon_A_deg,user_lon_B_deg));

for qA = 1:N_lat
    for qB = 1:N_lat
        value = Troute_emp_matrix(qA,qB);

        if isfinite(value)
            text(user_lat_deg(qB),user_lat_deg(qA), ...
                sprintf('%.1f',value), ...
                'HorizontalAlignment','center');
        end
    end
end

%% ============================================================
% CARTE DU NOMBRE MOYEN DE SAUTS
%% ============================================================

figure;
imagesc(user_lat_deg,user_lat_deg,MeanHopCount_matrix);
set(gca,'YDir','normal');
colorbar;

xlabel('\phi_B (deg)');
ylabel('\phi_A (deg)');
title(sprintf( ...
    'Nombre moyen empirique de sauts, \\Delta\\lambda=%.0f deg', ...
    Delta_lon_deg));

%% ============================================================
% CARTE DE DISPONIBILITE
%% ============================================================

figure;
imagesc(user_lat_deg,user_lat_deg,P_route_available_matrix);
set(gca,'YDir','normal');
colorbar;

xlabel('\phi_B (deg)');
ylabel('\phi_A (deg)');
title(sprintf( ...
    'Probabilite empirique de disponibilite, \\Delta\\lambda=%.0f deg', ...
    Delta_lon_deg));

%% ============================================================
% INFORMATIONS GEOMETRIQUES
%% ============================================================

fprintf('\n============================================================\n');
fprintf('LONGITUDES IMPOSEES\n');
fprintf('============================================================\n');
fprintf('Longitude utilisateur A : %.2f deg\n',user_lon_A_deg);
fprintf('Longitude utilisateur B : %.2f deg\n',user_lon_B_deg);
fprintf('Ecart de longitude       : %.2f deg\n',Delta_lon_deg);
fprintf('============================================================\n');

%% ============================================================
% SAUVEGARDE
%% ============================================================

save('T_route_emp_matrix_results.mat', ...
    'Troute_emp_matrix','Troute_emp_table', ...
    'Troute_median_matrix','Troute_std_matrix', ...
    'EpisodeCount_matrix','MeanHopCount_matrix', ...
    'P_both_assigned_matrix','P_route_available_matrix', ...
    'route_durations','route_hops','route_end_cause', ...
    'route_starts','user_lat_deg', ...
    'user_lon_A_deg','user_lon_B_deg','Delta_lon_deg', ...
    'lambda','inc_deg','dmax','elevation_min_deg', ...
    'N_realizations','dt','Tmax');

fprintf('\nResultats sauvegardes dans :\n');
fprintf('  T_route_emp_matrix_results.mat\n');
