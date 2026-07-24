clear; clc; close all;

%% ============================================================
% DUREE DE VIE EMPIRIQUE D'UNE ROUTE SOL--SATELLITES--SOL
% WALKER DELTA STOCHASTIQUE EN UNIFORMITE ORBITALE
%
% Pour chaque realisation :
%   1) deux utilisateurs au sol sont tires au hasard ;
%   2) a chaque instant, chacun est assigne au satellite visible
%      le plus proche ;
%   3) lorsqu'un chemin ISL existe entre les deux satellites d'acces,
%      le plus court chemin courant est memorise ;
%   4) cette route fixe disparait des que :
%        - l'assignation de l'utilisateur A change ;
%        - l'assignation de l'utilisateur B change ;
%        - au moins un ISL du chemin memorise disparait.
%
% Un nouvel episode peut ensuite commencer des qu'une nouvelle route
% complete existe.
%
% La moyenne empirique des durees est comparee, si les fichiers existent,
% a la fermeture theorique :
%
%   T_route_th = 1 / (1/T_assign,A + 1/T_assign,B
%                     + H_mean*beta_break_ISL).
%% ============================================================

rng(8);

%% Parametres physiques
R_earth = 6371;                  % km
h = 550;                         % km
R = R_earth+h;                   % km
mu = 398600;                     % km^3/s^2
omega_sat = sqrt(mu/R^3);        % rad/s
omega_earth = 2*pi/86164;        % rad/s

%% Constellation Delta aleatoire
lambda = 4e-7;                   % satellites/km^2
inc_deg = 58;
inc = deg2rad(inc_deg);

%% Liens
dmax = 1500;                     % km, portee ISL
elevation_min_deg = 20;
elevation_min = deg2rad(elevation_min_deg);

%% Simulation
N_realizations = 50;
dt = 5;                          % s
Tmax = 15000;                    % s
time_values = (0:dt:Tmax-dt).';
Nt = numel(time_values);

% Si true, les utilisateurs sont tires uniformement en surface dans
% la bande de latitude [-inc,+inc]. Si false, ils sont uniformes
% sur toute la Terre.
restrict_users_to_orbital_band = true;

% Nombre minimal d'echantillons constituant un episode valide.
min_episode_samples = 1;

%% Fichiers facultatifs pour comparaison theorique
script_dir = fileparts(mfilename('fullpath'));

assignment_file = fullfile( ...
    script_dir,'verification_duree_assignation_corrigee.mat');
shortest_path_file = fullfile( ...
    script_dir,'distribution_shortest_path_delta.mat');
pbreak_file = fullfile( ...
    script_dir,'pbreak_theorique_walker_delta_results.mat');

%% Stockage des episodes
route_durations = [];
route_hops_initial = [];
route_gamma_ground_deg = [];
route_end_cause = strings(0,1);
route_realization = [];
route_user_lat_A_deg = [];
route_user_lat_B_deg = [];

% Statistiques temporelles
total_samples = N_realizations*Nt;
samples_both_assigned = 0;
samples_route_exists = 0;
route_starts = 0;

%% ============================================================
% BOUCLE MONTE-CARLO
%% ============================================================
for r = 1:N_realizations

    %% Constellation
    N_mean = lambda*4*pi*R^2;
    N_sat = poissrnd(N_mean);

    if N_sat < 2
        warning('Realisation %d ignoree : moins de deux satellites.',r);
        continue;
    end

    Omega = 2*pi*rand(N_sat,1);
    u0 = 2*pi*rand(N_sat,1);

    %% Deux utilisateurs fixes sur la Terre tournante
    if restrict_users_to_orbital_band
        % Uniformite en surface dans la bande : sin(latitude) uniforme.
        sin_lat = (2*rand(2,1)-1)*sin(inc);
        user_lat = asin(sin_lat);
    else
        % Uniformite en surface sur toute la Terre.
        user_lat = asin(2*rand(2,1)-1);
    end

    user_lon0 = 2*pi*rand(2,1)-pi;

    % Separation angulaire entre les deux points au sol, constante.
    cos_gamma_ground = ...
        sin(user_lat(1))*sin(user_lat(2)) + ...
        cos(user_lat(1))*cos(user_lat(2))* ...
        cos(user_lon0(1)-user_lon0(2));
    cos_gamma_ground = max(-1,min(1,cos_gamma_ground));
    gamma_ground_deg = rad2deg(acos(cos_gamma_ground));

    %% Etat de la route courante
    route_active = false;
    route_start_time = NaN;
    route_path = [];
    access_A_initial = 0;
    access_B_initial = 0;
    route_hops = NaN;

    for k = 1:Nt
        t = time_values(k);

        %% Positions satellitaires
        u_t = mod(u0+omega_sat*t,2*pi);
        sat_pos = walker_delta_positions(R,inc,Omega,u_t);

        %% Graphe ISL courant
        D = squareform(pdist(sat_pos));
        A = (D <= dmax) & (D > 0);
        G = graph(A);

        %% Positions des utilisateurs
        lon_t = mod(user_lon0+omega_earth*t+pi,2*pi)-pi;
        user_pos = ground_user_positions(R_earth,user_lat,lon_t);

        %% Assignation au satellite visible le plus proche
        current_assignment = zeros(2,1);

        for q = 1:2
            rho = sat_pos-user_pos(q,:);
            dist = sqrt(sum(rho.^2,2));
            zenith = user_pos(q,:)/R_earth;

            sin_el = (rho*zenith.') ./ dist;
            el = asin(max(-1,min(1,sin_el)));

            visible = el >= elevation_min;

            if any(visible)
                ids = find(visible);
                [~,j] = min(dist(ids));
                current_assignment(q) = ids(j);
            end
        end

        both_assigned = all(current_assignment > 0);
        samples_both_assigned = samples_both_assigned+both_assigned;

        %% Si une route est active, tester sa persistance
        if route_active
            cause = "";

            if current_assignment(1) ~= access_A_initial
                cause = "assignment_A";
            elseif current_assignment(2) ~= access_B_initial
                cause = "assignment_B";
            elseif ~fixed_path_is_alive(A,route_path)
                cause = "ISL_break";
            end

            if strlength(cause) > 0
                duration = t-route_start_time;

                if duration >= min_episode_samples*dt
                    route_durations(end+1,1) = duration; %#ok<SAGROW>
                    route_hops_initial(end+1,1) = route_hops; %#ok<SAGROW>
                    route_gamma_ground_deg(end+1,1) = ...
                        gamma_ground_deg; %#ok<SAGROW>
                    route_end_cause(end+1,1) = cause; %#ok<SAGROW>
                    route_realization(end+1,1) = r; %#ok<SAGROW>
                    route_user_lat_A_deg(end+1,1) = ...
                        rad2deg(user_lat(1)); %#ok<SAGROW>
                    route_user_lat_B_deg(end+1,1) = ...
                        rad2deg(user_lat(2)); %#ok<SAGROW>
                end

                route_active = false;
                route_path = [];
            end
        end

        %% Si aucune route active, tenter d'en etablir une nouvelle
        if ~route_active && both_assigned
            a = current_assignment(1);
            b = current_assignment(2);

            if a == b
                % Les deux utilisateurs utilisent le meme satellite :
                % route sans ISL.
                path_now = a;
            else
                path_now = shortestpath(G,a,b,'Method','unweighted');
            end

            if ~isempty(path_now)
                samples_route_exists = samples_route_exists+1;

                route_active = true;
                route_start_time = t;
                route_path = path_now;
                access_A_initial = a;
                access_B_initial = b;
                route_hops = max(0,numel(path_now)-1);
                route_starts = route_starts+1;
            end
        elseif route_active
            samples_route_exists = samples_route_exists+1;
        end
    end

    % Episode encore actif a Tmax : censure a droite, donc exclu de la
    % moyenne empirique non censuree.
    fprintf(['Realisation %2d/%2d : N=%d, latitudes=(%.1f,%.1f) deg, ' ...
             'episodes cumules=%d\n'], ...
             r,N_realizations,N_sat,rad2deg(user_lat(1)), ...
             rad2deg(user_lat(2)),numel(route_durations));
end

%% ============================================================
% STATISTIQUES EMPIRIQUES
%% ============================================================
if isempty(route_durations)
    error(['Aucun episode complet de route n''a ete observe. ' ...
           'Augmenter Tmax/N_realizations ou verifier dmax et lambda.']);
end

Troute_emp_mean = mean(route_durations);
Troute_emp_median = median(route_durations);
Troute_emp_std = std(route_durations);
Troute_emp_sem = Troute_emp_std/sqrt(numel(route_durations));

ci95_low = Troute_emp_mean-1.96*Troute_emp_sem;
ci95_high = Troute_emp_mean+1.96*Troute_emp_sem;

H_route_emp_mean = mean(route_hops_initial);
P_both_assigned = samples_both_assigned/total_samples;
P_route_available = samples_route_exists/total_samples;

cause_categories = ["assignment_A","assignment_B","ISL_break"];
cause_counts = zeros(numel(cause_categories),1);
for j = 1:numel(cause_categories)
    cause_counts(j) = nnz(route_end_cause == cause_categories(j));
end
cause_prob = cause_counts/sum(cause_counts);

%% ============================================================
% COMPARAISON THEORIQUE FACULTATIVE
%% ============================================================
Troute_theory = NaN;
H_mean_theory_used = NaN;
Tassign_A_theory = NaN;
Tassign_B_theory = NaN;
beta_break_ISL = NaN;

if isfile(assignment_file) && isfile(shortest_path_file) && ...
        isfile(pbreak_file)

    ad = load(assignment_file);
    hd = load(shortest_path_file);
    bd = load(pbreak_file);

    % H moyen theorique s'il est sauvegarde, sinon H moyen empirique.
    if isfield(hd,'mean_H_theory')
        H_mean_theory_used = double(hd.mean_H_theory);
    elseif isfield(hd,'mean_H')
        H_mean_theory_used = double(hd.mean_H);
        warning('mean_H_theory absent : mean_H empirique utilise.');
    end

    % Durees d'assignation corrigees aux latitudes disponibles les plus
    % proches de la moyenne des latitudes tirees.
    if isfield(ad,'MeanAssign_corrected') && ...
            isfield(ad,'lambda_values') && ...
            isfield(ad,'user_lat_deg')

        [~,il] = min(abs(double(ad.lambda_values(:))-lambda));

        mean_lat_A = mean(route_user_lat_A_deg);
        mean_lat_B = mean(route_user_lat_B_deg);

        [~,qA] = min(abs(double(ad.user_lat_deg(:))-mean_lat_A));
        [~,qB] = min(abs(double(ad.user_lat_deg(:))-mean_lat_B));

        Tassign_A_theory = double(ad.MeanAssign_corrected(il,qA));
        Tassign_B_theory = double(ad.MeanAssign_corrected(il,qB));
    end

    if isfield(bd,'q_break_link') && isfield(bd,'Delta_t')
        q_break_link = double(bd.q_break_link);
        Delta_t_break = double(bd.Delta_t);
        beta_break_ISL = -log1p(-q_break_link)/Delta_t_break;
    end

    if all(isfinite([H_mean_theory_used,Tassign_A_theory, ...
            Tassign_B_theory,beta_break_ISL]))
        Troute_theory = 1/( ...
            1/Tassign_A_theory + ...
            1/Tassign_B_theory + ...
            H_mean_theory_used*beta_break_ISL);
    end
end

%% ============================================================
% TABLES
%% ============================================================
EpisodeResults = table( ...
    route_realization,route_durations,route_hops_initial, ...
    route_gamma_ground_deg,route_user_lat_A_deg, ...
    route_user_lat_B_deg,route_end_cause, ...
    'VariableNames',{ ...
    'Realization','RouteLifetime_s','InitialHopCount', ...
    'GroundAngularSeparation_deg','UserLatitudeA_deg', ...
    'UserLatitudeB_deg','EndCause'});

Summary = table( ...
    numel(route_durations),Troute_emp_mean,Troute_emp_median, ...
    Troute_emp_std,ci95_low,ci95_high,H_route_emp_mean, ...
    P_both_assigned,P_route_available, ...
    Troute_theory,H_mean_theory_used,Tassign_A_theory, ...
    Tassign_B_theory,beta_break_ISL, ...
    'VariableNames',{ ...
    'NumberOfCompleteEpisodes','MeanEmpiricalRouteLifetime_s', ...
    'MedianEmpiricalRouteLifetime_s','StdEmpiricalRouteLifetime_s', ...
    'CI95Low_s','CI95High_s','MeanInitialHopCount', ...
    'ProbabilityBothUsersAssigned','ProbabilityRouteAvailable', ...
    'TheoreticalRouteLifetime_s','TheoreticalMeanHopCountUsed', ...
    'TheoreticalAssignA_s','TheoreticalAssignB_s', ...
    'TheoreticalISLBreakRate_per_s'});

CauseStatistics = table( ...
    cause_categories.',cause_counts,cause_prob, ...
    'VariableNames',{'EndCause','Count','Probability'});

disp(Summary);
disp(CauseStatistics);

fprintf('\n============================================================\n');
fprintf('DUREE EMPIRIQUE DE ROUTE\n');
fprintf('============================================================\n');
fprintf('Episodes complets                 : %d\n',numel(route_durations));
fprintf('Troute empirique moyen            : %.3f s\n',Troute_emp_mean);
fprintf('Troute empirique median           : %.3f s\n',Troute_emp_median);
fprintf('IC 95 %% de la moyenne             : [%.3f, %.3f] s\n', ...
    ci95_low,ci95_high);
fprintf('H initial moyen des routes        : %.3f\n',H_route_emp_mean);
fprintf('P(deux utilisateurs assignes)     : %.4f\n',P_both_assigned);
fprintf('P(route disponible)               : %.4f\n',P_route_available);

if isfinite(Troute_theory)
    fprintf('Troute theorique compare          : %.3f s\n',Troute_theory);
    fprintf('Erreur relative theorie/empirique : %.2f %%\n', ...
        100*(Troute_theory-Troute_emp_mean)/Troute_emp_mean);
else
    fprintf('Comparaison theorique non disponible.\n');
end
fprintf('============================================================\n');

%% ============================================================
% GRAPHIQUES
%% ============================================================

% Distribution des durees.
figure;
histogram(route_durations,'Normalization','probability');
hold on;
xline(Troute_emp_mean,'--','LineWidth',1.6, ...
    'DisplayName',sprintf('Moyenne empirique = %.2f s', ...
    Troute_emp_mean));

if isfinite(Troute_theory)
    xline(Troute_theory,'-.','LineWidth',1.6, ...
        'DisplayName',sprintf('Theorie = %.2f s',Troute_theory));
end

grid on;
xlabel('Duree de vie de la route (s)');
ylabel('Probabilite');
title('Distribution empirique de T_{route}');
legend('Location','best');
hold off;

% Fonction de survie empirique.
sorted_T = sort(route_durations);
survival_emp = 1-(1:numel(sorted_T)).'/numel(sorted_T);

figure;
stairs(sorted_T,survival_emp,'LineWidth',1.7, ...
    'DisplayName','Survie empirique');
hold on;

if isfinite(Troute_theory)
    t_curve = linspace(0,max(sorted_T),500);
    plot(t_curve,exp(-t_curve/Troute_theory),'--', ...
        'LineWidth',1.7, ...
        'DisplayName','Exponentielle theorique');
end

grid on;
xlabel('Temps t (s)');
ylabel('P(T_{route}>t)');
title('Fonction de survie de la route');
legend('Location','best');
hold off;

% Causes de fin.
figure;
bar(categorical(cause_categories),cause_prob);
grid on;
ylabel('Proportion des episodes');
title('Causes empiriques de disparition de la route');

% Duree selon le nombre initial de sauts.
unique_h = unique(route_hops_initial);
mean_T_by_h = NaN(size(unique_h));
count_by_h = zeros(size(unique_h));

for j = 1:numel(unique_h)
    mask = route_hops_initial == unique_h(j);
    mean_T_by_h(j) = mean(route_durations(mask));
    count_by_h(j) = nnz(mask);
end

figure;
plot(unique_h,mean_T_by_h,'o-','LineWidth',1.6);
grid on;
xlabel('Nombre initial de sauts ISL');
ylabel('Duree moyenne de route (s)');
title('Duree empirique selon le nombre de sauts');

%% ============================================================
% SAUVEGARDE
%% ============================================================
writetable(EpisodeResults,'route_lifetime_empirical_episodes.csv');
writetable(Summary,'route_lifetime_empirical_summary.csv');
writetable(CauseStatistics,'route_lifetime_end_causes.csv');

save('route_lifetime_empirical_results.mat', ...
    'EpisodeResults','Summary','CauseStatistics', ...
    'route_durations','route_hops_initial', ...
    'route_gamma_ground_deg','route_end_cause', ...
    'Troute_emp_mean','Troute_emp_median','Troute_emp_std', ...
    'ci95_low','ci95_high','H_route_emp_mean', ...
    'P_both_assigned','P_route_available','Troute_theory', ...
    'lambda','inc_deg','dmax','elevation_min_deg', ...
    'N_realizations','dt','Tmax');

fprintf('\nResultats sauvegardes dans :\n');
fprintf('  route_lifetime_empirical_episodes.csv\n');
fprintf('  route_lifetime_empirical_summary.csv\n');
fprintf('  route_lifetime_end_causes.csv\n');
fprintf('  route_lifetime_empirical_results.mat\n');

%% ============================================================
% FONCTIONS LOCALES
%% ============================================================
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

function alive = fixed_path_is_alive(A,path)
    % Une route composee d'un seul satellite ne contient aucun ISL.
    if numel(path) <= 1
        alive = true;
        return;
    end

    idx = sub2ind(size(A),path(1:end-1),path(2:end));
    alive = all(A(idx));
end
