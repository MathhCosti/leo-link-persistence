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
% Le diagnostic ISL par lien est effectue independamment de cette
% classification : tous les ISL disparus sont comptes directement
% a chaque pas de temps, meme si une assignation change simultanement.
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
dt = 2;                          % s
Tmax = 15000;                    % s
time_values = (0:dt:Tmax-dt).';
Nt = numel(time_values);

% Si true, les utilisateurs sont tires uniformement en surface dans
% la bande de latitude [-inc,+inc]. Si false, ils sont uniformes
% sur toute la Terre.
restrict_users_to_orbital_band = true;

% Nombre minimal d'echantillons constituant un episode valide.
min_episode_samples = 1;

%% Fichier facultatif pour comparaison avec la theorie par quadrature
script_dir = fileparts(mfilename('fullpath'));

quadrature_file = fullfile( ...
    script_dir,'T_route_th_quadrature_delta_results.mat');

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

% Diagnostic ISL direct, independant de la cause de fin de route.
%
% A chaque transition temporelle d'une route active :
%   - chaque lien du chemin contribue dt secondes-lien d'exposition ;
%   - tous les liens disparus sont comptes, meme si une assignation
%     change simultanement.
direct_link_exposure = 0;          % [s.lien]
direct_ISL_break_count = 0;        % nombre total de liens disparus
direct_route_steps_with_ISL_break = 0;

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

            %% ------------------------------------------------
            % Diagnostic ISL DIRECT
            %
            % On examine les liens du chemin memorise AVANT de choisir
            % la cause de fin de route. Ainsi, une rupture ISL est
            % comptabilisee meme si une assignation change pendant
            % le meme intervalle temporel.
            %% ------------------------------------------------
            H_current = max(0,numel(route_path)-1);

            if H_current > 0

                % Tous les H liens ont ete exposes pendant l'intervalle
                % [t-dt,t].
                direct_link_exposure = ...
                    direct_link_exposure + H_current*dt;

                idx_links = sub2ind(size(A), ...
                    route_path(1:end-1), ...
                    route_path(2:end));

                links_alive_now = A(idx_links);
                n_broken_links_now = nnz(~links_alive_now);

                if n_broken_links_now > 0
                    direct_ISL_break_count = ...
                        direct_ISL_break_count + n_broken_links_now;

                    direct_route_steps_with_ISL_break = ...
                        direct_route_steps_with_ISL_break + 1;
                end
            end

            %% Cause de fin de route
            % Cette classification est conservee pour les statistiques
            % de route, mais n'intervient plus dans l'estimation directe
            % du taux de rupture par lien.
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
% DIAGNOSTIC DES TAUX EMPIRIQUES DE DISPARITION
%% ============================================================

% Temps total pendant lequel les episodes complets ont ete observes.
%
% Dans un modele de risques concurrents exponentiels :
%
%   beta_A   = N_A   / somme_k T_k
%   beta_B   = N_B   / somme_k T_k
%   beta_ISL = N_ISL / somme_k T_k
%
% et
%
%   beta_total = beta_A + beta_B + beta_ISL
%              = N_episodes / somme_k T_k
%              = 1 / mean(T_route).

total_route_exposure = sum(route_durations);

N_end_A = nnz(route_end_cause == "assignment_A");
N_end_B = nnz(route_end_cause == "assignment_B");
N_end_ISL = nnz(route_end_cause == "ISL_break");

beta_A_emp = N_end_A/total_route_exposure;
beta_B_emp = N_end_B/total_route_exposure;
beta_ISL_emp = N_end_ISL/total_route_exposure;

beta_GSL_emp = beta_A_emp + beta_B_emp;
beta_total_emp_from_causes = beta_GSL_emp + beta_ISL_emp;

Troute_from_empirical_hazards = 1/beta_total_emp_from_causes;

% Controle numerique :
% doit etre identique a mean(route_durations), a l'arrondi pres.
hazard_consistency_error = ...
    Troute_from_empirical_hazards - Troute_emp_mean;

%% ============================================================
% DIAGNOSTIC EMPIRIQUE PAR LIEN ISL
%% ============================================================

% ------------------------------------------------------------
% A. Ancien estimateur base sur la cause de fin de route
% ------------------------------------------------------------
%
% Cet estimateur est conserve uniquement comme comparaison :
%
%   beta_link_cause =
%       N(fin classee ISL) / somme_k H_k*T_k
%
% Il peut sous-estimer le taux ISL lorsque plusieurs evenements
% surviennent pendant le meme pas de temps.

link_exposure_per_episode = ...
    route_hops_initial .* route_durations;

total_link_exposure_cause = ...
    sum(link_exposure_per_episode);

if total_link_exposure_cause > 0
    beta_link_emp_cause = ...
        N_end_ISL/total_link_exposure_cause;
else
    beta_link_emp_cause = NaN;
end

% ------------------------------------------------------------
% B. Estimateur DIRECT retenu
% ------------------------------------------------------------
%
% A chaque intervalle dt pendant lequel une route est active :
%
%   exposition += H*dt
%
% puis chaque ISL du chemin qui n'existe plus au nouvel instant
% est compte comme une rupture, independamment des changements GSL.
%
%   beta_link_emp_direct =
%       N_ruptures_ISL_direct / exposition_directe

if direct_link_exposure > 0
    beta_link_emp_direct = ...
        direct_ISL_break_count/direct_link_exposure;
else
    beta_link_emp_direct = NaN;
end

if isfinite(beta_link_emp_direct) && beta_link_emp_direct > 0
    Tlink_emp_direct = 1/beta_link_emp_direct;
else
    Tlink_emp_direct = NaN;
end

if isfinite(beta_link_emp_cause) && beta_link_emp_cause > 0
    Tlink_emp_cause = 1/beta_link_emp_cause;
else
    Tlink_emp_cause = NaN;
end

% La valeur directe devient la valeur empirique principale par lien.
beta_link_emp = beta_link_emp_direct;
Tlink_emp = Tlink_emp_direct;
total_link_exposure = direct_link_exposure;

%% ============================================================
% COMPARAISON AVEC LA THEORIE PAR QUADRATURE
%% ============================================================

Troute_theory_quad = NaN;

beta_GSL_theory_quad = NaN;
beta_ISL_theory_quad = NaN;
beta_total_theory_quad = NaN;

% Les deux utilisateurs sont statistiquement identiques dans la moyenne
% globale par quadrature. On partage donc la contribution GSL en deux
% uniquement pour faciliter le diagnostic A/B.
beta_A_theory_quad = NaN;
beta_B_theory_quad = NaN;

H_mean_theory_quad = NaN;
beta_link_theory_quad = NaN;
Tlink_theory_quad = NaN;

if isfile(quadrature_file)

    qd = load(quadrature_file);

    if isfield(qd,'Troute_mean')
        Troute_theory_quad = double(qd.Troute_mean);
    end

    if isfield(qd,'MeanBetaRouteGSL')
        beta_GSL_theory_quad = double(qd.MeanBetaRouteGSL);
        beta_A_theory_quad = beta_GSL_theory_quad/2;
        beta_B_theory_quad = beta_GSL_theory_quad/2;
    end

    if isfield(qd,'MeanBetaRouteISL')
        beta_ISL_theory_quad = double(qd.MeanBetaRouteISL);
    end

    if isfield(qd,'H_mean_from_distribution')
        H_mean_theory_quad = double(qd.H_mean_from_distribution);
    end

    if isfinite(beta_ISL_theory_quad) && ...
            isfinite(H_mean_theory_quad) && H_mean_theory_quad > 0
        beta_link_theory_quad = ...
            beta_ISL_theory_quad/H_mean_theory_quad;

        if beta_link_theory_quad > 0
            Tlink_theory_quad = 1/beta_link_theory_quad;
        end
    end

    if isfinite(beta_GSL_theory_quad) && isfinite(beta_ISL_theory_quad)
        beta_total_theory_quad = ...
            beta_GSL_theory_quad + beta_ISL_theory_quad;
    end
end

%% Ecarts theorie / empirique
error_beta_A_percent = NaN;
error_beta_B_percent = NaN;
error_beta_GSL_percent = NaN;
error_beta_ISL_percent = NaN;
error_beta_link_percent = NaN;
error_beta_total_percent = NaN;

if isfinite(beta_A_theory_quad) && beta_A_emp > 0
    error_beta_A_percent = ...
        100*(beta_A_theory_quad-beta_A_emp)/beta_A_emp;
end

if isfinite(beta_B_theory_quad) && beta_B_emp > 0
    error_beta_B_percent = ...
        100*(beta_B_theory_quad-beta_B_emp)/beta_B_emp;
end

if isfinite(beta_GSL_theory_quad) && beta_GSL_emp > 0
    error_beta_GSL_percent = ...
        100*(beta_GSL_theory_quad-beta_GSL_emp)/beta_GSL_emp;
end

if isfinite(beta_ISL_theory_quad) && beta_ISL_emp > 0
    error_beta_ISL_percent = ...
        100*(beta_ISL_theory_quad-beta_ISL_emp)/beta_ISL_emp;
end

if isfinite(beta_link_theory_quad) && ...
        isfinite(beta_link_emp) && beta_link_emp > 0
    error_beta_link_percent = ...
        100*(beta_link_theory_quad-beta_link_emp)/beta_link_emp;
end

if isfinite(beta_total_theory_quad) && beta_total_emp_from_causes > 0
    error_beta_total_percent = ...
        100*(beta_total_theory_quad-beta_total_emp_from_causes) ...
        /beta_total_emp_from_causes;
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
    total_route_exposure, ...
    beta_A_emp,beta_B_emp,beta_GSL_emp,beta_ISL_emp, ...
    beta_total_emp_from_causes,Troute_from_empirical_hazards, ...
    Troute_theory_quad,beta_GSL_theory_quad, ...
    beta_ISL_theory_quad,beta_total_theory_quad, ...
    error_beta_GSL_percent,error_beta_ISL_percent, ...
    error_beta_total_percent, ...
    'VariableNames',{ ...
    'NumberOfCompleteEpisodes','MeanEmpiricalRouteLifetime_s', ...
    'MedianEmpiricalRouteLifetime_s','StdEmpiricalRouteLifetime_s', ...
    'CI95Low_s','CI95High_s','MeanInitialHopCount', ...
    'ProbabilityBothUsersAssigned','ProbabilityRouteAvailable', ...
    'TotalRouteExposure_s', ...
    'BetaAssignA_Emp_per_s','BetaAssignB_Emp_per_s', ...
    'BetaGSL_Emp_per_s','BetaISL_Emp_per_s', ...
    'BetaTotal_Emp_per_s','RouteLifetimeFromEmpiricalHazards_s', ...
    'RouteLifetimeTheoryQuadrature_s', ...
    'BetaGSL_TheoryQuadrature_per_s', ...
    'BetaISL_TheoryQuadrature_per_s', ...
    'BetaTotal_TheoryQuadrature_per_s', ...
    'BetaGSL_Error_percent','BetaISL_Error_percent', ...
    'BetaTotal_Error_percent'});

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

fprintf('\n--- TAUX EMPIRIQUES PAR CAUSE ---\n');
fprintf('Temps total d''exposition           : %.3f s\n', ...
    total_route_exposure);
fprintf('N fin assignation A                : %d\n',N_end_A);
fprintf('N fin assignation B                : %d\n',N_end_B);
fprintf('N fin rupture ISL                  : %d\n',N_end_ISL);

fprintf('\nbeta_A empirique                    : %.8e s^-1\n', ...
    beta_A_emp);
fprintf('beta_B empirique                    : %.8e s^-1\n', ...
    beta_B_emp);
fprintf('beta_GSL empirique                  : %.8e s^-1\n', ...
    beta_GSL_emp);
fprintf('beta_ISL empirique                  : %.8e s^-1\n', ...
    beta_ISL_emp);
fprintf('beta_total empirique                : %.8e s^-1\n', ...
    beta_total_emp_from_causes);
fprintf('1/beta_total empirique              : %.3f s\n', ...
    Troute_from_empirical_hazards);
fprintf('Erreur controle vs moyenne directe  : %.3e s\n', ...
    hazard_consistency_error);

fprintf('\n--- DIAGNOSTIC EMPIRIQUE PAR LIEN ---\n');

fprintf('Methode directe (retenue) :\n');
fprintf('Exposition directe des liens        : %.3f s.lien\n', ...
    direct_link_exposure);
fprintf('Nombre direct de ruptures ISL       : %d\n', ...
    direct_ISL_break_count);
fprintf('Pas avec au moins une rupture ISL   : %d\n', ...
    direct_route_steps_with_ISL_break);
fprintf('beta_link empirique direct          : %.8e s^-1\n', ...
    beta_link_emp_direct);
fprintf('Duree equivalente lien directe      : %.3f s\n', ...
    Tlink_emp_direct);

fprintf('\nAncienne methode par cause :\n');
fprintf('Exposition episodes complets        : %.3f s.lien\n', ...
    total_link_exposure_cause);
fprintf('Ruptures classees ISL               : %d\n', ...
    N_end_ISL);
fprintf('beta_link empirique par cause       : %.8e s^-1\n', ...
    beta_link_emp_cause);
fprintf('Duree equivalente lien par cause    : %.3f s\n', ...
    Tlink_emp_cause);

if isfinite(Troute_theory_quad)

    fprintf('\n--- THEORIE PAR QUADRATURE ---\n');
    fprintf('Troute theorique quadrature         : %.3f s\n', ...
        Troute_theory_quad);
    fprintf('beta_GSL theorie                    : %.8e s^-1\n', ...
        beta_GSL_theory_quad);
    fprintf('beta_ISL theorie                    : %.8e s^-1\n', ...
        beta_ISL_theory_quad);
    fprintf('beta_total theorie                  : %.8e s^-1\n', ...
        beta_total_theory_quad);
    fprintf('E[H] theorie quadrature             : %.6f\n', ...
        H_mean_theory_quad);
    fprintf('beta_link theorie                   : %.8e s^-1\n', ...
        beta_link_theory_quad);
    fprintf('Duree equivalente lien theorie      : %.3f s\n', ...
        Tlink_theory_quad);

    fprintf('\n--- ECARTS THEORIE / EMPIRIQUE ---\n');
    fprintf('Ecart beta_GSL                      : %+7.2f %%\n', ...
        error_beta_GSL_percent);
    fprintf('Ecart beta_ISL                      : %+7.2f %%\n', ...
        error_beta_ISL_percent);
    fprintf('Ecart beta_link                     : %+7.2f %%\n', ...
        error_beta_link_percent);
    fprintf('Ecart beta_total                    : %+7.2f %%\n', ...
        error_beta_total_percent);
    fprintf('Ecart Troute                        : %+7.2f %%\n', ...
        100*(Troute_theory_quad-Troute_emp_mean)/Troute_emp_mean);

else
    fprintf('\nComparaison quadrature non disponible :\n');
    fprintf('  %s introuvable.\n',quadrature_file);
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

if isfinite(Troute_theory_quad)
    xline(Troute_theory_quad,'-.','LineWidth',1.6, ...
        'DisplayName',sprintf('Theorie = %.2f s',Troute_theory_quad));
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

if isfinite(Troute_theory_quad)
    t_curve = linspace(0,max(sorted_T),500);
    plot(t_curve,exp(-t_curve/Troute_theory_quad),'--', ...
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
% DIAGNOSTIC GRAPHIQUE DES TAUX
%% ============================================================

figure;

if isfinite(beta_GSL_theory_quad) && isfinite(beta_ISL_theory_quad)

    beta_compare = [ ...
        beta_GSL_emp, beta_ISL_emp, beta_total_emp_from_causes; ...
        beta_GSL_theory_quad, beta_ISL_theory_quad, beta_total_theory_quad];

    bar(beta_compare.');

    set(gca,'XTickLabel',{'GSL','ISL','Total'});
    ylabel('Taux de disparition (s^{-1})');
    title('Diagnostic des contributions a la disparition de la route');
    legend({'Empirique','Theorie quadrature'},'Location','best');
    grid on;

else

    bar([beta_A_emp,beta_B_emp,beta_ISL_emp]);
    set(gca,'XTickLabel',{'Assignation A','Assignation B','ISL'});
    ylabel('Taux empirique (s^{-1})');
    title('Taux empiriques de disparition par cause');
    grid on;

end

%% ============================================================
% DIAGNOSTIC SPECIFIQUE DU TAUX DE RUPTURE PAR LIEN
%% ============================================================

if isfinite(beta_link_emp_direct)
    figure;

    if isfinite(beta_link_theory_quad)

        bar([ ...
            beta_link_emp_direct, ...
            beta_link_emp_cause, ...
            beta_link_theory_quad]);

        set(gca,'XTickLabel',{ ...
            'Empirique direct', ...
            'Empirique par cause', ...
            'Theorie quadrature'});

    else

        bar([beta_link_emp_direct,beta_link_emp_cause]);

        set(gca,'XTickLabel',{ ...
            'Empirique direct', ...
            'Empirique par cause'});
    end

    ylabel('Taux de rupture par lien (s^{-1})');
    title('Rupture d''un lien appartenant a une route');
    grid on;
end

%% ============================================================
% SAUVEGARDE
%% ============================================================

% writetable(EpisodeResults,'route_lifetime_empirical_episodes.csv');
% writetable(Summary,'route_lifetime_empirical_summary.csv');
% writetable(CauseStatistics,'route_lifetime_end_causes.csv');

save('T_route_emp_results.mat', ...
    'EpisodeResults','Summary','CauseStatistics', ...
    'route_durations','route_hops_initial', ...
    'route_gamma_ground_deg','route_end_cause', ...
    'Troute_emp_mean','Troute_emp_median','Troute_emp_std', ...
    'ci95_low','ci95_high','H_route_emp_mean', ...
    'P_both_assigned','P_route_available', ...
    'total_route_exposure', ...
    'beta_A_emp','beta_B_emp','beta_GSL_emp','beta_ISL_emp', ...
    'beta_total_emp_from_causes','Troute_from_empirical_hazards', ...
    'link_exposure_per_episode', ...
    'total_link_exposure','total_link_exposure_cause', ...
    'direct_link_exposure','direct_ISL_break_count', ...
    'direct_route_steps_with_ISL_break', ...
    'beta_link_emp','Tlink_emp', ...
    'beta_link_emp_direct','Tlink_emp_direct', ...
    'beta_link_emp_cause','Tlink_emp_cause', ...
    'Troute_theory_quad','beta_GSL_theory_quad', ...
    'beta_ISL_theory_quad','beta_total_theory_quad', ...
    'H_mean_theory_quad','beta_link_theory_quad','Tlink_theory_quad', ...
    'error_beta_GSL_percent','error_beta_ISL_percent', ...
    'error_beta_link_percent','error_beta_total_percent', ...
    'lambda','inc_deg','dmax','elevation_min_deg', ...
    'N_realizations','dt','Tmax');

fprintf('\nResultats sauvegardes dans :\n');
fprintf('  route_lifetime_empirical_episodes.csv\n');
fprintf('  route_lifetime_empirical_summary.csv\n');
fprintf('  route_lifetime_end_causes.csv\n');
fprintf('  T_route_emp_results.mat\n');

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
