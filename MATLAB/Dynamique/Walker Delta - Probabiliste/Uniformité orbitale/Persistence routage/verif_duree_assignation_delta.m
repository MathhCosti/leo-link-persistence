clear; clc; close all;

%% ============================================================
% VERIFICATION EMPIRIQUE / HEURISTIQUE CORRIGEE :
% DUREE D'ASSIGNATION AU SATELLITE VISIBLE LE PLUS PROCHE
%
% A chaque instant, chaque utilisateur choisit le satellite admissible
% minimisant la distance sol-satellite.
%
% Deux heuristiques sont comparees :
%
% (1) Heuristique initiale de regime dense :
%       1/T_assign ~= 1/T_vis + 1/T_nearest
%
% (2) Heuristique corrigee de regime peu dense :
%       T_assign ~= p_cov / (nu_reconnexion + nu_HO_corr)
%
% avec
%       p_cov = 1-exp(-m)
%       nu_reconnexion ~= m exp(-m)/T_vis
%       nu_HO_corr ~= P(Nvis>=2 | Nvis>=1) nu_HO_base
%       nu_HO_base = (4/pi) v_rel sqrt(rho_eff)
%
% et
%       P(Nvis>=2 | Nvis>=1)
%       = [1-exp(-m)(1+m)]/[1-exp(-m)].
%
% m est ici le nombre moyen THEORIQUE de satellites visibles, calcule
% avec la probabilite exacte de visibilite du modele Delta orbital.
%% ============================================================

rng(3);

%% Parametres physiques
R_earth = 6371;                    % km
h = 550;                           % km
R = R_earth+h;                     % km
mu = 398600;                       % km^3/s^2
omega_sat = sqrt(mu/R^3);          % rad/s
omega_earth = 2*pi/86164;          % rad/s

%% Modele orbital Delta aleatoire
inc_deg = 58;
inc = deg2rad(inc_deg);

%% Critere de visibilite GSL
elevation_min_deg = 20;
elevation_min = deg2rad(elevation_min_deg);
psi_max = acos((R_earth/R)*cos(elevation_min))-elevation_min;

% Aire de la calotte visible sur la Terre [km^2].
A_visible = 2*pi*R_earth^2*(1-cos(psi_max));

%% Densites testees
lambda_values = [1e-7 2e-7 4e-7 8e-7 1.2e-6].';
N_lambda = numel(lambda_values);

%% Utilisateurs
user_lat_deg = [0 30 45].';
user_lat = deg2rad(user_lat_deg);
N_users = numel(user_lat);
user_lon0 = 2*pi*rand(N_users,1)-pi;

%% Simulation Monte-Carlo
N_realizations = 30;
dt = 5;
Tmax = 12000;
time_values = (0:dt:Tmax-dt).';
Nt = numel(time_values);

%% Stockage
MeanAssign_emp = NaN(N_lambda,N_users);
MeanVisible_emp = NaN(N_lambda,N_users);
Outage_emp = NaN(N_lambda,N_users);

MeanVisible_theory = NaN(N_lambda,N_users);
Outage_theory = NaN(N_lambda,N_users);
Tvis_theory = NaN(N_lambda,N_users);
Tnearest_theory = NaN(N_lambda,N_users);
MeanAssign_old = NaN(N_lambda,N_users);
MeanAssign_corrected = NaN(N_lambda,N_users);

Pmulti_cond = NaN(N_lambda,N_users);
NuReconnection = NaN(N_lambda,N_users);
NuHO_base = NaN(N_lambda,N_users);
NuHO_corrected = NaN(N_lambda,N_users);

%% Probabilite exacte de visibilite selon la latitude
pvis_exact = zeros(N_users,1);
for q = 1:N_users
    pvis_exact(q) = visibility_probability_delta_exact( ...
        user_lat(q),inc,psi_max);
end

%% ============================================================
% BOUCLE PRINCIPALE
%% ============================================================
for il = 1:N_lambda
    lambda = lambda_values(il);
    N_mean = lambda*4*pi*R^2;

    durations = cell(N_users,1);
    total_visible = zeros(N_users,1);
    total_outage = zeros(N_users,1);
    total_samples = N_realizations*Nt;

    for r = 1:N_realizations
        N_sat = poissrnd(N_mean);
        Omega = 2*pi*rand(N_sat,1);
        u0 = 2*pi*rand(N_sat,1);

        previous_assignment = zeros(N_users,1);
        episode_start = zeros(N_users,1);

        for k = 1:Nt
            t = time_values(k);

            %% Satellites dans le repere inertiel
            u_t = mod(u0+omega_sat*t,2*pi);
            sat_pos = walker_delta_positions(R,inc,Omega,u_t);

            %% Utilisateurs fixes sur la Terre tournante
            lon_t = mod(user_lon0+omega_earth*t+pi,2*pi)-pi;
            user_pos = ground_user_positions(R_earth,user_lat,lon_t);

            current_assignment = zeros(N_users,1);

            for q = 1:N_users
                rho = sat_pos-user_pos(q,:);
                dist = sqrt(sum(rho.^2,2));
                zenith = user_pos(q,:)/R_earth;

                sin_el = (rho*zenith.') ./ dist;
                el = asin(max(-1,min(1,sin_el)));

                admissible = el >= elevation_min;
                nvis = nnz(admissible);

                total_visible(q) = total_visible(q)+nvis;
                total_outage(q) = total_outage(q)+(nvis == 0);

                if nvis > 0
                    ids = find(admissible);
                    [~,j] = min(dist(ids));
                    current_assignment(q) = ids(j);
                end
            end

            %% Extraction des episodes continus d'assignation
            if k == 1
                previous_assignment = current_assignment;
                episode_start(:) = t;
            else
                changed = current_assignment ~= previous_assignment;

                for q = find(changed).'
                    if previous_assignment(q) > 0
                        durations{q}(end+1,1) = ...
                            t-episode_start(q); %#ok<SAGROW>
                    end
                    episode_start(q) = t;
                end

                previous_assignment = current_assignment;
            end
        end

        % Les episodes encore ouverts a Tmax sont exclus afin de ne pas
        % introduire une duree tronquee dans la moyenne empirique.
    end

    %% Statistiques et modeles analytiques
    for q = 1:N_users
        if ~isempty(durations{q})
            MeanAssign_emp(il,q) = mean(durations{q});
        end

        MeanVisible_emp(il,q) = total_visible(q)/total_samples;
        Outage_emp(il,q) = total_outage(q)/total_samples;

        phi = user_lat(q);

        %% Nombre visible et outage theoriques
        m = N_mean*pvis_exact(q);
        MeanVisible_theory(il,q) = m;
        Outage_theory(il,q) = exp(-m);
        p_cov = 1-exp(-m);

        %% Duree analytique d'un passage visible
        nu_rel = sqrt(omega_sat^2+(omega_earth*cos(phi))^2);
        v_rel = R_earth*nu_rel;
        Tvis = pi*R_earth*psi_max/(2*v_rel);
        Tvis_theory(il,q) = Tvis;

        %% Densite effective moyenne dans la calotte visible
        % Cette definition est plus robuste pres des bords de la bande
        % orbitale que la densite locale ponctuelle divergente.
        rho_eff = m/A_visible;

        %% Heuristique initiale Poisson-Voronoi
        if rho_eff > 0
            Tnearest = pi/(4*v_rel*sqrt(rho_eff));
            nu_HO0 = 1/Tnearest;
        else
            Tnearest = Inf;
            nu_HO0 = 0;
        end

        Tnearest_theory(il,q) = Tnearest;
        NuHO_base(il,q) = nu_HO0;
        MeanAssign_old(il,q) = 1/(1/Tvis+nu_HO0);

        %% Correction pour le regime peu dense
        if p_cov > eps
            p_multi = (1-exp(-m)*(1+m))/p_cov;
        else
            p_multi = 0;
        end

        % Taux de reprise de couverture apres outage.
        nu_reconnect = m*exp(-m)/Tvis;

        % Un handover ne peut avoir lieu que si au moins deux satellites
        % sont visibles conditionnellement a la couverture.
        nu_HO_corr = p_multi*nu_HO0;

        Pmulti_cond(il,q) = p_multi;
        NuReconnection(il,q) = nu_reconnect;
        NuHO_corrected(il,q) = nu_HO_corr;

        nu_start = nu_reconnect+nu_HO_corr;

        if nu_start > 0
            MeanAssign_corrected(il,q) = p_cov/nu_start;
        else
            MeanAssign_corrected(il,q) = NaN;
        end
    end
end

%% ============================================================
% TABLE DE RESULTATS
%% ============================================================
rows = N_lambda*N_users;
Lambda = zeros(rows,1);
Latitude_deg = zeros(rows,1);
Nmean = zeros(rows,1);
MeanVisibleEmp = zeros(rows,1);
MeanVisibleTheory = zeros(rows,1);
OutageEmp = zeros(rows,1);
OutageTheory = zeros(rows,1);
AssignEmp_s = zeros(rows,1);
TvisTheory_s = zeros(rows,1);
TnearestTheory_s = zeros(rows,1);
PmultiGivenCovered = zeros(rows,1);
NuReconnect_per_s = zeros(rows,1);
NuHOBase_per_s = zeros(rows,1);
NuHOCorrected_per_s = zeros(rows,1);
AssignOld_s = zeros(rows,1);
AssignCorrected_s = zeros(rows,1);
ErrorOld_pct = zeros(rows,1);
ErrorCorrected_pct = zeros(rows,1);

idx = 0;
for il = 1:N_lambda
    for q = 1:N_users
        idx = idx+1;

        Lambda(idx) = lambda_values(il);
        Latitude_deg(idx) = user_lat_deg(q);
        Nmean(idx) = lambda_values(il)*4*pi*R^2;

        MeanVisibleEmp(idx) = MeanVisible_emp(il,q);
        MeanVisibleTheory(idx) = MeanVisible_theory(il,q);
        OutageEmp(idx) = Outage_emp(il,q);
        OutageTheory(idx) = Outage_theory(il,q);

        AssignEmp_s(idx) = MeanAssign_emp(il,q);
        TvisTheory_s(idx) = Tvis_theory(il,q);
        TnearestTheory_s(idx) = Tnearest_theory(il,q);
        PmultiGivenCovered(idx) = Pmulti_cond(il,q);
        NuReconnect_per_s(idx) = NuReconnection(il,q);
        NuHOBase_per_s(idx) = NuHO_base(il,q);
        NuHOCorrected_per_s(idx) = NuHO_corrected(il,q);
        AssignOld_s(idx) = MeanAssign_old(il,q);
        AssignCorrected_s(idx) = MeanAssign_corrected(il,q);

        ErrorOld_pct(idx) = 100*abs(AssignOld_s(idx)-AssignEmp_s(idx)) ...
            /max(AssignEmp_s(idx),eps);
        ErrorCorrected_pct(idx) = ...
            100*abs(AssignCorrected_s(idx)-AssignEmp_s(idx)) ...
            /max(AssignEmp_s(idx),eps);
    end
end

Results = table( ...
    Lambda,Latitude_deg,Nmean, ...
    MeanVisibleEmp,MeanVisibleTheory,OutageEmp,OutageTheory, ...
    AssignEmp_s,TvisTheory_s,TnearestTheory_s, ...
    PmultiGivenCovered,NuReconnect_per_s, ...
    NuHOBase_per_s,NuHOCorrected_per_s, ...
    AssignOld_s,AssignCorrected_s, ...
    ErrorOld_pct,ErrorCorrected_pct);

disp(Results);

%% ============================================================
% GRAPHIQUES
%% ============================================================
figure;
hold on;
for q = 1:N_users
    plot(lambda_values,MeanAssign_emp(:,q),'o-', ...
        'LineWidth',1.5, ...
        'DisplayName',sprintf('Empirique, lat %.0f deg', ...
        user_lat_deg(q)));

    plot(lambda_values,MeanAssign_corrected(:,q),'--', ...
        'LineWidth',1.5, ...
        'DisplayName',sprintf('Heuristique corrigee, lat %.0f deg', ...
        user_lat_deg(q)));
end
grid on;
xlabel('Densite \lambda (satellites/km^2)');
ylabel('Duree moyenne d''assignation (s)');
title('Assignation au plus proche : empirique et heuristique corrigee');
legend('Location','best');

figure;
hold on;
for q = 1:N_users
    plot(lambda_values,MeanAssign_emp(:,q),'o-', ...
        'LineWidth',1.5, ...
        'DisplayName',sprintf('Empirique, lat %.0f deg', ...
        user_lat_deg(q)));

    plot(lambda_values,MeanAssign_old(:,q),':', ...
        'LineWidth',1.5, ...
        'DisplayName',sprintf('Ancienne heuristique, lat %.0f deg', ...
        user_lat_deg(q)));

    plot(lambda_values,MeanAssign_corrected(:,q),'--', ...
        'LineWidth',1.5, ...
        'DisplayName',sprintf('Heuristique corrigee, lat %.0f deg', ...
        user_lat_deg(q)));
end
grid on;
xlabel('Densite \lambda (satellites/km^2)');
ylabel('Duree moyenne d''assignation (s)');
title('Comparaison des deux heuristiques');
legend('Location','best');

figure;
hold on;
for q = 1:N_users
    plot(lambda_values,MeanVisible_emp(:,q),'o-', ...
        'LineWidth',1.5, ...
        'DisplayName',sprintf('Visible empirique, lat %.0f deg', ...
        user_lat_deg(q)));

    plot(lambda_values,MeanVisible_theory(:,q),'--', ...
        'LineWidth',1.5, ...
        'DisplayName',sprintf('Visible theorique, lat %.0f deg', ...
        user_lat_deg(q)));
end
grid on;
xlabel('Densite \lambda (satellites/km^2)');
ylabel('Nombre moyen de satellites visibles');
title('Nombre moyen visible utilise dans la correction');
legend('Location','best');

figure;
hold on;
for q = 1:N_users
    plot(lambda_values,Pmulti_cond(:,q),'o-', ...
        'LineWidth',1.5, ...
        'DisplayName',sprintf('lat %.0f deg',user_lat_deg(q)));
end
grid on;
xlabel('Densite \lambda (satellites/km^2)');
ylabel('P(N_{vis}\geq2\mid N_{vis}\geq1)');
title('Facteur correctif associe a la multiplicite de couverture');
legend('Location','best');

writetable(Results,'verification_duree_assignation_corrigee.csv');
save('verification_duree_assignation_corrigee.mat');

fprintf('\nResultats sauvegardes dans :\n');
fprintf('  verification_duree_assignation_corrigee.csv\n');
fprintf('  verification_duree_assignation_corrigee.mat\n');

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

function pvis = visibility_probability_delta_exact(phi_user,inc,psi_max)
    % Probabilite qu'un satellite du modele Delta aleatoire soit visible
    % depuis un utilisateur de latitude phi_user.
    %
    % u est uniforme sur [0,2*pi] et le RAAN uniforme implique une
    % longitude satellitaire uniforme conditionnellement a u.

    integrand = @(u) longitude_visible_fraction( ...
        u,phi_user,inc,psi_max);

    pvis = integral(integrand,0,2*pi, ...
        'ArrayValued',true, ...
        'RelTol',1e-9, ...
        'AbsTol',1e-11)/(2*pi);
end

function fraction = longitude_visible_fraction(u,phi_user,inc,psi_max)
    phi_sat = asin(sin(inc).*sin(u));

    denominator = cos(phi_user).*cos(phi_sat);
    numerator = cos(psi_max)-sin(phi_user).*sin(phi_sat);

    q = numerator./denominator;

    delta_theta = zeros(size(q));
    delta_theta(q <= -1) = pi;

    middle = (q > -1)&(q < 1);
    delta_theta(middle) = acos(q(middle));

    % Pour une longitude uniforme sur [0,2*pi], la fraction visible vaut
    % 2*delta_theta/(2*pi) = delta_theta/pi.
    fraction = delta_theta/pi;
end
