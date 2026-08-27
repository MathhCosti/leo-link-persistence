clear; clc; close all;

%% ============================================================
% DUREE D'ASSIGNATION AU SATELLITE VISIBLE LE PLUS PROCHE
% WALKER DELTA - UNIFORMITE SPATIALE
%
% Theorie locale quasi stationnaire :
%
%   T_assign(t,phi)
%       ~= p_cov(t,phi) /
%          [nu_reconnexion(t,phi) + nu_HO_corr(t,phi)]
%
% avec
%
%   m(t,phi) = E[N_vis(t,phi)]
%   p_cov(t,phi) = 1-exp[-m(t,phi)]
%
%   nu_reconnexion(t,phi)
%       ~= m(t,phi) exp[-m(t,phi)] / T_vis(phi)
%
%   nu_HO_corr(t,phi)
%       ~= chi_HO(t,phi) * (4/pi) v_rel(phi)
%           sqrt(rho_sat(t,phi))
%
%   chi_HO(t,phi)
%       = P(N_vis>=2 | N_vis>=1)
%
% Dans le modele spatial :
%
%   f_U(u,0) = |cos(u)|/4,
%
% puis
%
%   f_U(u,t) = |cos(u-omega_sat*t)|/4.
%
% Les resultats theoriques et empiriques sont donc calcules en fonction
% de la phase temporelle dans une periode orbitale.
%% ============================================================

rng(3);

%% ============================================================
% 1. PARAMETRES PHYSIQUES
%% ============================================================

R_earth = 6371;                    % [km]
h = 550;                           % [km]
R = R_earth+h;                     % [km]

mu = 398600;                       % [km^3/s^2]
omega_sat = sqrt(mu/R^3);          % [rad/s]
omega_earth = 2*pi/86164;          % [rad/s]

T_orbit = 2*pi/omega_sat;          % [s]

%% ============================================================
% 2. MODELE WALKER DELTA
%% ============================================================

inc_deg = 58;
inc = deg2rad(inc_deg);

%% ============================================================
% 3. CRITERE DE VISIBILITE
%% ============================================================

elevation_min_deg = 20;
elevation_min = deg2rad(elevation_min_deg);

psi_max = ...
    acos((R_earth/R)*cos(elevation_min))-elevation_min;

% Aire de la calotte visible projetee sur la Terre.
A_visible = ...
    2*pi*R_earth^2*(1-cos(psi_max));

%% ============================================================
% 4. DENSITES TESTEES
%% ============================================================

lambda_values = ...
    [1e-7 4e-7].';

N_lambda = numel(lambda_values);

%% ============================================================
% 5. UTILISATEURS
%% ============================================================

user_lat_deg = ...
    [-55 -50 -45 -40 -35 -30 -25 -20 -15 -10 -5 ...
      0 5 10 15 20 25 30 35 40 45 50 55].';

user_lat = deg2rad(user_lat_deg);

N_users = numel(user_lat);

% Une longitude initiale fixe par utilisateur, commune a toutes
% les realisations et a toutes les densites.
user_lon0 = 2*pi*rand(N_users,1)-pi;

%% ============================================================
% 6. DISCRETISATION TEMPORELLE
%% ============================================================

N_realizations = 200;

dt = 5;                           % [s]
Tmax = 12000;                     % [s]

time_values = (0:dt:Tmax-dt).';
Nt = numel(time_values);

% La theorie spatiale est periodique de periode T_orbit.
% On regroupe les donnees empiriques selon la phase temporelle
% modulo T_orbit.
N_time_bins = 24;

time_edges = linspace(0,T_orbit,N_time_bins+1).';
time_phase = ...
    0.5*(time_edges(1:end-1)+time_edges(2:end));

%% ============================================================
% 7. STOCKAGE
%
% Dimensions :
%   (lambda, temps, latitude)
%% ============================================================

MeanAssign_emp = ...
    NaN(N_lambda,N_time_bins,N_users);

MeanVisible_emp = ...
    NaN(N_lambda,N_time_bins,N_users);

Outage_emp = ...
    NaN(N_lambda,N_time_bins,N_users);

MeanVisible_theory = ...
    NaN(N_lambda,N_time_bins,N_users);

Outage_theory = ...
    NaN(N_lambda,N_time_bins,N_users);

Tvis_theory = ...
    NaN(N_lambda,N_time_bins,N_users);

MeanAssign_theory = ...
    NaN(N_lambda,N_time_bins,N_users);

Pmulti_cond = ...
    NaN(N_lambda,N_time_bins,N_users);

NuReconnection = ...
    NaN(N_lambda,N_time_bins,N_users);

NuHO_base = ...
    NaN(N_lambda,N_time_bins,N_users);

NuHO_corrected = ...
    NaN(N_lambda,N_time_bins,N_users);

RhoSat_theory = ...
    NaN(N_lambda,N_time_bins,N_users);

%% ============================================================
% 8. PRECALCUL DE LA GEOMETRIE DE VISIBILITE
%% ============================================================

Nu = 12000;

u_grid = linspace(0,2*pi,Nu);

phi_sat_grid = ...
    asin(sin(inc).*sin(u_grid));

longitude_fraction = ...
    zeros(N_users,Nu);

for q = 1:N_users

    phi_u = user_lat(q);

    denominator = ...
        cos(phi_u).*cos(phi_sat_grid);

    numerator = ...
        cos(psi_max) ...
        - sin(phi_u).*sin(phi_sat_grid);

    qq = numerator./denominator;

    fraction = zeros(size(u_grid));

    fraction(qq <= -1) = 1;
    fraction(qq >= 1) = 0;

    middle = (qq > -1) & (qq < 1);

    fraction(middle) = ...
        acos(qq(middle))/pi;

    singular = ...
        abs(denominator) < 1e-14;

    if any(singular)

        angular_distance = ...
            abs(phi_sat_grid(singular)-phi_u);

        fraction(singular) = ...
            double(angular_distance <= psi_max);
    end

    longitude_fraction(q,:) = fraction;
end

%% ============================================================
% 9. BOUCLE PRINCIPALE SUR LES DENSITES
%% ============================================================

for il = 1:N_lambda

    lambda = lambda_values(il);

    N_mean = ...
        lambda*4*pi*R^2;

    %% --------------------------------------------------------
    % Stockage empirique par classe temporelle
    %% --------------------------------------------------------

    durations_by_bin = ...
        cell(N_time_bins,N_users);

    visible_sum_bin = ...
        zeros(N_time_bins,N_users);

    outage_sum_bin = ...
        zeros(N_time_bins,N_users);

    sample_count_bin = ...
        zeros(N_time_bins,N_users);

    %% --------------------------------------------------------
    % Monte-Carlo
    %% --------------------------------------------------------

    for r = 1:N_realizations

        N_sat = poissrnd(N_mean);

        if N_sat < 1
            continue;
        end

        Omega = ...
            2*pi*rand(N_sat,1);

        % UNIFORMITE SPATIALE A t=0 :
        %
        %   f_U(u,0)=|cos(u)|/4.
        u0 = ...
            sample_u_spatial(N_sat);

        previous_assignment = ...
            zeros(N_users,1);

        episode_start = ...
            NaN(N_users,1);

        for k = 1:Nt

            t = time_values(k);

            %% Satellites
            u_t = ...
                mod(u0+omega_sat*t,2*pi);

            sat_pos = ...
                walker_delta_positions( ...
                R,inc,Omega,u_t);

            %% Utilisateurs
            lon_t = ...
                mod(user_lon0 ...
                + omega_earth*t ...
                + pi,2*pi)-pi;

            user_pos = ...
                ground_user_positions( ...
                R_earth,user_lat,lon_t);

            current_assignment = ...
                zeros(N_users,1);

            nvis_current = ...
                zeros(N_users,1);

            %% Assignation de chaque utilisateur
            for q = 1:N_users

                rho = ...
                    sat_pos-user_pos(q,:);

                dist = ...
                    sqrt(sum(rho.^2,2));

                zenith = ...
                    user_pos(q,:)/R_earth;

                sin_el = ...
                    (rho*zenith.') ./ dist;

                el = ...
                    asin(max(-1,min(1,sin_el)));

                admissible = ...
                    el >= elevation_min;

                nvis = nnz(admissible);

                nvis_current(q) = nvis;

                if nvis > 0

                    ids = find(admissible);

                    [~,j] = ...
                        min(dist(ids));

                    current_assignment(q) = ...
                        ids(j);
                end
            end

            %% Classe temporelle modulo la periode orbitale
            phase_t = ...
                mod(t,T_orbit);

            b_now = ...
                discretize(phase_t,time_edges);

            % Cas numerique eventuel a la borne droite.
            if isnan(b_now)
                b_now = N_time_bins;
            end

            for q = 1:N_users

                visible_sum_bin(b_now,q) = ...
                    visible_sum_bin(b_now,q) ...
                    + nvis_current(q);

                outage_sum_bin(b_now,q) = ...
                    outage_sum_bin(b_now,q) ...
                    + (nvis_current(q)==0);

                sample_count_bin(b_now,q) = ...
                    sample_count_bin(b_now,q)+1;
            end

            %% Episodes continus d'assignation
            if k == 1

                previous_assignment = ...
                    current_assignment;

                for q = 1:N_users

                    if current_assignment(q) > 0
                        episode_start(q) = t;
                    else
                        episode_start(q) = NaN;
                    end
                end

            else

                changed = ...
                    current_assignment ...
                    ~= previous_assignment;

                for q = find(changed).'

                    % Fermeture de l'ancien episode.
                    if previous_assignment(q) > 0 ...
                            && isfinite(episode_start(q))

                        duration = ...
                            t-episode_start(q);

                        % On classe l'episode suivant sa phase de debut.
                        phase_start = ...
                            mod(episode_start(q),T_orbit);

                        b_start = ...
                            discretize(phase_start,time_edges);

                        if isnan(b_start)
                            b_start = N_time_bins;
                        end

                        durations_by_bin{b_start,q}(end+1,1) = ...
                            duration; %#ok<SAGROW>
                    end

                    % Debut eventuel d'un nouvel episode.
                    if current_assignment(q) > 0
                        episode_start(q) = t;
                    else
                        episode_start(q) = NaN;
                    end
                end

                previous_assignment = ...
                    current_assignment;
            end
        end

        % Les episodes encore ouverts a Tmax sont exclus
        % afin d'eviter une censure a droite.
    end

    %% --------------------------------------------------------
    % Statistiques empiriques
    %% --------------------------------------------------------

    for b = 1:N_time_bins
        for q = 1:N_users

            if sample_count_bin(b,q) > 0

                MeanVisible_emp(il,b,q) = ...
                    visible_sum_bin(b,q) ...
                    / sample_count_bin(b,q);

                Outage_emp(il,b,q) = ...
                    outage_sum_bin(b,q) ...
                    / sample_count_bin(b,q);
            end

            d = durations_by_bin{b,q};

            if ~isempty(d)

                MeanAssign_emp(il,b,q) = ...
                    mean(d);
            end
        end
    end

    %% --------------------------------------------------------
    % Theorie spatiale
    %% --------------------------------------------------------

    for b = 1:N_time_bins

        t = time_phase(b);

        % Loi instantanee de phase :
        %
        %   f_U(u,t)=|cos(u-omega*t)|/4.
        fU_t = ...
            abs(cos(u_grid-omega_sat*t))/4;

        for q = 1:N_users

            phi = user_lat(q);

            %% Nombre moyen visible exact
            p_vis = ...
                trapz( ...
                u_grid, ...
                fU_t .* longitude_fraction(q,:));

            m = ...
                N_mean*p_vis;

            MeanVisible_theory(il,b,q) = m;

            Outage_theory(il,b,q) = ...
                exp(-m);

            p_cov = ...
                1-exp(-m);

            %% Duree theorique de visibilite
            nu_rel = ...
                sqrt( ...
                omega_sat^2 ...
                + (omega_earth*cos(phi))^2);

            v_rel = ...
                R_earth*nu_rel;

            Tvis = ...
                pi*R_earth*psi_max ...
                /(2*v_rel);

            Tvis_theory(il,b,q) = ...
                Tvis;

            %% Densite surfacique locale instantanee
            fphi_local = ...
                latitude_pdf_spatial( ...
                phi,t,inc,omega_sat);

            rho_sat = ...
                N_mean*fphi_local ...
                /(2*pi*R_earth^2 ...
                * max(cos(phi),eps));

            RhoSat_theory(il,b,q) = ...
                rho_sat;

            %% Taux de handover Poisson-Voronoi de base
            if rho_sat > 0

                nu_HO0 = ...
                    (4/pi)*v_rel*sqrt(rho_sat);

            else

                nu_HO0 = 0;
            end

            NuHO_base(il,b,q) = ...
                nu_HO0;

            %% Facteur faible densite
            if p_cov > eps

                p_multi = ...
                    (1-exp(-m)*(1+m)) ...
                    / p_cov;

            else

                p_multi = 0;
            end

            %% Taux de reconnexion
            if Tvis > 0

                nu_reconnect = ...
                    m*exp(-m)/Tvis;

            else

                nu_reconnect = 0;
            end

            %% Taux de handover corrige
            nu_HO_corr = ...
                p_multi*nu_HO0;

            Pmulti_cond(il,b,q) = ...
                p_multi;

            NuReconnection(il,b,q) = ...
                nu_reconnect;

            NuHO_corrected(il,b,q) = ...
                nu_HO_corr;

            %% Duree moyenne d'assignation
            nu_start = ...
                nu_reconnect ...
                + nu_HO_corr;

            if nu_start > 0

                MeanAssign_theory(il,b,q) = ...
                    p_cov/nu_start;

            else

                MeanAssign_theory(il,b,q) = ...
                    NaN;
            end
        end
    end

    fprintf( ...
        'Densite %d/%d terminee : lambda = %.2e\n', ...
        il,N_lambda,lambda);
end

%% ============================================================
% 10. TABLE DE RESULTATS
%% ============================================================

rows = ...
    N_lambda*N_time_bins*N_users;

Lambda = zeros(rows,1);
Time_s = zeros(rows,1);
Latitude_deg = zeros(rows,1);

MeanVisibleEmp = NaN(rows,1);
MeanVisibleTheory = NaN(rows,1);

OutageEmp = NaN(rows,1);
OutageTheory = NaN(rows,1);

AssignEmp_s = NaN(rows,1);
AssignTheory_s = NaN(rows,1);

TvisTheory_s = NaN(rows,1);

RhoSat_per_km2 = NaN(rows,1);

PmultiGivenCovered = NaN(rows,1);

NuReconnect_per_s = NaN(rows,1);
NuHOBase_per_s = NaN(rows,1);
NuHOCorrected_per_s = NaN(rows,1);

ErrorTheory_pct = NaN(rows,1);

idx = 0;

for il = 1:N_lambda
    for b = 1:N_time_bins
        for q = 1:N_users

            idx = idx+1;

            Lambda(idx) = ...
                lambda_values(il);

            Time_s(idx) = ...
                time_phase(b);

            Latitude_deg(idx) = ...
                user_lat_deg(q);

            MeanVisibleEmp(idx) = ...
                MeanVisible_emp(il,b,q);

            MeanVisibleTheory(idx) = ...
                MeanVisible_theory(il,b,q);

            OutageEmp(idx) = ...
                Outage_emp(il,b,q);

            OutageTheory(idx) = ...
                Outage_theory(il,b,q);

            AssignEmp_s(idx) = ...
                MeanAssign_emp(il,b,q);

            AssignTheory_s(idx) = ...
                MeanAssign_theory(il,b,q);

            TvisTheory_s(idx) = ...
                Tvis_theory(il,b,q);

            RhoSat_per_km2(idx) = ...
                RhoSat_theory(il,b,q);

            PmultiGivenCovered(idx) = ...
                Pmulti_cond(il,b,q);

            NuReconnect_per_s(idx) = ...
                NuReconnection(il,b,q);

            NuHOBase_per_s(idx) = ...
                NuHO_base(il,b,q);

            NuHOCorrected_per_s(idx) = ...
                NuHO_corrected(il,b,q);

            if isfinite(AssignEmp_s(idx)) ...
                    && AssignEmp_s(idx) > 0

                ErrorTheory_pct(idx) = ...
                    100*abs( ...
                    AssignTheory_s(idx) ...
                    - AssignEmp_s(idx)) ...
                    / AssignEmp_s(idx);
            end
        end
    end
end

Results = table( ...
    Lambda,Time_s,Latitude_deg, ...
    MeanVisibleEmp,MeanVisibleTheory, ...
    OutageEmp,OutageTheory, ...
    AssignEmp_s,AssignTheory_s, ...
    TvisTheory_s,RhoSat_per_km2, ...
    PmultiGivenCovered, ...
    NuReconnect_per_s, ...
    NuHOBase_per_s, ...
    NuHOCorrected_per_s, ...
    ErrorTheory_pct);

disp(Results);

%% ============================================================
% 11. DIAGNOSTICS GLOBAUX
%% ============================================================

valid_assign = ...
    isfinite(AssignEmp_s) ...
    & AssignEmp_s > 0 ...
    & isfinite(AssignTheory_s);

fprintf('\n============================================================\n');
fprintf('T_ASSIGNATION - DELTA UNIFORMITE SPATIALE\n');
fprintf('============================================================\n');
fprintf('Periode orbitale                     : %.2f s\n', ...
    T_orbit);
fprintf('Nombre de classes temporelles        : %d\n', ...
    N_time_bins);

if any(valid_assign)
    fprintf('Erreur relative moyenne T_assign     : %.2f %%\n', ...
        mean(ErrorTheory_pct(valid_assign),'omitnan'));
end

fprintf('============================================================\n');

%% ============================================================
% 12. GRAPHIQUES
%% ============================================================

% Densite de reference pour les cartes.
[~,il_ref] = ...
    min(abs(lambda_values-4e-7));

%% Carte empirique T_assign(t,phi)
figure;

imagesc( ...
    user_lat_deg, ...
    time_phase, ...
    squeeze(MeanAssign_emp(il_ref,:,:)));

set(gca,'YDir','normal');
colorbar;

xlabel('Latitude utilisateur (deg)');
ylabel('Temps dans la periode orbitale (s)');

title(sprintf( ...
    'T_{assign}^{emp}(t,phi), lambda = %.1e', ...
    lambda_values(il_ref)));

%% Carte theorique T_assign(t,phi)
figure;

imagesc( ...
    user_lat_deg, ...
    time_phase, ...
    squeeze(MeanAssign_theory(il_ref,:,:)));

set(gca,'YDir','normal');
colorbar;

xlabel('Latitude utilisateur (deg)');
ylabel('Temps dans la periode orbitale (s)');

title(sprintf( ...
    'T_{assign}^{th}(t,phi), lambda = %.1e', ...
    lambda_values(il_ref)));

%% Carte du nombre moyen visible theorique
figure;

imagesc( ...
    user_lat_deg, ...
    time_phase, ...
    squeeze(MeanVisible_theory(il_ref,:,:)));

set(gca,'YDir','normal');
colorbar;

xlabel('Latitude utilisateur (deg)');
ylabel('Temps dans la periode orbitale (s)');

title(sprintf( ...
    'N_{vis}^{th}(t,phi), lambda = %.1e', ...
    lambda_values(il_ref)));

%% Comparaison a quelques instants
snapshot_idx = ...
    unique(round(linspace(1,N_time_bins,4)));

figure;
hold on;

for s = 1:numel(snapshot_idx)

    b = snapshot_idx(s);

    plot( ...
        user_lat_deg, ...
        squeeze(MeanAssign_emp(il_ref,b,:)), ...
        'o-', ...
        'LineWidth',1.3, ...
        'DisplayName', ...
        sprintf('Emp t=%.0f s',time_phase(b)));

    plot( ...
        user_lat_deg, ...
        squeeze(MeanAssign_theory(il_ref,b,:)), ...
        '--', ...
        'LineWidth',1.3, ...
        'DisplayName', ...
        sprintf('Th t=%.0f s',time_phase(b)));
end

grid on;

xlabel('Latitude utilisateur (deg)');
ylabel('Duree moyenne d''assignation (s)');

title('T_{assign}(t,phi) : empirique / theorie');

legend('Location','best');

%% ============================================================
% 13. SAUVEGARDE
%% ============================================================

save('T_assignation_results.mat', ...
    'Results', ...
    'lambda_values', ...
    'user_lat_deg', ...
    'time_phase','time_edges', ...
    'MeanAssign_emp', ...
    'MeanAssign_theory', ...
    'MeanVisible_emp', ...
    'MeanVisible_theory', ...
    'Outage_emp', ...
    'Outage_theory', ...
    'Tvis_theory', ...
    'RhoSat_theory', ...
    'Pmulti_cond', ...
    'NuReconnection', ...
    'NuHO_base', ...
    'NuHO_corrected', ...
    'R_earth','R','h', ...
    'inc_deg','elevation_min_deg', ...
    'psi_max','A_visible', ...
    'omega_sat','omega_earth', ...
    'T_orbit', ...
    'N_realizations','dt','Tmax', ...
    'N_time_bins');

fprintf('\nResultats sauvegardes dans :\n');
fprintf('  T_assignation_spatial_results.mat\n');

%% ============================================================
% FONCTIONS LOCALES
%% ============================================================

function u = sample_u_spatial(N)

    % Echantillonnage de
    %
    %   f_U(u,0)=|cos(u)|/4,
    %
    % equivalent a une repartition initialement uniforme en
    % surface dans la bande |phi|<=i.

    s = 2*rand(N,1)-1;

    a = asin(s);

    branch = rand(N,1)<0.5;

    u = a;

    u(~branch) = ...
        pi-a(~branch);

    u = mod(u,2*pi);
end

function positions = ...
    walker_delta_positions(R,inc,Omega,u)

    x = ...
        R*( ...
        cos(Omega).*cos(u) ...
        - sin(Omega).*sin(u).*cos(inc));

    y = ...
        R*( ...
        sin(Omega).*cos(u) ...
        + cos(Omega).*sin(u).*cos(inc));

    z = ...
        R*(sin(u).*sin(inc));

    positions = [x y z];
end

function positions = ...
    ground_user_positions(R,lat,lon)

    positions = [ ...
        R*cos(lat).*cos(lon), ...
        R*cos(lat).*sin(lon), ...
        R*sin(lat)];
end

function fphi = ...
    latitude_pdf_spatial(phi,t,inc,omega)

    % Loi instantanee de latitude :
    %
    % f_phi(t,phi)
    %   = cos(phi)/sqrt(sin^2(i)-sin^2(phi))
    %     * [f_U(u_+,t)+f_U(u_-,t)].

    if abs(phi) >= inc

        fphi = 0;

        return;
    end

    x = ...
        sin(phi)/sin(inc);

    x = ...
        max(-1,min(1,x));

    u_plus = ...
        asin(x);

    u_minus = ...
        pi-u_plus;

    fU_plus = ...
        abs(cos(u_plus-omega*t))/4;

    fU_minus = ...
        abs(cos(u_minus-omega*t))/4;

    denominator = ...
        sqrt( ...
        max( ...
        sin(inc)^2 ...
        - sin(phi)^2, ...
        eps));

    fphi = ...
        cos(phi)/denominator ...
        * (fU_plus+fU_minus);
end
