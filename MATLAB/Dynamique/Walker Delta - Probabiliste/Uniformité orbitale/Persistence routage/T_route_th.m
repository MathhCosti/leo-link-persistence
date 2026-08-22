clear; clc; close all;

%% ============================================================
% DUREE MOYENNE THEORIQUE DE ROUTE PAR QUADRATURE DETERMINISTE
% WALKER DELTA - UNIFORMITE ORBITALE
%
% Objectif :
%   - ne plus utiliser H_mean dans la formule de T_route ;
%   - ne plus tirer aleatoirement les positions pour obtenir H ;
%   - calculer deterministement la loi P(H=h) ;
%   - calculer directement :
%
%       E[T_route] = somme_q w_q T_route(q)
%
%     ou chaque point q de quadrature correspond a une geometrie
%     (u_A,u_B,Delta_lambda).
%
% Modele Delta orbital :
%
%   u_A,u_B ~ U(0,2*pi)
%   phi = asin(sin(i)*sin(u))
%   Delta_lambda ~ U(0,pi)
%
% L'utilisation de [0,pi] pour Delta_lambda est exacte ici car
% gamma ne depend de Delta_lambda que par cos(Delta_lambda).
%
% Pour chaque geometrie :
%
%   cos(gamma) =
%       sin(phi_A)sin(phi_B)
%       + cos(phi_A)cos(phi_B)cos(Delta_lambda)
%
%   H = max(1,ceil(gamma/alpha_max))
%
%   alpha_max = 2 asin(dmax/(2R))
%
% Puis les H liens sont repartis lineairement entre phi_A et phi_B.
%
% Pour chaque lien j :
%
%   p_break(phi_j)
%       = (4*v_orb*Delta_t)/(pi*dmax)
%         * sqrt(sin(i)^2-sin(phi_j)^2)/cos(phi_j)
%
%   beta_ISL(phi_j)
%       = -log(1-p_break(phi_j))/Delta_t
%
% Enfin :
%
%   T_route(q) =
%       1 / [ beta_assign(phi_A)
%             + beta_assign(phi_B)
%             + somme_j beta_ISL(phi_j) ]
%
% et la moyenne finale est obtenue directement par quadrature.
%
% IMPORTANT :
%   Cette moyenne globale suppose que les extremites de la route suivent
%   la loi spatiale Delta orbitale. Comme dans le modele precedent,
%   on assimile la latitude du satellite assigne a celle de l'utilisateur
%   pour recuperer T_assign(phi).
%% ============================================================

%% ============================================================
% 0. REGLAGES
%% ============================================================

assignment_source = 'theoretical';  % 'theoretical' ou 'empirical'

% Ordres de quadrature.
% 40 x 40 x 60 = 96 000 geometries deterministes.
Nu = 40;
Ndlon = 60;

Delta_t_break = 10;    % [s]
mu = 398600;           % [km^3/s^2]

%% ============================================================
% 1. FICHIERS
%% ============================================================

script_dir = fileparts(mfilename('fullpath'));

assignment_file = fullfile(script_dir,'T_assignation_results.mat');
shortest_path_file = fullfile(script_dir,'H_jumps_results.mat');

% Resultats empiriques de T_route utilises uniquement pour les
% remplacements successifs H_emp puis beta_link_emp.
empirical_route_file = fullfile(script_dir,'T_route_emp_results.mat');

assert(isfile(assignment_file), ...
    'Fichier d''assignation introuvable : %s',assignment_file);

assert(isfile(shortest_path_file), ...
    'Fichier H_jumps introuvable : %s',shortest_path_file);

assignment_data = load(assignment_file);
path_data = load(shortest_path_file);

has_empirical_route_data = isfile(empirical_route_file);

H_mean_emp = NaN;
beta_link_emp = NaN;

if has_empirical_route_data
    emp_route_data = load(empirical_route_file);

    % Nombre moyen empirique de sauts des routes effectivement observees.
    if isfield(emp_route_data,'H_route_emp_mean')
        H_mean_emp = double(emp_route_data.H_route_emp_mean);
    elseif isfield(path_data,'mean_H')
        H_mean_emp = double(path_data.mean_H);
    end

    % Taux empirique DIRECT de rupture d'un lien de route.
    if isfield(emp_route_data,'beta_link_emp_direct')
        beta_link_emp = double(emp_route_data.beta_link_emp_direct);
    elseif isfield(emp_route_data,'beta_link_emp')
        beta_link_emp = double(emp_route_data.beta_link_emp);
    end
elseif isfield(path_data,'mean_H')
    % On peut au moins recuperer H empirique depuis H_jumps.
    H_mean_emp = double(path_data.mean_H);
end

%% ============================================================
% 1.b PROBABILITE EMPIRIQUE DE CONNEXION SELON GAMMA
%% ============================================================

% On utilise GammaStatistics produit par H_jumps.m :
%
%   P_conn(gamma) = P(H < Inf | gamma)
%
% Cette courbe sert uniquement a conditionner la moyenne theorique
% de T_route sur les geometries pour lesquelles un chemin existe.

% Le conditionnement empirique est conserve uniquement pour comparaison.
% Le conditionnement analytique ajoute plus bas ne depend pas de
% GammaStatistics.

has_empirical_connectivity = false;
has_empirical_H_by_gamma = false;

gamma_conn_deg = [];
Pconn_gamma = [];

% Courbe empirique du nombre de sauts des ROUTES effectivement observees :
%
%   H_route_emp(gamma)
%       = E[H_initial | gamma_ground ~= gamma, route observee]
%
% Elle sera construite directement depuis T_route_emp_results.mat.
gamma_H_emp_deg = [];
H_emp_by_gamma = [];
H_emp_bin_counts = [];

if isfield(path_data,'GammaStatistics')

    GammaStatistics = path_data.GammaStatistics;

    % ----- P_conn(gamma) empirique -----
    required_conn_fields = { ...
        'GammaCenter_deg','ConnectionProbability'};

    has_conn_fields = all(ismember(required_conn_fields, ...
        GammaStatistics.Properties.VariableNames));

    if has_conn_fields
        gamma_conn_deg = double(GammaStatistics.GammaCenter_deg(:));
        Pconn_gamma = double(GammaStatistics.ConnectionProbability(:));

        valid_conn = isfinite(gamma_conn_deg) & isfinite(Pconn_gamma);

        gamma_conn_deg = gamma_conn_deg(valid_conn);
        Pconn_gamma = Pconn_gamma(valid_conn);

        [gamma_conn_deg,idx_conn] = sort(gamma_conn_deg);
        Pconn_gamma = Pconn_gamma(idx_conn);

        Pconn_gamma = min(max(Pconn_gamma,0),1);

        has_empirical_connectivity = numel(gamma_conn_deg) >= 2;
    end

end

%% ============================================================
% 1.c H EMPIRIQUE DES ROUTES EN FONCTION DE GAMMA
%% ============================================================

% IMPORTANT :
% on n'utilise plus GammaStatistics.MeanHopCount, qui correspond aux
% plus courts chemins entre deux satellites generiques connectes.
%
% On construit ici la courbe directement avec les episodes de route :
%
%   route_gamma_ground_deg
%   route_hops_initial
%
% contenus dans T_route_emp_results.mat.

if has_empirical_route_data && ...
        isfield(emp_route_data,'route_gamma_ground_deg') && ...
        isfield(emp_route_data,'route_hops_initial')

    gamma_route_emp = double(emp_route_data.route_gamma_ground_deg(:));
    H_route_emp = double(emp_route_data.route_hops_initial(:));

    valid_route_H = ...
        isfinite(gamma_route_emp) & ...
        isfinite(H_route_emp) & ...
        H_route_emp >= 0;

    gamma_route_emp = gamma_route_emp(valid_route_H);
    H_route_emp = H_route_emp(valid_route_H);

    % Classes angulaires. 10 deg donne en general assez d'episodes
    % par classe tout en conservant la dependance en gamma.
    gamma_H_bin_width_deg = 10;
    gamma_H_edges_deg = ...
        (0:gamma_H_bin_width_deg:180).';

    if gamma_H_edges_deg(end) < 180
        gamma_H_edges_deg(end+1) = 180;
    end

    gamma_H_centers_all = ...
        0.5*(gamma_H_edges_deg(1:end-1)+ ...
             gamma_H_edges_deg(2:end));

    N_H_bins = numel(gamma_H_centers_all);

    H_emp_mean_all = NaN(N_H_bins,1);
    H_emp_bin_counts_all = zeros(N_H_bins,1);

    for jb = 1:N_H_bins

        if jb < N_H_bins
            mask_bin = ...
                gamma_route_emp >= gamma_H_edges_deg(jb) & ...
                gamma_route_emp < gamma_H_edges_deg(jb+1);
        else
            mask_bin = ...
                gamma_route_emp >= gamma_H_edges_deg(jb) & ...
                gamma_route_emp <= gamma_H_edges_deg(jb+1);
        end

        H_emp_bin_counts_all(jb) = nnz(mask_bin);

        if H_emp_bin_counts_all(jb) > 0
            H_emp_mean_all(jb) = mean(H_route_emp(mask_bin));
        end
    end

    valid_H_bins = ...
        isfinite(H_emp_mean_all) & ...
        H_emp_bin_counts_all > 0;

    gamma_H_emp_deg = ...
        gamma_H_centers_all(valid_H_bins);

    H_emp_by_gamma = ...
        H_emp_mean_all(valid_H_bins);

    H_emp_bin_counts = ...
        H_emp_bin_counts_all(valid_H_bins);

    has_empirical_H_by_gamma = ...
        numel(gamma_H_emp_deg) >= 2;

else
    gamma_H_bin_width_deg = NaN;
    gamma_H_edges_deg = [];
end

%% ============================================================
% 2. PARAMETRES DU MODELE
%% ============================================================

assert(isfield(path_data,'R'), ...
    'Le fichier H_jumps ne contient pas R.');
assert(isfield(path_data,'dmax'), ...
    'Le fichier H_jumps ne contient pas dmax.');
assert(isfield(path_data,'inc_deg'), ...
    'Le fichier H_jumps ne contient pas inc_deg.');

R_orbit = double(path_data.R);
dmax = double(path_data.dmax);
inc_deg = double(path_data.inc_deg);
inc = deg2rad(inc_deg);

v_orb = sqrt(mu/R_orbit);

alpha_max = 2*asin(min(1,dmax/(2*R_orbit)));

Hmax = ceil(pi/alpha_max);

%% ============================================================
% 3. DUREE D'ASSIGNATION T_assign(phi)
%% ============================================================

assert(isfield(assignment_data,'lambda_values'), ...
    'Champ lambda_values manquant.');
assert(isfield(assignment_data,'user_lat_deg'), ...
    'Champ user_lat_deg manquant.');

lambda_values = double(assignment_data.lambda_values(:));
user_lat_deg = double(assignment_data.user_lat_deg(:));

switch lower(assignment_source)

    case 'theoretical'
        field_assign = 'MeanAssign_theory';
        assignment_label = 'theorique';

    case 'empirical'
        field_assign = 'MeanAssign_emp';
        assignment_label = 'empirique';

    otherwise
        error('assignment_source doit valoir theoretical ou empirical.');
end

assert(isfield(assignment_data,field_assign), ...
    'Le fichier d''assignation ne contient pas %s.',field_assign);

Tassign_matrix = double(assignment_data.(field_assign));

% Densite correspondant a celle du modele H_jumps.
if isfield(path_data,'lambda')
    lambda_target = double(path_data.lambda);
else
    lambda_target = lambda_values(1);
end

[~,il] = min(abs(lambda_values-lambda_target));

Tassign_lat = Tassign_matrix(il,:).';

% Tri pour interpolation.
[user_lat_deg_sorted,idx_sort] = sort(user_lat_deg);
Tassign_lat_sorted = Tassign_lat(idx_sort);

% Elimination des points non finis.
valid_assign = isfinite(Tassign_lat_sorted) & Tassign_lat_sorted > 0;
lat_assign_deg = user_lat_deg_sorted(valid_assign);
Tassign_values = Tassign_lat_sorted(valid_assign);

assert(numel(lat_assign_deg) >= 2, ...
    'Pas assez de valeurs valides de T_assign pour interpoler.');

lat_min_assign = min(lat_assign_deg);
lat_max_assign = max(lat_assign_deg);

%% ============================================================
% 3.b CONDITIONNEMENT ANALYTIQUE APPROXIMATIF
%% ============================================================

% Approximation :
% pour progresser d'un saut vers la destination, on demande au moins
% un satellite relais dans la moitie "utile" de la calotte ISL.
%
% Aire totale de voisinage ISL :
%   A_cap = 2*pi*R^2*(1-cos(alpha_max))
%
% Aire de progression :
%   A_prog ~= A_cap/2
%
% Sous une approximation de Poisson locale homogene :
%   p_step = 1-exp(-lambda*A_prog)
%
% Pour H sauts, H-1 relais intermediaires sont necessaires :
%   P_conn_ana(H) ~= p_step^(H-1)
%
% Cette fermeture ne depend d'aucune mesure empirique de connectivite.

lambda_conn = lambda_values(il);

A_cap_ISL = ...
    2*pi*R_orbit^2*(1-cos(alpha_max));

forward_area_fraction = 0.5;
A_progress = forward_area_fraction*A_cap_ISL;

mean_forward_relays = ...
    lambda_conn*A_progress;

p_step_analytic = ...
    1-exp(-mean_forward_relays);

%% ============================================================
% 4. QUADRATURE DE GAUSS-LEGENDRE
%% ============================================================

% u_A,u_B uniformes sur [0,2*pi].
[xu,wu] = gauss_legendre(Nu,0,2*pi);

% Delta_lambda peut etre ramene a [0,pi] :
% E[f(cos Delta_lambda)] = (1/pi) integral_0^pi f(cos d) dd.
[xd,wd] = gauss_legendre(Ndlon,0,pi);

% Normalisation des mesures uniformes.
wu = wu/(2*pi);
wd = wd/pi;

% Latitudes Delta associees aux noeuds de quadrature u.
phi_u = asin(sin(inc).*sin(xu));
phi_u_deg = rad2deg(phi_u);

%% ============================================================
% 5. INTEGRATION DETERMINISTE
%% ============================================================

P_H = zeros(Hmax,1);

Troute_weighted_sum = 0;
weight_sum = 0;

H_weighted_sum = 0;
Gamma_weighted_sum = 0;

% Version conditionnee EMPIRIQUEMENT par P_conn(gamma) :
Troute_cond_emp_weighted_sum = 0;
weight_cond_emp_sum = 0;
H_cond_emp_weighted_sum = 0;
Gamma_cond_emp_weighted_sum = 0;
BetaGSL_cond_emp_weighted_sum = 0;
BetaISL_cond_emp_weighted_sum = 0;

% Comparaisons successives SOUS LE MEME CONDITIONNEMENT EMPIRIQUE :
%
%   1) H_th(gamma) + pbreak_th       -> Troute_cond_emp
%   2) H_route,emp(gamma) + pbreak_th
%   3) H_route,emp(gamma) + pbreak_emp
%
% Le poids w_cond_emp est identique dans les trois cas.
Troute_cond_emp_HgammaEmp_weighted_sum = 0;
Troute_cond_emp_HgammaEmp_PbreakEmp_weighted_sum = 0;
Hgamma_emp_under_emp_conditioning_weighted_sum = 0;

% Version conditionnee ANALYTIQUEMENT :
Troute_cond_ana_weighted_sum = 0;
weight_cond_ana_sum = 0;
H_cond_ana_weighted_sum = 0;
Gamma_cond_ana_weighted_sum = 0;
BetaGSL_cond_ana_weighted_sum = 0;
BetaISL_cond_ana_weighted_sum = 0;

% Diagnostics secondaires sous conditionnement analytique.
% Les comparaisons principales d'approximations seront effectuees plus bas
% avec le MEME conditionnement empirique pour toutes les variantes.
%
% H_emp(gamma) = E[H_initial | gamma_ground, route observee],
% construit directement depuis T_route_emp_results.mat.
Troute_cond_ana_HgammaEmp_weighted_sum = 0;
Troute_cond_ana_HgammaEmp_PbreakEmp_weighted_sum = 0;

Hgamma_emp_cond_weighted_sum = 0;

% Diagnostics :
Tassign_A_weighted_sum = 0;
Tassign_B_weighted_sum = 0;
BetaGSL_weighted_sum = 0;
BetaISL_weighted_sum = 0;

for ia = 1:Nu

    phi_A = phi_u(ia);
    phi_A_deg = phi_u_deg(ia);

    % On borne seulement pour l'interpolation si la grille T_assign
    % ne couvre pas exactement +/- inclinaison.
    phi_A_assign_deg = min(max(phi_A_deg,lat_min_assign),lat_max_assign);

    Tassign_A = interp1( ...
        lat_assign_deg,Tassign_values,phi_A_assign_deg,'linear');

    beta_assign_A = 1/Tassign_A;

    for ib = 1:Nu

        phi_B = phi_u(ib);
        phi_B_deg = phi_u_deg(ib);

        phi_B_assign_deg = min(max(phi_B_deg,lat_min_assign),lat_max_assign);

        Tassign_B = interp1( ...
            lat_assign_deg,Tassign_values,phi_B_assign_deg,'linear');

        beta_assign_B = 1/Tassign_B;

        beta_route_GSL = beta_assign_A + beta_assign_B;

        % Poids commun aux deux phases orbitales.
        w_u = wu(ia)*wu(ib);

        for id = 1:Ndlon

            Delta_lambda = xd(id);

            %% Separation angulaire
            cos_gamma = ...
                sin(phi_A)*sin(phi_B) ...
                + cos(phi_A)*cos(phi_B)*cos(Delta_lambda);

            cos_gamma = max(-1,min(1,cos_gamma));
            gamma = acos(cos_gamma);
            gamma_deg = rad2deg(gamma);

            %% Nombre de sauts theorique
            H = max(1,ceil(gamma/alpha_max));
            H = min(H,Hmax);

            %% Latitudes des H liens
            phi_nodes = linspace(phi_A,phi_B,H+1);

            phi_links = ...
                0.5*(phi_nodes(1:end-1)+phi_nodes(2:end));

            %% Probabilites de rupture ISL
            latitude_factor = ...
                sqrt(max(sin(inc)^2-sin(phi_links).^2,0)) ...
                ./ max(cos(phi_links),eps);

            p_break_links = ...
                (4*v_orb*Delta_t_break/(pi*dmax)) ...
                .* latitude_factor;

            p_break_links = ...
                min(max(p_break_links,0),1-eps);

            beta_break_links = ...
                -log1p(-p_break_links)/Delta_t_break;

            beta_route_ISL = sum(beta_break_links);

            % Taux theorique moyen par lien pour CETTE geometrie.
            beta_link_theory_local = ...
                beta_route_ISL/max(H,1);

            %% Duree conditionnelle de la route
            beta_route_total = ...
                beta_route_GSL + beta_route_ISL;

            Troute = 1/beta_route_total;

            %% ------------------------------------------------
            % Comparaison 1 : remplacement de H par H_emp(gamma)
            %
            % On utilise la moyenne empirique du plus court chemin
            % CONDITIONNELLE a la separation angulaire gamma :
            %
            %   H_emp(gamma) = E[H | gamma, route connectee].
            %
            % Cela conserve la dependance geometrique et evite d'injecter
            % une moyenne globale incompatible avec le conditionnement.
            %% ------------------------------------------------
            if has_empirical_H_by_gamma

                % On borne gamma a la plage effectivement observee
                % afin de ne pas extrapoler artificiellement H.
                gamma_for_H = min(max(gamma_deg, ...
                    min(gamma_H_emp_deg)),max(gamma_H_emp_deg));

                H_emp_gamma = interp1( ...
                    gamma_H_emp_deg,H_emp_by_gamma,gamma_for_H, ...
                    'linear');

                % Une route sans ISL peut avoir H=0.
                H_emp_gamma = max(H_emp_gamma,0);

                % On conserve le taux theorique moyen par lien local.
                beta_route_ISL_HgammaEmp = ...
                    H_emp_gamma*beta_link_theory_local;

                Troute_HgammaEmp = ...
                    1/(beta_route_GSL+beta_route_ISL_HgammaEmp);

            else
                H_emp_gamma = NaN;
                Troute_HgammaEmp = NaN;
            end

            %% ------------------------------------------------
            % Comparaison 2 :
            % H_emp(gamma) + beta_link empirique direct
            %% ------------------------------------------------
            if isfinite(H_emp_gamma) && ...
                    isfinite(beta_link_emp) && beta_link_emp >= 0

                beta_route_ISL_HgammaEmp_PbreakEmp = ...
                    H_emp_gamma*beta_link_emp;

                Troute_HgammaEmp_PbreakEmp = ...
                    1/(beta_route_GSL+ ...
                    beta_route_ISL_HgammaEmp_PbreakEmp);

            else
                Troute_HgammaEmp_PbreakEmp = NaN;
            end

            %% Poids de quadrature
            w = w_u*wd(id);

            %% Probabilites de connexion pour cette geometrie

            % ----- Conditionnement analytique approximatif -----
            %
            % H=1 : aucun relais intermediaire necessaire.
            % H>1 : H-1 succes de progression.
            P_conn_ana = p_step_analytic^(max(H-1,0));
            P_conn_ana = min(max(P_conn_ana,0),1);

            w_cond_ana = w*P_conn_ana;

            % ----- Conditionnement empirique, pour comparaison -----
            if has_empirical_connectivity
                P_conn_emp = interp1( ...
                    gamma_conn_deg,Pconn_gamma,gamma_deg, ...
                    'linear','extrap');

                P_conn_emp = min(max(P_conn_emp,0),1);
                w_cond_emp = w*P_conn_emp;
            else
                P_conn_emp = NaN;
                w_cond_emp = 0;
            end

            %% Accumulation brute
            P_H(H) = P_H(H) + w;

            Troute_weighted_sum = ...
                Troute_weighted_sum + w*Troute;

            H_weighted_sum = ...
                H_weighted_sum + w*H;

            Gamma_weighted_sum = ...
                Gamma_weighted_sum + w*gamma;

            Tassign_A_weighted_sum = ...
                Tassign_A_weighted_sum + w*Tassign_A;

            Tassign_B_weighted_sum = ...
                Tassign_B_weighted_sum + w*Tassign_B;

            BetaGSL_weighted_sum = ...
                BetaGSL_weighted_sum + w*beta_route_GSL;

            BetaISL_weighted_sum = ...
                BetaISL_weighted_sum + w*beta_route_ISL;

            %% Accumulation conditionnee analytique
            Troute_cond_ana_weighted_sum = ...
                Troute_cond_ana_weighted_sum + w_cond_ana*Troute;

            H_cond_ana_weighted_sum = ...
                H_cond_ana_weighted_sum + w_cond_ana*H;

            Gamma_cond_ana_weighted_sum = ...
                Gamma_cond_ana_weighted_sum + w_cond_ana*gamma;

            BetaGSL_cond_ana_weighted_sum = ...
                BetaGSL_cond_ana_weighted_sum + ...
                w_cond_ana*beta_route_GSL;

            BetaISL_cond_ana_weighted_sum = ...
                BetaISL_cond_ana_weighted_sum + ...
                w_cond_ana*beta_route_ISL;

            % Remplacements empiriques, avec LE MEME poids de
            % conditionnement analytique.
            if isfinite(Troute_HgammaEmp)
                Troute_cond_ana_HgammaEmp_weighted_sum = ...
                    Troute_cond_ana_HgammaEmp_weighted_sum + ...
                    w_cond_ana*Troute_HgammaEmp;

                Hgamma_emp_cond_weighted_sum = ...
                    Hgamma_emp_cond_weighted_sum + ...
                    w_cond_ana*H_emp_gamma;
            end

            if isfinite(Troute_HgammaEmp_PbreakEmp)
                Troute_cond_ana_HgammaEmp_PbreakEmp_weighted_sum = ...
                    Troute_cond_ana_HgammaEmp_PbreakEmp_weighted_sum + ...
                    w_cond_ana*Troute_HgammaEmp_PbreakEmp;
            end

            %% Accumulation conditionnee empirique
            if has_empirical_connectivity
                Troute_cond_emp_weighted_sum = ...
                    Troute_cond_emp_weighted_sum + w_cond_emp*Troute;

                H_cond_emp_weighted_sum = ...
                    H_cond_emp_weighted_sum + w_cond_emp*H;

                Gamma_cond_emp_weighted_sum = ...
                    Gamma_cond_emp_weighted_sum + w_cond_emp*gamma;

                BetaGSL_cond_emp_weighted_sum = ...
                    BetaGSL_cond_emp_weighted_sum + ...
                    w_cond_emp*beta_route_GSL;

                BetaISL_cond_emp_weighted_sum = ...
                    BetaISL_cond_emp_weighted_sum + ...
                    w_cond_emp*beta_route_ISL;

                % ---------------------------------------------
                % Remplacement 1 sous conditionnement empirique :
                % H_route,emp(gamma), pbreak theorique
                % ---------------------------------------------
                if isfinite(Troute_HgammaEmp)
                    Troute_cond_emp_HgammaEmp_weighted_sum = ...
                        Troute_cond_emp_HgammaEmp_weighted_sum + ...
                        w_cond_emp*Troute_HgammaEmp;

                    Hgamma_emp_under_emp_conditioning_weighted_sum = ...
                        Hgamma_emp_under_emp_conditioning_weighted_sum + ...
                        w_cond_emp*H_emp_gamma;
                end

                % ---------------------------------------------
                % Remplacement 2 sous conditionnement empirique :
                % H_route,emp(gamma), pbreak empirique
                % ---------------------------------------------
                if isfinite(Troute_HgammaEmp_PbreakEmp)
                    Troute_cond_emp_HgammaEmp_PbreakEmp_weighted_sum = ...
                        Troute_cond_emp_HgammaEmp_PbreakEmp_weighted_sum + ...
                        w_cond_emp*Troute_HgammaEmp_PbreakEmp;
                end
            end

            weight_sum = weight_sum + w;
            weight_cond_ana_sum = weight_cond_ana_sum + w_cond_ana;

            if has_empirical_connectivity
                weight_cond_emp_sum = weight_cond_emp_sum + w_cond_emp;
            end
        end
    end
end

%% ============================================================
% 6. NORMALISATION ET RESULTATS
%% ============================================================

P_H = P_H/weight_sum;

Troute_mean = Troute_weighted_sum/weight_sum;
H_mean_from_distribution = H_weighted_sum/weight_sum;
Gamma_mean_deg = rad2deg(Gamma_weighted_sum/weight_sum);

MeanTassign_A = Tassign_A_weighted_sum/weight_sum;
MeanTassign_B = Tassign_B_weighted_sum/weight_sum;

MeanBetaRouteGSL = BetaGSL_weighted_sum/weight_sum;
MeanBetaRouteISL = BetaISL_weighted_sum/weight_sum;

%% Version conditionnee ANALYTIQUEMENT
assert(weight_cond_ana_sum > 0, ...
    'Poids conditionne analytique nul.');

Troute_mean_cond_ana = ...
    Troute_cond_ana_weighted_sum/weight_cond_ana_sum;

H_mean_cond_ana = ...
    H_cond_ana_weighted_sum/weight_cond_ana_sum;

Gamma_mean_cond_ana_deg = ...
    rad2deg(Gamma_cond_ana_weighted_sum/weight_cond_ana_sum);

MeanBetaRouteGSL_cond_ana = ...
    BetaGSL_cond_ana_weighted_sum/weight_cond_ana_sum;

MeanBetaRouteISL_cond_ana = ...
    BetaISL_cond_ana_weighted_sum/weight_cond_ana_sum;

P_connected_analytic = ...
    weight_cond_ana_sum/weight_sum;

conditioning_gain_ana_s = ...
    Troute_mean_cond_ana-Troute_mean;

conditioning_gain_ana_percent = ...
    100*conditioning_gain_ana_s/Troute_mean;

%% ============================================================
% REMPLACEMENTS SUCCESSIFS PAR LES VALEURS EMPIRIQUES
%% ============================================================

Troute_mean_cond_ana_HgammaEmp = NaN;
Troute_mean_cond_ana_HgammaEmp_PbreakEmp = NaN;
Hgamma_emp_mean_under_ana_conditioning = NaN;

if has_empirical_H_by_gamma
    Troute_mean_cond_ana_HgammaEmp = ...
        Troute_cond_ana_HgammaEmp_weighted_sum/ ...
        weight_cond_ana_sum;

    Hgamma_emp_mean_under_ana_conditioning = ...
        Hgamma_emp_cond_weighted_sum/weight_cond_ana_sum;
end

if has_empirical_H_by_gamma && ...
        isfinite(beta_link_emp) && beta_link_emp >= 0

    Troute_mean_cond_ana_HgammaEmp_PbreakEmp = ...
        Troute_cond_ana_HgammaEmp_PbreakEmp_weighted_sum/ ...
        weight_cond_ana_sum;
end

%% Version conditionnee EMPIRIQUEMENT, pour comparaison
Troute_mean_cond_emp = NaN;
H_mean_cond_emp = NaN;
Gamma_mean_cond_emp_deg = NaN;
MeanBetaRouteGSL_cond_emp = NaN;
MeanBetaRouteISL_cond_emp = NaN;
P_connected_empirical = NaN;
conditioning_gain_emp_s = NaN;
conditioning_gain_emp_percent = NaN;

Troute_mean_cond_emp_HgammaEmp = NaN;
Troute_mean_cond_emp_HgammaEmp_PbreakEmp = NaN;
Hgamma_emp_mean_under_emp_conditioning = NaN;

if has_empirical_connectivity && weight_cond_emp_sum > 0

    Troute_mean_cond_emp = ...
        Troute_cond_emp_weighted_sum/weight_cond_emp_sum;

    H_mean_cond_emp = ...
        H_cond_emp_weighted_sum/weight_cond_emp_sum;

    Gamma_mean_cond_emp_deg = ...
        rad2deg(Gamma_cond_emp_weighted_sum/weight_cond_emp_sum);

    MeanBetaRouteGSL_cond_emp = ...
        BetaGSL_cond_emp_weighted_sum/weight_cond_emp_sum;

    MeanBetaRouteISL_cond_emp = ...
        BetaISL_cond_emp_weighted_sum/weight_cond_emp_sum;

    P_connected_empirical = ...
        weight_cond_emp_sum/weight_sum;

    conditioning_gain_emp_s = ...
        Troute_mean_cond_emp-Troute_mean;

    conditioning_gain_emp_percent = ...
        100*conditioning_gain_emp_s/Troute_mean;

    if has_empirical_H_by_gamma
        Troute_mean_cond_emp_HgammaEmp = ...
            Troute_cond_emp_HgammaEmp_weighted_sum / ...
            weight_cond_emp_sum;

        Hgamma_emp_mean_under_emp_conditioning = ...
            Hgamma_emp_under_emp_conditioning_weighted_sum / ...
            weight_cond_emp_sum;
    end

    if has_empirical_H_by_gamma && ...
            isfinite(beta_link_emp) && beta_link_emp >= 0
        Troute_mean_cond_emp_HgammaEmp_PbreakEmp = ...
            Troute_cond_emp_HgammaEmp_PbreakEmp_weighted_sum / ...
            weight_cond_emp_sum;
    end
end

% Controle de la loi de H.
P_H_sum = sum(P_H);

h_values = (1:Hmax).';
cdf_H = cumsum(P_H);

HopDistribution = table( ...
    h_values,P_H,cdf_H, ...
    'VariableNames',{'HopCount','Probability','CDF'});

%% ============================================================
% 7. COMPARAISON AVEC L'ANCIENNE APPROXIMATION PAR H MOYEN
%% ============================================================

% Cette comparaison sert uniquement a quantifier la difference entre
% E[1/beta_route] et 1/E[beta_route].
%
% Elle ne remplace pas le resultat principal.
Troute_from_mean_rates = ...
    1/(MeanBetaRouteGSL + MeanBetaRouteISL);

relative_difference_percent = ...
    100*(Troute_mean-Troute_from_mean_rates)/Troute_mean;

%% ============================================================
% 8. AFFICHAGE
%% ============================================================

fprintf('\n============================================================\n');
fprintf('T_ROUTE MOYEN - QUADRATURE DETERMINISTE DELTA ORBITALE\n');
fprintf('============================================================\n');

fprintf('Ordre quadrature u                  : %d\n',Nu);
fprintf('Ordre quadrature Delta_lambda       : %d\n',Ndlon);
fprintf('Nombre de geometries                : %d\n',Nu*Nu*Ndlon);

fprintf('\nInclinaison                         : %.2f deg\n',inc_deg);
fprintf('dmax                                : %.2f km\n',dmax);
fprintf('alpha_max                           : %.4f deg\n',rad2deg(alpha_max));
fprintf('Source T_assign                     : %s\n',assignment_label);
fprintf('lambda                              : %.4e km^-2\n',lambda_values(il));

fprintf('\nSomme P(H=h)                        : %.12f\n',P_H_sum);
fprintf('E[H] par quadrature                 : %.6f\n',H_mean_from_distribution);
fprintf('E[gamma]                            : %.6f deg\n',Gamma_mean_deg);

fprintf('\nE[T_assign,A]                       : %.6f s\n',MeanTassign_A);
fprintf('E[T_assign,B]                       : %.6f s\n',MeanTassign_B);

fprintf('\nE[beta_route,GSL]                   : %.10e s^-1\n', ...
    MeanBetaRouteGSL);
fprintf('E[beta_route,ISL]                   : %.10e s^-1\n', ...
    MeanBetaRouteISL);

fprintf('\nRESULTATS PRINCIPAUX :\n');
fprintf('E[T_route] brute                    : %.6f s\n', ...
    Troute_mean);

fprintf('\n--- CONDITIONNEMENT ANALYTIQUE ---\n');
fprintf('Aire de progression                 : %.6e km^2\n', ...
    A_progress);
fprintf('Nombre moyen de relais utiles       : %.6f\n', ...
    mean_forward_relays);
fprintf('p_step analytique                   : %.6f\n', ...
    p_step_analytic);
fprintf('P(route existe) analytique          : %.6f\n', ...
    P_connected_analytic);
fprintf('E[T_route | route existe] analytique: %.6f s\n', ...
    Troute_mean_cond_ana);
fprintf('Gain analytique                     : %+.6f s\n', ...
    conditioning_gain_ana_s);
fprintf('Gain analytique relatif             : %+.3f %%\n', ...
    conditioning_gain_ana_percent);
fprintf('E[H | route] analytique             : %.6f\n', ...
    H_mean_cond_ana);
fprintf('E[gamma | route] analytique         : %.6f deg\n', ...
    Gamma_mean_cond_ana_deg);
fprintf('E[beta_GSL | route] analytique      : %.10e s^-1\n', ...
    MeanBetaRouteGSL_cond_ana);
fprintf('E[beta_ISL | route] analytique      : %.10e s^-1\n', ...
    MeanBetaRouteISL_cond_ana);

fprintf('\n--- COMPARAISON DES APPROXIMATIONS ---\n');

if has_empirical_connectivity
    fprintf(['Meme conditionnement empirique P_conn(gamma) ' ...
             'pour les trois lignes :\n']);

    fprintf(['1) H_th + pbreak_th                  : ' ...
             '%.6f s\n'], ...
        Troute_mean_cond_emp);

    if has_empirical_H_by_gamma
        fprintf('   E[H_route,emp(gamma) | cond.emp]   : %.6f\n', ...
            Hgamma_emp_mean_under_emp_conditioning);

        fprintf(['2) H_route,emp(gamma) + pbreak_th    : ' ...
                 '%.6f s\n'], ...
            Troute_mean_cond_emp_HgammaEmp);
    else
        fprintf('2) H_route,emp(gamma) non disponible.\n');
    end

    if has_empirical_H_by_gamma && isfinite(beta_link_emp)
        fprintf('   beta_link_emp direct               : %.10e s^-1\n', ...
            beta_link_emp);

        fprintf(['3) H_route,emp(gamma) + pbreak_emp   : ' ...
                 '%.6f s\n'], ...
            Troute_mean_cond_emp_HgammaEmp_PbreakEmp);
    else
        fprintf('3) pbreak_emp non disponible.\n');
    end

    if has_empirical_route_data && isfield(emp_route_data,'Troute_emp_mean')
        fprintf(['\nReference simulation complete          : ' ...
                 '%.6f s\n'], ...
            double(emp_route_data.Troute_emp_mean));
    end
else
    fprintf(['Conditionnement empirique indisponible : ' ...
             'comparaison successive impossible.\n']);
end

fprintf('\nPour comparaison seulement :\n');
fprintf('1 / E[beta_route]                   : %.6f s\n', ...
    Troute_from_mean_rates);
fprintf('Ecart relatif                       : %.3f %%\n', ...
    relative_difference_percent);

fprintf('============================================================\n');

disp(HopDistribution);

%% ============================================================
% 9. PROBABILITE DE CONNEXION ET DISTRIBUTION THEORIQUE DE H
%% ============================================================

%% Comparaison P_conn analytique / empirique selon gamma
gamma_plot_deg = linspace(0,180,361).';
H_plot = max(1,ceil(deg2rad(gamma_plot_deg)/alpha_max));
Pconn_analytic_plot = ...
    p_step_analytic.^(max(H_plot-1,0));

figure;
plot(gamma_plot_deg,Pconn_analytic_plot,'LineWidth',1.7, ...
    'DisplayName','Analytique approx.');
hold on;

if has_empirical_connectivity
    plot(gamma_conn_deg,Pconn_gamma,'o-','LineWidth',1.3, ...
        'DisplayName','Empirique');
end

grid on;
xlabel('Separation angulaire gamma (deg)');
ylabel('P_{conn}(gamma)');
title('Probabilite de connexion');
ylim([0 1.05]);
legend('Location','best');
hold off;

%% Comparaison H theorique / H empirique conditionnel a gamma
if has_empirical_H_by_gamma
    figure;

    gamma_H_plot = linspace( ...
        min(gamma_H_emp_deg),max(gamma_H_emp_deg),400).';

    H_th_plot = max(1,ceil(deg2rad(gamma_H_plot)/alpha_max));

    H_emp_plot = interp1( ...
        gamma_H_emp_deg,H_emp_by_gamma,gamma_H_plot, ...
        'linear','extrap');

    plot(gamma_H_plot,H_th_plot,'LineWidth',1.6, ...
        'DisplayName','H theorique');
    hold on;
    plot(gamma_H_emp_deg,H_emp_by_gamma,'o-','LineWidth',1.5, ...
        'DisplayName','E[H_{route,emp} | gamma]');

    grid on;
    xlabel('Separation angulaire gamma (deg)');
    ylabel('Nombre de sauts');
    title('Nombre de sauts des routes selon la separation angulaire');
    legend('Location','best');
    hold off;

    figure;
    bar(gamma_H_emp_deg,H_emp_bin_counts);
    grid on;
    xlabel('Separation angulaire gamma (deg)');
    ylabel('Nombre d''episodes');
    title('Nombre d''episodes utilises pour H_{route,emp}(gamma)');
end

%% Comparaison propre des approximations
if has_empirical_connectivity
    figure;

    values_compare = Troute_mean_cond_emp;
    labels_compare = {'Cond. emp + H_{th} + pbreak_{th}'};

    if isfinite(Troute_mean_cond_emp_HgammaEmp)
        values_compare(end+1) = Troute_mean_cond_emp_HgammaEmp;
        labels_compare{end+1} = ...
            'Cond. emp + H_{route,emp}(gamma) + pbreak_{th}';
    end

    if isfinite(Troute_mean_cond_emp_HgammaEmp_PbreakEmp)
        values_compare(end+1) = ...
            Troute_mean_cond_emp_HgammaEmp_PbreakEmp;
        labels_compare{end+1} = ...
            'Cond. emp + H_{route,emp}(gamma) + pbreak_{emp}';
    end

    if has_empirical_route_data && isfield(emp_route_data,'Troute_emp_mean')
        values_compare(end+1) = double(emp_route_data.Troute_emp_mean);
        labels_compare{end+1} = 'Simulation empirique';
    end

    bar(values_compare);
    set(gca,'XTick',1:numel(labels_compare));
    set(gca,'XTickLabel',labels_compare);
    xtickangle(18);

    ylabel('T_{route} moyen (s)');
    title('Impact successif des approximations');
    grid on;
end

figure;

bar(h_values,P_H);

grid on;
xlabel('Nombre de sauts H');
ylabel('P(H=h)');
title('Distribution theorique de H - quadrature Delta orbitale');

xlim([0.5 Hmax+0.5]);

%% CDF
figure;

stairs(h_values,cdf_H,'LineWidth',1.6);

grid on;
xlabel('Nombre de sauts h');
ylabel('P(H \leq h)');
title('CDF theorique de H - quadrature Delta orbitale');

xlim([0.5 Hmax+0.5]);
ylim([0 1.02]);

%% ============================================================
% 10. TABLE GLOBALE
%% ============================================================

Results = table( ...
    lambda_values(il),inc_deg,dmax,R_orbit, ...
    Nu,Ndlon, ...
    H_mean_from_distribution,Gamma_mean_deg, ...
    Troute_mean, ...
    A_progress,mean_forward_relays,p_step_analytic, ...
    P_connected_analytic,Troute_mean_cond_ana, ...
    H_mean_emp,beta_link_emp, ...
    Hgamma_emp_mean_under_ana_conditioning, ...
    Troute_mean_cond_ana_HgammaEmp, ...
    Troute_mean_cond_ana_HgammaEmp_PbreakEmp, ...
    H_mean_cond_ana,Gamma_mean_cond_ana_deg, ...
    MeanBetaRouteGSL_cond_ana,MeanBetaRouteISL_cond_ana, ...
    P_connected_empirical,Troute_mean_cond_emp, ...
    Troute_mean_cond_emp_HgammaEmp, ...
    Troute_mean_cond_emp_HgammaEmp_PbreakEmp, ...
    Hgamma_emp_mean_under_emp_conditioning, ...
    H_mean_cond_emp,Gamma_mean_cond_emp_deg, ...
    MeanBetaRouteGSL_cond_emp,MeanBetaRouteISL_cond_emp, ...
    'VariableNames',{ ...
    'Lambda','Inclination_deg','Dmax_km','OrbitalRadius_km', ...
    'QuadratureOrderU','QuadratureOrderDeltaLongitude', ...
    'MeanHopCountRaw','MeanAngularSeparationRaw_deg', ...
    'MeanRouteLifetimeRaw_s', ...
    'ForwardArea_km2','MeanForwardRelays','AnalyticStepProbability', ...
    'AnalyticConnectionProbability','MeanRouteLifetimeAnalyticCond_s', ...
    'EmpiricalGlobalMeanHopCount','EmpiricalBetaLink_per_s', ...
    'MeanRouteEmpiricalHopCountGammaUnderAnalyticConditioning', ...
    'MeanRouteLifetimeAnalyticCond_HrouteGammaEmp_s', ...
    'MeanRouteLifetimeAnalyticCond_HrouteGammaEmp_PbreakEmp_s', ...
    'MeanHopCountAnalyticCond','MeanAngularSeparationAnalyticCond_deg', ...
    'MeanBetaGSLAnalyticCond_per_s','MeanBetaISLAnalyticCond_per_s', ...
    'EmpiricalConnectionProbability','MeanRouteLifetimeEmpiricalCond_Hth_PbreakTh_s', ...
    'MeanRouteLifetimeEmpiricalCond_HrouteEmp_PbreakTh_s', ...
    'MeanRouteLifetimeEmpiricalCond_HrouteEmp_PbreakEmp_s', ...
    'MeanRouteEmpHopCountUnderEmpiricalConditioning', ...
    'MeanHopCountTheoreticalUnderEmpiricalConditioning', ...
    'MeanAngularSeparationEmpiricalCond_deg', ...
    'MeanBetaGSLEmpiricalCond_per_s','MeanBetaISLEmpiricalCond_per_s'});

disp(Results);

%% ============================================================
% 11. SAUVEGARDE
%% ============================================================

save('T_route_th_results.mat', ...
    'Results','HopDistribution', ...
    'h_values','P_H','cdf_H', ...
    'Troute_mean', ...
    'Troute_mean_cond_ana','Troute_mean_cond_emp', ...
    'Troute_mean_cond_emp_HgammaEmp', ...
    'Troute_mean_cond_emp_HgammaEmp_PbreakEmp', ...
    'Hgamma_emp_mean_under_emp_conditioning', ...
    'Troute_mean_cond_ana_HgammaEmp', ...
    'Troute_mean_cond_ana_HgammaEmp_PbreakEmp', ...
    'Hgamma_emp_mean_under_ana_conditioning', ...
    'gamma_H_emp_deg','H_emp_by_gamma','H_emp_bin_counts', ...
    'gamma_H_bin_width_deg','gamma_H_edges_deg', ...
    'has_empirical_H_by_gamma', ...
    'H_mean_emp','beta_link_emp', ...
    'H_mean_from_distribution','H_mean_cond_ana','H_mean_cond_emp', ...
    'Gamma_mean_deg','Gamma_mean_cond_ana_deg','Gamma_mean_cond_emp_deg', ...
    'MeanTassign_A','MeanTassign_B', ...
    'MeanBetaRouteGSL','MeanBetaRouteISL', ...
    'MeanBetaRouteGSL_cond_ana','MeanBetaRouteISL_cond_ana', ...
    'MeanBetaRouteGSL_cond_emp','MeanBetaRouteISL_cond_emp', ...
    'P_connected_analytic','P_connected_empirical', ...
    'A_cap_ISL','A_progress','forward_area_fraction', ...
    'mean_forward_relays','p_step_analytic', ...
    'conditioning_gain_ana_s','conditioning_gain_ana_percent', ...
    'conditioning_gain_emp_s','conditioning_gain_emp_percent', ...
    'has_empirical_connectivity','gamma_conn_deg','Pconn_gamma', ...
    'Troute_from_mean_rates', ...
    'relative_difference_percent', ...
    'Nu','Ndlon', ...
    'lambda_values','lambda_target','il', ...
    'inc_deg','inc','dmax','R_orbit','v_orb', ...
    'alpha_max','Delta_t_break', ...
    'assignment_source');

fprintf('\nResultats sauvegardes dans :\n');
fprintf('  T_route_th_results.mat\n');

%% ============================================================
% FONCTION LOCALE : QUADRATURE DE GAUSS-LEGENDRE
%% ============================================================

function [x,w] = gauss_legendre(n,a,b)

    % Noeuds et poids de Gauss-Legendre sur [a,b].
    %
    % Methode de Golub-Welsch, sans toolbox specifique.

    k = (1:n-1).';

    beta = k ./ sqrt(4*k.^2-1);

    J = diag(beta,1) + diag(beta,-1);

    [V,D] = eig(J);

    x0 = diag(D);

    [x0,idx] = sort(x0);
    V = V(:,idx);

    w0 = 2*(V(1,:).^2).';

    % Passage de [-1,1] a [a,b].
    x = (a+b)/2 + (b-a)/2*x0;
    w = (b-a)/2*w0;
end
