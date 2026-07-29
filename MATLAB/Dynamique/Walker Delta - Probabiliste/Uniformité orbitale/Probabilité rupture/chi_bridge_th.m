%% chi_bridge_uniformite_orbitale.m
% Approximation analytique de chi_bridge et du nombre moyen de ponts
% par composante pour un Walker Delta à uniformité orbitale.
%
% Hypothèses :
% - argument de latitude uniforme ;
% - RAAN uniforme ;
% - inclinaison commune i ;
% - approximation locale plane ;
% - un lien est localement critique s'il n'a aucun voisin commun.
%
% La densité orbitale idéale diverge aux latitudes de retournement ±i.
% Une régularisation epsilon_lat_deg est donc introduite.

clearvars;
clc; close all;

%% Paramètres
R_earth = 6371;       % km
h       = 550;        % km
R       = R_earth + h;

N       = 204;        % nombre de satellites
i_deg   = 58;         % inclinaison orbitale en degrés
dmax    = 1500;       % portée de communication en km

% Fichiers contenant les grandeurs analytiques du Walker Delta
script_dir = fileparts(mfilename('fullpath'));
moments_file = fullfile(script_dir, '..', 'Paramètres', 'betti_results.mat');
plink_file = fullfile(script_dir, '..', 'Paramètres', 'plink_results.mat');

% Si true, le script utilise :
% beta0 = C_macro + E[N1] + E[N2] + E[N3]
use_N123_file = true;

% Régularisation près de ±i
epsilon_lat_deg = 0.10;

% Résolution numérique
n_phi = 4000;
n_r   = 1500;

%% Vérifications
validateattributes(R, {'numeric'}, {'scalar','real','finite','positive'});
validateattributes(N, {'numeric'}, {'scalar','integer','>=',2});
validateattributes(i_deg, {'numeric'}, {'scalar','real','>',0,'<=',90});
validateattributes(dmax, {'numeric'}, {'scalar','real','positive'});

if dmax >= 2*R
    error('dmax doit être strictement inférieur à 2R.');
end

%% Densité orbitale locale
incl = deg2rad(i_deg);
eps_lat = deg2rad(epsilon_lat_deg);
phi_max = incl - eps_lat;

if phi_max <= 0
    error('epsilon_lat_deg est trop grand.');
end

phi = linspace(-phi_max, phi_max, n_phi).';

denom = sqrt(max(sin(incl)^2 - sin(phi).^2, realmin));
lambda_orb = N ./ (2*pi^2*R^2 .* denom);  % satellites/km^2

% Élément de surface intégré sur la longitude
surface_weight = 2*pi*R^2 .* cos(phi);

% Renormalisation pour conserver exactement N satellites
N_before_norm = trapz(phi, lambda_orb .* surface_weight);
lambda_orb = lambda_orb .* (N / N_before_norm);
N_check = trapz(phi, lambda_orb .* surface_weight);

%% Géométrie de communication
alpha_max = 2 * asin(dmax / (2*R));
A_comm_sphere = 2*pi*R^2 * (1 - cos(alpha_max));

r = linspace(0, dmax, n_r).';

A_inter = ...
    2*dmax^2 .* acos(min(max(r ./ (2*dmax), -1), 1)) ...
    - 0.5 .* r .* sqrt(max(4*dmax^2 - r.^2, 0));

f_r_given_link = 2*r ./ dmax^2;

%% Fraction locale de liens sans voisin commun
chi_bridge_local = zeros(size(phi));

for k = 1:numel(phi)
    p_no_common = exp(-lambda_orb(k) .* A_inter);
    chi_bridge_local(k) = trapz(r, p_no_common .* f_r_given_link);
end

chi_bridge_local = min(max(chi_bridge_local, 0), 1);

%% Moyenne conditionnée par l'existence d'un lien
% Le nombre local de liens est proportionnel à lambda_orb^2.
link_weight = lambda_orb.^2 .* surface_weight;

chi_bridge = ...
    trapz(phi, chi_bridge_local .* link_weight) ...
    / trapz(phi, link_weight);

chi_bridge = min(max(chi_bridge, 0), 1);

%% Nombre moyen total de liens à partir de p_link Delta

if ~isfile(plink_file)
    error('Fichier introuvable : %s', plink_file);
end

P = load(plink_file);

if ~isfield(P, 'p_link_delta')
    error('La variable p_link_delta est absente de %s.', plink_file);
end

p_link_delta = double(P.p_link_delta);

% Vérifications de cohérence entre le script et le fichier p_link
if isfield(P,'N') && abs(double(P.N) - N) > 0
    warning('N du script = %d, alors que N du fichier p_link = %.0f.', ...
        N, double(P.N));
end

if isfield(P,'R') && abs(double(P.R) - R) > 1e-9
    warning('R du script = %.6f km, alors que R du fichier p_link = %.6f km.', ...
        R, double(P.R));
end

if isfield(P,'dmax') && abs(double(P.dmax) - dmax) > 1e-9
    warning('dmax du script = %.6f km, alors que dmax du fichier p_link = %.6f km.', ...
        dmax, double(P.dmax));
end

if isfield(P,'inc_deg') && abs(double(P.inc_deg) - i_deg) > 1e-9
    warning(['Inclinaison du script = %.6f deg, alors que celle du ' ...
             'fichier p_link = %.6f deg.'], ...
        i_deg, double(P.inc_deg));
end

validateattributes(p_link_delta, {'numeric'}, ...
    {'scalar','real','finite','>=',0,'<=',1});

% Formule théorique exacte conditionnellement à N :
% E[|E| | N] = C(N,2) p_link_delta
E_edges = N * (N - 1) / 2 * p_link_delta;

% Contrôle avec E_delta si cette variable est enregistrée dans le .mat
if isfield(P, 'E_delta')
    E_edges_mat = double(P.E_delta);

    if abs(E_edges - E_edges_mat) > 1e-8
        warning(['E_edges reconstruit = %.12f, alors que E_delta dans ' ...
                 'le .mat = %.12f.'], E_edges, E_edges_mat);
    end
else
    E_edges_mat = NaN;
end

%% Calcul de beta0 à partir de N1, N2 et N3

% Valeur locale de N1 conservée uniquement à titre de comparaison
E_N1_local = trapz(phi, ...
    lambda_orb .* exp(-lambda_orb .* A_comm_sphere) ...
    .* surface_weight);

if use_N123_file
    if ~isfile(moments_file)
        error('Fichier introuvable : %s', moments_file);
    end

    S = load(moments_file);

    required_fields = {'N1_th_local','N2_th_geom','N3_th_geom'};
    for kk = 1:numel(required_fields)
        if ~isfield(S, required_fields{kk})
            error('Variable manquante dans le .mat : %s', required_fields{kk});
        end
    end

    % Vérification de cohérence entre les paramètres du script et du .mat
    if isfield(S,'N') && abs(double(S.N) - N) > 0
        warning('N du script = %d, alors que N du .mat = %.0f.', ...
            N, double(S.N));
    end

    if isfield(S,'R') && abs(double(S.R) - R) > 1e-9
        warning('R du script = %.6f km, alors que R du .mat = %.6f km.', ...
            R, double(S.R));
    end

    if isfield(S,'dmax') && abs(double(S.dmax) - dmax) > 1e-9
        warning('dmax du script = %.6f km, alors que dmax du .mat = %.6f km.', ...
            dmax, double(S.dmax));
    end

    if isfield(S,'inc_deg') && abs(double(S.inc_deg) - i_deg) > 1e-9
        warning('Inclinaison du script = %.6f deg, alors que celle du .mat = %.6f deg.', ...
            i_deg, double(S.inc_deg));
    end

    E_N1 = double(S.N1_th_local);
    E_N2 = double(S.N2_th_geom);
    E_N3 = double(S.N3_th_geom);

    if isfield(S,'C_macro')
        C_macro = double(S.C_macro);
    else
        C_macro = 2;
    end

    beta0_used = C_macro + E_N1 + E_N2 + E_N3;
    beta0_source = 'C_{macro} + E[N1] + E[N2] + E[N3] depuis le .mat';

    % Contrôle avec beta0_th_123, s'il est présent
    if isfield(S,'beta0_th_123')
        beta0_mat = double(S.beta0_th_123);
        if abs(beta0_used - beta0_mat) > 1e-8
            warning(['beta0 reconstruit = %.12f, mais beta0_th_123 ' ...
                     'dans le .mat = %.12f.'], beta0_used, beta0_mat);
        end
    end
else
    E_N1 = E_N1_local;
    E_N2 = 0;
    E_N3 = 0;
    C_macro = 1;

    beta0_used = C_macro + E_N1;
    beta0_source = 'approximation locale C_{macro} + E[N1]';
end

validateattributes(beta0_used, {'numeric'}, ...
    {'scalar','real','finite','>=',1,'<=',N});

%% Liens et ponts par composante
mean_links_per_component = E_edges / beta0_used;
mean_bridges_per_component = ...
    mean_links_per_component * chi_bridge;

%% Affichage
fprintf('\n');
fprintf('====================================================================\n');
fprintf(' Walker Delta - uniformité orbitale\n');
fprintf('====================================================================\n');
fprintf('N                                        : %d\n', N);
fprintf('R                                        : %.3f km\n', R);
fprintf('Inclinaison                              : %.3f deg\n', i_deg);
fprintf('dmax                                     : %.3f km\n', dmax);
fprintf('Régularisation latitude                  : %.3f deg\n', epsilon_lat_deg);
fprintf('Contrôle intégrale densité               : %.6f satellites\n', N_check);
fprintf('--------------------------------------------------------------------\n');
fprintf('p_link Delta utilisé                     : %.9f\n', p_link_delta);
fprintf('E[|E|] = N(N-1)p_link/2                  : %.6f\n', E_edges);
if ~isnan(E_edges_mat)
    fprintf('E_delta lu dans le .mat                  : %.6f\n', E_edges_mat);
end
fprintf('E[N1] utilisé                            : %.6f\n', E_N1);
fprintf('E[N2] utilisé                            : %.6f\n', E_N2);
fprintf('E[N3] utilisé                            : %.6f\n', E_N3);
fprintf('C_macro                                  : %.6f\n', C_macro);
fprintf('E[N1] local (comparaison)                : %.6f\n', E_N1_local);
fprintf('E[beta0] utilisé                         : %.6f\n', beta0_used);
fprintf('Source de E[beta0]                       : %s\n', beta0_source);
fprintf('--------------------------------------------------------------------\n');
fprintf('chi_bridge orbital                       : %.6f\n', chi_bridge);
fprintf('Nombre moyen de liens / composante       : %.6f\n', mean_links_per_component);
fprintf('Nombre moyen de ponts / composante       : %.6f\n', mean_bridges_per_component);
fprintf('====================================================================\n\n');

%% Graphiques
figure;
plot(rad2deg(phi), lambda_orb, 'LineWidth', 1.5);
grid on;
xlabel('Latitude \phi (deg)');
ylabel('\lambda_{orb}(\phi) (satellites/km^2)');
title('Densité orbitale locale régularisée');

figure;
plot(rad2deg(phi), chi_bridge_local, 'LineWidth', 1.5);
grid on;
xlabel('Latitude \phi (deg)');
ylabel('\chi_{bridge}(\phi)');
title('Fraction locale de liens sans voisin commun');

%% Sauvegarde
results = struct();
results.N = N;
results.R = R;
results.i_deg = i_deg;
results.dmax = dmax;
results.epsilon_lat_deg = epsilon_lat_deg;
results.phi = phi;
results.lambda_orb = lambda_orb;
results.chi_bridge_local = chi_bridge_local;
results.chi_bridge = chi_bridge;
results.p_link_delta = p_link_delta;
results.E_edges = E_edges;
results.E_edges_mat = E_edges_mat;
results.E_N1 = E_N1;
results.E_N2 = E_N2;
results.E_N3 = E_N3;
results.C_macro = C_macro;
results.E_N1_local = E_N1_local;
results.beta0_used = beta0_used;
results.mean_links_per_component = mean_links_per_component;
results.mean_bridges_per_component = mean_bridges_per_component;

save('chi_bridge_results.mat', 'results');

fprintf('Résultats sauvegardés dans chi_bridge_results.mat\n');
