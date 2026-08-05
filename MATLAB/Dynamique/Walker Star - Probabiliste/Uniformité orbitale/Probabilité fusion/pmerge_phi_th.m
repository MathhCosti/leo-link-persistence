%% pmerge_phi_th.m
% Probabilite theorique locale de fusion p_merge^Delta(phi)
% Walker Delta a uniformite orbitale.
%
% Formule locale :
%
% p_merge(phi) = 1 - exp( -4*dmax*v_orb*dt*lambda(phi)
%                         * [N^Delta(phi)/B_0^Delta(phi)]
%                         * sqrt(sin(i)^2-sin(phi)^2)/cos(phi)
%                         * eta_sweep(phi) )
%
% avec :
%   N^Delta(phi)       : densite de satellites par radian de latitude ;
%   B_0^Delta(phi)     : densite locale de composantes par radian ;
%   lambda(phi)        : densite surfacique locale ;
%   eta_sweep(phi)     : facteur local de redondance spatiale.
%
% La probabilite globale est une moyenne par composante :
%
%   f_{Phi_C}(phi) =
%       B_0^Delta(phi) / integral B_0^Delta(psi) dpsi
%
% puis :
%
%   p_merge =
%       integral p_merge(phi) B_0^Delta(phi) dphi
%       ------------------------------------------------
%              integral B_0^Delta(phi) dphi.
%
% Numeriquement, chaque tranche est ponderee par le nombre theorique
% de composantes fourni par betti_phi_results.mat.
%
% Entrees :
%   N_phi_results.mat
%   betti_phi_results.mat
%   densite_phi_results.mat
%   eta_sweep_phi_results.mat
%
% Sortie :
%   pmerge_phi_th_results.mat

clear; clc; close all;

%% ============================================================
%  1. Parametres
%% ============================================================

% Utilise uniquement si aucun fichier analysis_temp_results.mat
% contenant dt n'est trouve.
dt_default = 10; % s

% Choix du facteur eta_sweep :
%   'theory'             : eta obtenu avec lambda theorique ;
%   'empirical_density'  : eta obtenu avec lambda empirique ;
%   'direct'             : mesure directe de eta_sweep.
eta_source = 'theory';

% Choix de la densite surfacique :
%   'theory'    : lambda_theory_bins ;
%   'empirical' : lambda_empirical.
lambda_source = 'theory';

%% ============================================================
%  2. Chargement
%% ============================================================

script_dir = fileparts(mfilename('fullpath'));
N_file = fullfile(script_dir, '..', '..', '..', '..', 'Statique', 'Uniformité orbitale', 'Valeurs locales', 'N_phi_results.mat');

betti_file = fullfile(script_dir, '..', '..', '..', '..', 'Statique', 'Uniformité orbitale', 'Betti', 'betti_phi_results.mat');

density_file = fullfile(script_dir, '..', 'Paramètres', 'densite_phi_results.mat');

eta_file = fullfile(script_dir, '..', 'Paramètres', 'eta_sweep_phi_results.mat');

SN = load(N_file);
SB = load(betti_file);
SD = load(density_file);
SE = load(eta_file);

require_fields(SN,{ ...
    'N','phi_vals','phi_bin_probability', ...
    'satellites_density_phi_th'},N_file);

require_fields(SB,{ ...
    'N','R','inc','dmax','phi_vals', ...
    'betti0_density_phi_th'},betti_file);

require_fields(SD,{ ...
    'N','phi_centers','lambda_theory_bins','lambda_empirical'}, ...
    density_file);

require_fields(SE,{ ...
    'phi_vals','eta_sweep_phi_th'},eta_file);

%% ============================================================
%  3. Grille locale de reference
%
% On utilise la grille commune de N(phi) et beta_0(phi).
%% ============================================================

phi_vals = row_vector(SN.phi_vals);
phi_bin_probability = row_vector(SN.phi_bin_probability);

phi_betti = row_vector(SB.phi_vals);

if numel(phi_vals) ~= numel(phi_betti) || ...
        max(abs(phi_vals-phi_betti)) > 1e-10
    error(['N_phi_results.mat et betti_phi_results.mat doivent ', ...
           'utiliser la meme grille de latitude.']);
end

Nb = numel(phi_vals);

satellites_density_phi = ...
    row_vector(SN.satellites_density_phi_th);

betti0_density_phi = ...
    row_vector(SB.betti0_density_phi_th);

if numel(satellites_density_phi) ~= Nb || ...
        numel(betti0_density_phi) ~= Nb
    error('Dimensions incoherentes pour N(phi) ou beta_0(phi).');
end

% Nombre theorique de composantes dans chaque tranche.
if isfield(SB,'betti0_per_bin_th')
    betti0_per_bin_th = row_vector(SB.betti0_per_bin_th);
elseif isfield(SB,'betti0_bin_th')
    betti0_per_bin_th = row_vector(SB.betti0_bin_th);
else
    phi_edges_betti = zeros(1,Nb+1);
    phi_edges_betti(1) = -inc;
    phi_edges_betti(end) = inc;
    phi_edges_betti(2:end-1) = ...
        0.5*(phi_vals(1:end-1)+phi_vals(2:end));

    dphi_betti = diff(phi_edges_betti);
    betti0_per_bin_th = betti0_density_phi .* dphi_betti;
end

if numel(betti0_per_bin_th) ~= Nb
    error(['Le nombre de valeurs beta0 par tranche est incoherent ', ...
           'avec la grille phi_vals.']);
end

betti0_per_bin_th = max(betti0_per_bin_th,0);

N_local = double(SN.N);
N_betti = double(SB.N);

R = double(SB.R);
inc = double(SB.inc);
dmax = double(SB.dmax);

%% ============================================================
%  4. Recuperation de dt et de la vitesse orbitale
%% ============================================================

analysis_candidates = {
    fullfile(script_dir,'analysis_temp_results.mat')
    fullfile(script_dir,'analysis_temp_results(8).mat')
    fullfile(script_dir,'analysis_temp_results(7).mat')
    fullfile(script_dir,'..','analysis_temp_results.mat')
    fullfile(script_dir,'..','analysis_temp_results(8).mat')
    fullfile(script_dir,'..','analysis_temp_results(7).mat')
};

analysis_file = '';
dt = dt_default;
mu = 398600; % km^3/s^2

for k = 1:numel(analysis_candidates)
    if isfile(analysis_candidates{k})
        analysis_file = analysis_candidates{k};
        SA = load(analysis_file);

        if isfield(SA,'dt')
            dt = double(SA.dt);
        end

        if isfield(SA,'mu')
            mu = double(SA.mu);
        end

        break;
    end
end

v_orb = sqrt(mu/R);

%% ============================================================
%  5. Densite surfacique lambda(phi)
%
% La grille du fichier de densite peut differer. On interpole les
% valeurs moyennees par bande sur la grille locale de reference.
%% ============================================================

phi_density = row_vector(SD.phi_centers);

switch lower(lambda_source)
    case 'theory'
        lambda_density_source = ...
            row_vector(SD.lambda_theory_bins);

    case 'empirical'
        lambda_density_source = ...
            row_vector(SD.lambda_empirical);

    otherwise
        error('lambda_source doit valoir ''theory'' ou ''empirical''.');
end

lambda_phi = interp1( ...
    phi_density,lambda_density_source,phi_vals, ...
    'pchip','extrap');

% Protection contre une extrapolation numerique negative
lambda_phi = max(lambda_phi,0);

%% ============================================================
%  6. Facteur eta_sweep(phi)
%% ============================================================

phi_eta = row_vector(SE.phi_vals);

switch lower(eta_source)
    case 'theory'
        eta_source_values = ...
            row_vector(SE.eta_sweep_phi_th);

    case 'empirical_density'
        if ~isfield(SE,'eta_sweep_phi_from_empirical_density')
            error(['Le fichier eta_sweep ne contient pas ', ...
                   'eta_sweep_phi_from_empirical_density.']);
        end

        eta_source_values = ...
            row_vector(SE.eta_sweep_phi_from_empirical_density);

    case 'direct'
        if isfield(SE,'eta_sweep_phi_emp_direct')
            eta_source_values = ...
                row_vector(SE.eta_sweep_phi_emp_direct);
        elseif isfield(SE,'eta_sweep_phi_emp')
            eta_source_values = ...
                row_vector(SE.eta_sweep_phi_emp);
        else
            error(['Le fichier eta_sweep ne contient aucune ', ...
                   'mesure directe reconnue.']);
        end

    otherwise
        error(['eta_source doit valoir ''theory'', ', ...
               '''empirical_density'' ou ''direct''.']);
end

eta_sweep_phi = interp1( ...
    phi_eta,eta_source_values,phi_vals, ...
    'pchip','extrap');

eta_sweep_phi = min(max(eta_sweep_phi,0),1);

% Mesure empirique directe de eta_sweep, utilisee uniquement pour
% la version corrigee empiriquement.
if isfield(SE,'eta_sweep_phi_emp_direct')
    eta_sweep_phi_emp_source = ...
        row_vector(SE.eta_sweep_phi_emp_direct);
elseif isfield(SE,'eta_sweep_phi_emp')
    eta_sweep_phi_emp_source = ...
        row_vector(SE.eta_sweep_phi_emp);
else
    error(['Le fichier eta_sweep_phi_results.mat doit contenir ', ...
           'eta_sweep_phi_emp_direct ou eta_sweep_phi_emp pour ', ...
           'calculer la version corrigee empiriquement.']);
end

eta_sweep_phi_emp = interp1( ...
    phi_eta,eta_sweep_phi_emp_source,phi_vals, ...
    'pchip','extrap');

eta_sweep_phi_emp = min(max(eta_sweep_phi_emp,0),1);

%% ============================================================
%  7. Mise a l'echelle d'une composante
%
% Nbar_C(phi) = N^Delta(phi) / B_0^Delta(phi)
%% ============================================================

mean_satellites_per_component_phi = nan(1,Nb);

valid_topology = ...
    isfinite(satellites_density_phi) ...
    & isfinite(betti0_density_phi) ...
    & betti0_density_phi > 0;

mean_satellites_per_component_phi(valid_topology) = ...
    satellites_density_phi(valid_topology) ...
    ./ betti0_density_phi(valid_topology);

%% ============================================================
%  7.b Betti local empirique issu de analysis_temp_results.mat
%
% Chaque composante de taille s apporte un poids 1/s a chacun de
% ses satellites. La somme de ses contributions vaut donc 1.
%% ============================================================

if isempty(analysis_file)
    error(['analysis_temp_results.mat est necessaire pour calculer ', ...
           'les Betti locaux empiriques.']);
end

if ~isfield(SA,'Positions') || ~isfield(SA,'Adjacency')
    error(['Le fichier %s doit contenir Positions et Adjacency ', ...
           'pour reconstruire beta0(phi) empiriquement.'], ...
           analysis_file);
end

phi_edges = zeros(1,Nb+1);
phi_edges(1) = -inc;
phi_edges(end) = inc;
phi_edges(2:end-1) = ...
    0.5*(phi_vals(1:end-1)+phi_vals(2:end));

dphi_bins = diff(phi_edges);

betti0_per_bin_emp = ...
    empirical_betti_per_bin( ...
        SA.Positions,SA.Adjacency,phi_edges,Nb);

betti0_density_phi_emp = ...
    betti0_per_bin_emp ./ dphi_bins;

valid_topology_emp = ...
    isfinite(satellites_density_phi) ...
    & isfinite(betti0_density_phi_emp) ...
    & betti0_density_phi_emp > 0;

mean_satellites_per_component_phi_emp = nan(1,Nb);

mean_satellites_per_component_phi_emp(valid_topology_emp) = ...
    satellites_density_phi(valid_topology_emp) ...
    ./ betti0_density_phi_emp(valid_topology_emp);

%% ============================================================
%  8. Vitesse relative locale
%
% v_rel(phi) =
%   2*v_orb*sqrt(sin(i)^2-sin(phi)^2)/cos(phi)
%% ============================================================

geometry_factor = nan(1,Nb);

radicand = sin(inc)^2-sin(phi_vals).^2;

valid_geometry = ...
    abs(phi_vals) < inc ...
    & radicand >= 0 ...
    & abs(cos(phi_vals)) > 1e-12;

geometry_factor(valid_geometry) = ...
    sqrt(max(radicand(valid_geometry),0)) ...
    ./ cos(phi_vals(valid_geometry));

v_rel_phi = 2*v_orb*geometry_factor;

%% ============================================================
%  9. Probabilite locale de fusion
%% ============================================================

merge_exponent_phi = nan(1,Nb);

valid = ...
    valid_topology ...
    & valid_geometry ...
    & isfinite(lambda_phi) ...
    & isfinite(eta_sweep_phi);

merge_exponent_phi(valid) = ...
    4*dmax*v_orb*dt ...
    .* lambda_phi(valid) ...
    .* mean_satellites_per_component_phi(valid) ...
    .* geometry_factor(valid) ...
    .* eta_sweep_phi(valid);

merge_exponent_phi = max(merge_exponent_phi,0);

p_merge_phi_th = nan(1,Nb);
p_merge_phi_th(valid) = ...
    1-exp(-merge_exponent_phi(valid));

p_merge_phi_th = min(max(p_merge_phi_th,0),1);

%% ============================================================
%  10. Integration globale par composante
%% ============================================================

valid_integration = ...
    isfinite(p_merge_phi_th) ...
    & isfinite(betti0_per_bin_th) ...
    & betti0_per_bin_th >= 0;

component_mass_used = ...
    sum(betti0_per_bin_th(valid_integration));

if component_mass_used <= 0
    error('La masse totale de composantes issue du fichier Betti est nulle.');
end

component_probability_bin = zeros(1,Nb);
component_probability_bin(valid_integration) = ...
    betti0_per_bin_th(valid_integration) ...
    ./ component_mass_used;

p_merge_th = ...
    sum(p_merge_phi_th(valid_integration) ...
    .* component_probability_bin(valid_integration));

% Ancienne moyenne par satellite, gardee seulement pour comparaison.
valid_satellite_integration = ...
    isfinite(p_merge_phi_th) ...
    & isfinite(phi_bin_probability);

satellite_probability_mass_used = ...
    sum(phi_bin_probability(valid_satellite_integration));

p_merge_th_satellite_weighting = ...
    sum(p_merge_phi_th(valid_satellite_integration) ...
    .* phi_bin_probability(valid_satellite_integration));

if satellite_probability_mass_used > 0 && ...
        abs(satellite_probability_mass_used-1) > 1e-10
    p_merge_th_satellite_weighting = ...
        p_merge_th_satellite_weighting ...
        / satellite_probability_mass_used;
end

% Disparition d'une barre H0 lors d'une fusion binaire.
p_disp_fusion_th = 0.5*p_merge_th;

%% ============================================================
%  10.b Version corrigee avec Betti et eta_sweep empiriques
%% ============================================================

merge_exponent_phi_corrected = nan(1,Nb);

valid_corrected = ...
    valid_topology_emp ...
    & valid_geometry ...
    & isfinite(lambda_phi) ...
    & isfinite(eta_sweep_phi_emp);

merge_exponent_phi_corrected(valid_corrected) = ...
    4*dmax*v_orb*dt ...
    .* lambda_phi(valid_corrected) ...
    .* mean_satellites_per_component_phi_emp(valid_corrected) ...
    .* geometry_factor(valid_corrected) ...
    .* eta_sweep_phi_emp(valid_corrected);

merge_exponent_phi_corrected = ...
    max(merge_exponent_phi_corrected,0);

p_merge_phi_th_corrected = nan(1,Nb);
p_merge_phi_th_corrected(valid_corrected) = ...
    1-exp(-merge_exponent_phi_corrected(valid_corrected));

p_merge_phi_th_corrected = ...
    min(max(p_merge_phi_th_corrected,0),1);

valid_integration_corrected = ...
    isfinite(p_merge_phi_th_corrected) ...
    & isfinite(betti0_per_bin_emp) ...
    & betti0_per_bin_emp >= 0;

component_mass_emp_used = ...
    sum(betti0_per_bin_emp(valid_integration_corrected));

if component_mass_emp_used <= 0
    error('La masse empirique totale de composantes est nulle.');
end

component_probability_bin_emp = zeros(1,Nb);
component_probability_bin_emp(valid_integration_corrected) = ...
    betti0_per_bin_emp(valid_integration_corrected) ...
    ./ component_mass_emp_used;

p_merge_th_corrected = ...
    sum(p_merge_phi_th_corrected(valid_integration_corrected) ...
    .* component_probability_bin_emp(valid_integration_corrected));

p_disp_fusion_th_corrected = ...
    0.5*p_merge_th_corrected;

%% Facteurs correctifs locaux, donnes a titre de diagnostic
correction_betti_phi = nan(1,Nb);
valid_betti_correction = ...
    betti0_density_phi > 0 ...
    & betti0_density_phi_emp > 0;

correction_betti_phi(valid_betti_correction) = ...
    betti0_density_phi(valid_betti_correction) ...
    ./ betti0_density_phi_emp(valid_betti_correction);

correction_eta_phi = nan(1,Nb);
valid_eta_correction = eta_sweep_phi > 0;

correction_eta_phi(valid_eta_correction) = ...
    eta_sweep_phi_emp(valid_eta_correction) ...
    ./ eta_sweep_phi(valid_eta_correction);

%% ============================================================
%  11. Graphe
%% ============================================================

figure;
hold on;
plot(rad2deg(phi_vals),p_merge_phi_th, ...
    'LineWidth',2, ...
    'DisplayName','Theorie');
plot(rad2deg(phi_vals),p_merge_phi_th_corrected, ...
    '--','LineWidth',2, ...
    'DisplayName','Correction Betti + eta empiriques');
grid on;
xlabel('Latitude \phi (deg)');
ylabel('p_{merge}^{\Delta}(\phi)');
title('Probabilite locale de fusion');
legend('Location','best');
ylim([0,max(1.05*max([p_merge_phi_th, ...
    p_merge_phi_th_corrected],[],'omitnan'),1e-3)]);
hold off;

figure;
plot(rad2deg(phi_vals),merge_exponent_phi,'LineWidth',1.8);
grid on;
xlabel('Latitude \phi (deg)');
ylabel('\Lambda_{merge}^{\Delta}(\phi)');
title('Exposant local du modele de fusion');


figure;
hold on;
plot(rad2deg(phi_vals),component_probability_bin, ...
    'LineWidth',1.8, ...
    'DisplayName','Betti theorique');
plot(rad2deg(phi_vals),component_probability_bin_emp, ...
    '--','LineWidth',1.8, ...
    'DisplayName','Betti empirique');
grid on;
xlabel('Latitude \phi (deg)');
ylabel('P(\Phi_C dans la tranche)');
title('Loi de latitude des composantes');
legend('Location','best');
hold off;

%% ============================================================
%  12. Affichage console
%% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' p_merge(phi) THEORIQUE - WALKER DELTA ORBITAL\n');
fprintf('============================================================\n');
fprintf('N utilise dans N(phi)               : %.8f\n',N_local);
fprintf('N utilise dans beta0(phi)           : %.8f\n',N_betti);
fprintf('N utilise dans lambda(phi)          : %.8f\n',double(SD.N));
fprintf('Rayon orbital R                     : %.8f km\n',R);
fprintf('Inclinaison                         : %.8f deg\n',rad2deg(inc));
fprintf('dmax                                : %.8f km\n',dmax);
fprintf('dt                                  : %.8f s\n',dt);
fprintf('v_orb                               : %.8f km/s\n',v_orb);
fprintf('Source de lambda                    : %s\n',lambda_source);
fprintf('Source de eta_sweep                 : %s\n',eta_source);
fprintf('Nombre theorique de composantes     : %.12f\n', ...
    component_mass_used);
fprintf('Somme des poids de composantes      : %.12f\n', ...
    sum(component_probability_bin));
fprintf('------------------------------------------------------------\n');
fprintf('p_merge pondere par composantes     : %.12f\n',p_merge_th);
fprintf('p_disp fusion theorique             : %.12f\n', ...
    p_disp_fusion_th);
fprintf('p_merge pondere par satellites      : %.12f\n', ...
    p_merge_th_satellite_weighting);
fprintf('------------------------------------------------------------\n');
fprintf('Nombre empirique moyen composantes  : %.12f\n', ...
    component_mass_emp_used);
fprintf('p_merge corrige Betti + eta emp     : %.12f\n', ...
    p_merge_th_corrected);
fprintf('p_disp fusion corrige               : %.12f\n', ...
    p_disp_fusion_th_corrected);
fprintf('============================================================\n');

if abs(N_local-N_betti) > 1e-8 || ...
        abs(N_local-double(SD.N)) > 1e-8
    warning(['Les fichiers locaux n''utilisent pas tous le meme N. ', ...
             'Regenerer densite_phi, N_phi, betti_phi et eta_sweep_phi ', ...
             'avec le N de analysis_temp_results.mat pour une ', ...
             'comparaison strictement coherente.']);
end

%% ============================================================
%  13. Sauvegarde
%% ============================================================

output_file = fullfile(script_dir,'pmerge_phi_th_results.mat');

save(output_file, ...
    'N_local','N_betti','R','inc','dmax','dt','mu','v_orb', ...
    'phi_vals','phi_bin_probability', ...
    'satellites_density_phi','betti0_density_phi', ...
    'betti0_per_bin_th','component_probability_bin', ...
    'mean_satellites_per_component_phi', ...
    'lambda_phi','lambda_source', ...
    'eta_sweep_phi','eta_source', ...
    'geometry_factor','v_rel_phi', ...
    'merge_exponent_phi','p_merge_phi_th', ...
    'component_mass_used','satellite_probability_mass_used', ...
    'p_merge_th','p_disp_fusion_th', ...
    'p_merge_th_satellite_weighting', ...
    'betti0_per_bin_emp','betti0_density_phi_emp', ...
    'component_probability_bin_emp','component_mass_emp_used', ...
    'mean_satellites_per_component_phi_emp', ...
    'eta_sweep_phi_emp', ...
    'merge_exponent_phi_corrected', ...
    'p_merge_phi_th_corrected', ...
    'p_merge_th_corrected', ...
    'p_disp_fusion_th_corrected', ...
    'correction_betti_phi','correction_eta_phi', ...
    'N_file','betti_file','density_file','eta_file', ...
    'analysis_file');

fprintf('Resultats sauvegardes dans %s\n',output_file);

%% ============================================================
%  Fonctions locales
%% ============================================================

function path_out = find_result_file(script_dir,file_names)
    search_dirs = {
        script_dir
        fullfile(script_dir,'..')
        fullfile(script_dir,'..','Valeurs locales')
    };

    path_out = '';

    for d = 1:numel(search_dirs)
        for k = 1:numel(file_names)
            candidate = fullfile(search_dirs{d},file_names{k});

            if isfile(candidate)
                path_out = candidate;
                return;
            end
        end
    end

    error('Aucun fichier trouve parmi : %s.', ...
        strjoin(file_names,', '));
end

function betti0_per_bin_emp = empirical_betti_per_bin( ...
        Positions,Adjacency,phi_edges,Nb)

    Nt = min(numel(Positions),numel(Adjacency));
    betti_sum = zeros(1,Nb);
    n_valid_times = 0;

    for t = 1:Nt
        pos = Positions{t};
        A = Adjacency{t};

        if isempty(pos) || isempty(A)
            continue;
        end

        pos = double(pos);

        if size(pos,2) ~= 3
            error('Positions{%d} doit etre une matrice N x 3.',t);
        end

        % Latitude geocentrique de chaque satellite.
        radius = sqrt(sum(pos.^2,2));
        phi_sat = asin(max(min(pos(:,3)./radius,1),-1));

        % Identifiants des composantes connexes.
        G = graph(sparse(A),'upper');
        component_id = conncomp(G).';

        bin_id = discretize(phi_sat,phi_edges);

        local_mass = zeros(1,Nb);
        component_labels = unique(component_id);

        for c = component_labels.'
            nodes = find(component_id == c);
            s = numel(nodes);

            valid_nodes = nodes(~isnan(bin_id(nodes)));

            if s > 0 && ~isempty(valid_nodes)
                local_mass = local_mass + ...
                    accumarray( ...
                        bin_id(valid_nodes), ...
                        ones(numel(valid_nodes),1)/s, ...
                        [Nb,1],@sum,0).';
            end
        end

        betti_sum = betti_sum + local_mass;
        n_valid_times = n_valid_times + 1;
    end

    if n_valid_times == 0
        error('Aucun instant valide pour reconstruire beta0(phi).');
    end

    betti0_per_bin_emp = betti_sum/n_valid_times;
end

function require_fields(S,names,file_name)
    for k = 1:numel(names)
        if ~isfield(S,names{k})
            error('Le fichier %s doit contenir %s.', ...
                file_name,names{k});
        end
    end
end

function x = row_vector(x)
    x = double(x(:).');
end
