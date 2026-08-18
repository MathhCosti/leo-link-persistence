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
%  7.c Vraie probabilite empirique qu'un nouveau lien relie
%      deux composantes differentes
%
% p_diffcomp_phi_emp(phi)
%   = P(C_i(t) ~= C_j(t) | nouveau lien entre t et t+dt, phi)
%
% Un nouveau lien est defini par :
%   A_ij(t) = 0 et A_ij(t+dt) = 1.
%
% On teste ensuite si ses deux extremites appartenaient a des
% composantes differentes dans G(t), AVANT apparition du lien.
%
% Pour etre coherent avec la theorie conditionnee par la latitude
% d'un satellite, chaque lien apporte un poids 1/2 dans la tranche
% de chacune de ses deux extremites a l'instant t.
%% ============================================================

[p_diffcomp_phi_emp, ...
 n_new_links_phi, ...
 n_merge_links_phi] = ...
    empirical_diffcomp_probability_phi( ...
        SA.Positions,SA.Adjacency,phi_edges,Nb);


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
%  8.b Diagnostic du flux brut de nouveaux liens
%
% On compare ici le nombre de nouveaux liens observe entre t et t+dt
% au nombre predit par le terme geometrique du modele, AVANT toute
% correction topologique.
%
% Pour un satellite localise en phi :
%
%   mu_new_link,sat(phi)
%      = 4*dmax*v_orb*dt*lambda(phi)*geometry_factor(phi)
%
% Le nombre theorique de nouveaux liens par tranche et par pas de temps
% est calcule par quadrature de f_Phi(phi)*mu_new_link(phi) dans chaque
% tranche, puis multiplie par N/2 pour ne compter chaque paire qu'une fois.
%% ============================================================

mu_new_link_per_sat_phi = nan(1,Nb);

valid_flux = ...
    valid_geometry ...
    & isfinite(lambda_phi) ...
    & lambda_phi >= 0;

mu_new_link_per_sat_phi(valid_flux) = ...
    4*dmax*v_orb*dt ...
    .* lambda_phi(valid_flux) ...
    .* geometry_factor(valid_flux);

% Masse theorique totale de latitude, pour controle.
phi_probability_mass = sum(phi_bin_probability,'omitnan');

if abs(phi_probability_mass-1) > 1e-8
    warning(['La somme de phi_bin_probability vaut %.12f au lieu de 1. ', ...
             'Le diagnostic utilise tout de meme la masse fournie.'], ...
             phi_probability_mass);
end

% Nombre moyen de satellites par tranche, conserve pour diagnostic.
N_per_bin_th = ...
    N_local .* phi_bin_probability;

% -----------------------------------------------------------------
% Quadrature dans chaque tranche
%
% Au lieu d'utiliser :
%
%   N_b * mu_new_link(phi_b)
%
% on calcule directement :
%
%   N_new,b^th
%     = N/2 * integral_b f_Phi(phi) mu_new_link(phi) dphi
%
% avec la loi exacte de latitude pour une phase orbitale uniforme :
%
%   f_Phi(phi)
%     = cos(phi)/(pi*sqrt(sin(i)^2-sin(phi)^2)),
%     |phi| < i.
%
% Cela evite l'artefact de discretisation pres de |phi| = i.
% -----------------------------------------------------------------

n_new_links_per_step_phi_th = zeros(1,Nb);

% Interpolation de lambda(phi) pour la quadrature.
lambda_interp = @(x) max( ...
    interp1(phi_density,lambda_density_source,x,'pchip','extrap'),0);

for b = 1:Nb

    a = max(phi_edges(b),-inc);
    c = min(phi_edges(b+1),inc);

    if c <= a
        n_new_links_per_step_phi_th(b) = 0;
        continue;
    end

    % On evite d'evaluer exactement aux bornes singulieres.
    eps_edge = 1e-10;
    aa = a;
    cc = c;

    if abs(aa + inc) < 1e-12
        aa = aa + eps_edge;
    end

    if abs(cc - inc) < 1e-12
        cc = cc - eps_edge;
    end

    integrand = @(phi) ...
        0.5 .* N_local ...
        .* latitude_pdf_delta(phi,inc) ...
        .* (4*dmax*v_orb*dt) ...
        .* lambda_interp(phi) ...
        .* geometry_factor_delta(phi,inc);

    n_new_links_per_step_phi_th(b) = ...
        integral(integrand,aa,cc, ...
                 'ArrayValued',true, ...
                 'RelTol',1e-7, ...
                 'AbsTol',1e-10);
end

% Nombre empirique moyen par transition.
Nt_transitions = max(min(numel(SA.Positions),numel(SA.Adjacency))-1,1);

n_new_links_per_step_phi_emp = ...
    n_new_links_phi ./ Nt_transitions;

% Totaux globaux moyens par transition.
n_new_links_per_step_th = ...
    sum(n_new_links_per_step_phi_th,'omitnan');

n_new_links_per_step_emp = ...
    sum(n_new_links_per_step_phi_emp,'omitnan');

% Facteur de surestimation du flux brut.
new_link_flux_ratio = ...
    n_new_links_per_step_th / max(n_new_links_per_step_emp,eps);

% Correction empirique locale du flux brut.
new_link_flux_correction_phi = nan(1,Nb);
valid_flux_ratio = ...
    isfinite(n_new_links_per_step_phi_th) ...
    & n_new_links_per_step_phi_th > 0 ...
    & isfinite(n_new_links_per_step_phi_emp);

new_link_flux_correction_phi(valid_flux_ratio) = ...
    n_new_links_per_step_phi_emp(valid_flux_ratio) ...
    ./ n_new_links_per_step_phi_th(valid_flux_ratio);


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
%  10.b Version avec beta0 empirique uniquement
%% ============================================================

merge_exponent_phi_betti_emp = nan(1,Nb);

valid_betti_emp = ...
    valid_topology_emp ...
    & valid_geometry ...
    & isfinite(lambda_phi) ...
    & isfinite(eta_sweep_phi);

% Meme modele que la theorie, mais N/beta0 est calcule avec
% beta0(phi) empirique. eta_sweep reste theorique.
merge_exponent_phi_betti_emp(valid_betti_emp) = ...
    4*dmax*v_orb*dt ...
    .* lambda_phi(valid_betti_emp) ...
    .* mean_satellites_per_component_phi_emp(valid_betti_emp) ...
    .* geometry_factor(valid_betti_emp) ...
    .* eta_sweep_phi(valid_betti_emp);

merge_exponent_phi_betti_emp = ...
    max(merge_exponent_phi_betti_emp,0);

p_merge_phi_th_betti_emp = nan(1,Nb);
p_merge_phi_th_betti_emp(valid_betti_emp) = ...
    1-exp(-merge_exponent_phi_betti_emp(valid_betti_emp));

p_merge_phi_th_betti_emp = ...
    min(max(p_merge_phi_th_betti_emp,0),1);

valid_integration_betti_emp = ...
    isfinite(p_merge_phi_th_betti_emp) ...
    & isfinite(betti0_per_bin_emp) ...
    & betti0_per_bin_emp >= 0;

component_mass_emp_used = ...
    sum(betti0_per_bin_emp(valid_integration_betti_emp));

if component_mass_emp_used <= 0
    error('La masse empirique totale de composantes est nulle.');
end

component_probability_bin_emp = zeros(1,Nb);
component_probability_bin_emp(valid_integration_betti_emp) = ...
    betti0_per_bin_emp(valid_integration_betti_emp) ...
    ./ component_mass_emp_used;

p_merge_th_betti_emp = ...
    sum(p_merge_phi_th_betti_emp(valid_integration_betti_emp) ...
    .* component_probability_bin_emp(valid_integration_betti_emp));

p_disp_fusion_th_betti_emp = ...
    0.5*p_merge_th_betti_emp;

%% ============================================================
%  10.c Version avec beta0 + eta_sweep empiriques
%% ============================================================

merge_exponent_phi_betti_eta_emp = nan(1,Nb);

valid_betti_eta_emp = ...
    valid_topology_emp ...
    & valid_geometry ...
    & isfinite(lambda_phi) ...
    & isfinite(eta_sweep_phi_emp);

merge_exponent_phi_betti_eta_emp(valid_betti_eta_emp) = ...
    4*dmax*v_orb*dt ...
    .* lambda_phi(valid_betti_eta_emp) ...
    .* mean_satellites_per_component_phi_emp(valid_betti_eta_emp) ...
    .* geometry_factor(valid_betti_eta_emp) ...
    .* eta_sweep_phi_emp(valid_betti_eta_emp);

merge_exponent_phi_betti_eta_emp = ...
    max(merge_exponent_phi_betti_eta_emp,0);

p_merge_phi_th_betti_eta_emp = nan(1,Nb);
p_merge_phi_th_betti_eta_emp(valid_betti_eta_emp) = ...
    1-exp(-merge_exponent_phi_betti_eta_emp(valid_betti_eta_emp));

p_merge_phi_th_betti_eta_emp = ...
    min(max(p_merge_phi_th_betti_eta_emp,0),1);

valid_integration_betti_eta_emp = ...
    isfinite(p_merge_phi_th_betti_eta_emp) ...
    & isfinite(betti0_per_bin_emp) ...
    & betti0_per_bin_emp >= 0;

component_mass_emp_eta_used = ...
    sum(betti0_per_bin_emp(valid_integration_betti_eta_emp));

if component_mass_emp_eta_used <= 0
    error('La masse empirique totale de composantes est nulle.');
end

component_probability_bin_emp_eta = zeros(1,Nb);
component_probability_bin_emp_eta(valid_integration_betti_eta_emp) = ...
    betti0_per_bin_emp(valid_integration_betti_eta_emp) ...
    ./ component_mass_emp_eta_used;

p_merge_th_betti_eta_emp = ...
    sum(p_merge_phi_th_betti_eta_emp(valid_integration_betti_eta_emp) ...
    .* component_probability_bin_emp_eta(valid_integration_betti_eta_emp));

p_disp_fusion_th_betti_eta_emp = ...
    0.5*p_merge_th_betti_eta_emp;


%% ============================================================
%  10.d Version avec beta0 empirique + vraie probabilite empirique
%
% eta_sweep(phi) est remplace par :
%
%   p_diffcomp_phi_emp(phi)
%     = P(C_i(t) ~= C_j(t) | nouveau lien, phi)
%
% Cette version conserve lambda(phi), N(phi), la vitesse relative
% et la structure du modele theorique, mais utilise les deux
% quantites topologiques empiriques beta0 et p_diffcomp.
%% ============================================================

merge_exponent_phi_betti_true_emp = nan(1,Nb);

valid_betti_true_emp = ...
    valid_topology_emp ...
    & valid_geometry ...
    & isfinite(lambda_phi) ...
    & isfinite(p_diffcomp_phi_emp);

merge_exponent_phi_betti_true_emp(valid_betti_true_emp) = ...
    4*dmax*v_orb*dt ...
    .* lambda_phi(valid_betti_true_emp) ...
    .* mean_satellites_per_component_phi_emp(valid_betti_true_emp) ...
    .* geometry_factor(valid_betti_true_emp) ...
    .* p_diffcomp_phi_emp(valid_betti_true_emp);

merge_exponent_phi_betti_true_emp = ...
    max(merge_exponent_phi_betti_true_emp,0);

p_merge_phi_th_betti_true_emp = nan(1,Nb);
p_merge_phi_th_betti_true_emp(valid_betti_true_emp) = ...
    1-exp(-merge_exponent_phi_betti_true_emp(valid_betti_true_emp));

p_merge_phi_th_betti_true_emp = ...
    min(max(p_merge_phi_th_betti_true_emp,0),1);

valid_integration_betti_true_emp = ...
    isfinite(p_merge_phi_th_betti_true_emp) ...
    & isfinite(betti0_per_bin_emp) ...
    & betti0_per_bin_emp >= 0;

component_mass_emp_true_used = ...
    sum(betti0_per_bin_emp(valid_integration_betti_true_emp));

if component_mass_emp_true_used <= 0
    error('La masse empirique totale de composantes est nulle.');
end

component_probability_bin_emp_true = zeros(1,Nb);
component_probability_bin_emp_true(valid_integration_betti_true_emp) = ...
    betti0_per_bin_emp(valid_integration_betti_true_emp) ...
    ./ component_mass_emp_true_used;

p_merge_th_betti_true_emp = ...
    sum(p_merge_phi_th_betti_true_emp(valid_integration_betti_true_emp) ...
    .* component_probability_bin_emp_true(valid_integration_betti_true_emp));

p_disp_fusion_th_betti_true_emp = ...
    0.5*p_merge_th_betti_true_emp;


%% ============================================================
%  10.e Version utilisant le VRAI flux empirique de nouveaux liens
%
% Le terme geometrique
%
%   4*dmax*v_orb*dt*lambda(phi)*g(phi)
%
% est remplace par le nombre reel de nouveaux liens observe entre
% t et t+dt dans chaque tranche.
%
% n_new_links_per_step_phi_emp compte chaque lien une seule fois :
% chaque extremite apporte un poids 1/2. Pour revenir au nombre
% d'incidences de nouveaux liens par composante, on multiplie donc
% par 2 avant de diviser par beta0 empirique :
%
%   mu_new,comp^emp(phi)
%      = 2*N_new^emp(phi) / beta0^emp(phi).
%
% Deux variantes sont calculees :
%   1) vrai flux de nouveaux liens + eta_sweep theorique ;
%   2) vrai flux de nouveaux liens + vraie P(C_i ~= C_j | nouveau lien).
%% ============================================================

mu_new_links_per_component_phi_emp = nan(1,Nb);

valid_newlinks_component = ...
    isfinite(n_new_links_per_step_phi_emp) ...
    & isfinite(betti0_per_bin_emp) ...
    & betti0_per_bin_emp > 0;

mu_new_links_per_component_phi_emp(valid_newlinks_component) = ...
    2 .* n_new_links_per_step_phi_emp(valid_newlinks_component) ...
    ./ betti0_per_bin_emp(valid_newlinks_component);

mu_new_links_per_component_phi_emp = ...
    max(mu_new_links_per_component_phi_emp,0);

% ------------------------------------------------------------
% 10.e.1 Vrai flux empirique + eta_sweep theorique
% ------------------------------------------------------------

merge_exponent_phi_newlinks_emp_eta_th = nan(1,Nb);

valid_newlinks_eta_th = ...
    valid_newlinks_component ...
    & isfinite(eta_sweep_phi);

merge_exponent_phi_newlinks_emp_eta_th(valid_newlinks_eta_th) = ...
    mu_new_links_per_component_phi_emp(valid_newlinks_eta_th) ...
    .* eta_sweep_phi(valid_newlinks_eta_th);

merge_exponent_phi_newlinks_emp_eta_th = ...
    max(merge_exponent_phi_newlinks_emp_eta_th,0);

p_merge_phi_newlinks_emp_eta_th = nan(1,Nb);
p_merge_phi_newlinks_emp_eta_th(valid_newlinks_eta_th) = ...
    1-exp(-merge_exponent_phi_newlinks_emp_eta_th(valid_newlinks_eta_th));

p_merge_phi_newlinks_emp_eta_th = ...
    min(max(p_merge_phi_newlinks_emp_eta_th,0),1);

valid_integration_newlinks_eta_th = ...
    isfinite(p_merge_phi_newlinks_emp_eta_th) ...
    & isfinite(betti0_per_bin_emp) ...
    & betti0_per_bin_emp >= 0;

mass_newlinks_eta_th = ...
    sum(betti0_per_bin_emp(valid_integration_newlinks_eta_th));

if mass_newlinks_eta_th <= 0
    error('Masse empirique nulle pour la version nouveaux liens + eta theorique.');
end

component_probability_bin_newlinks_eta_th = zeros(1,Nb);
component_probability_bin_newlinks_eta_th(valid_integration_newlinks_eta_th) = ...
    betti0_per_bin_emp(valid_integration_newlinks_eta_th) ...
    ./ mass_newlinks_eta_th;

p_merge_newlinks_emp_eta_th = ...
    sum(p_merge_phi_newlinks_emp_eta_th(valid_integration_newlinks_eta_th) ...
    .* component_probability_bin_newlinks_eta_th(valid_integration_newlinks_eta_th));

% ------------------------------------------------------------
% 10.e.2 Vrai flux empirique + vraie probabilite topologique
%
% C'est la version la plus directement corrigee :
%
%   mu_merge^emp(phi)
%      = [2*N_new^emp(phi)/beta0^emp(phi)]
%        * P(C_i ~= C_j | nouveau lien,phi).
% ------------------------------------------------------------

merge_exponent_phi_newlinks_true_emp = nan(1,Nb);

valid_newlinks_true_emp = ...
    valid_newlinks_component ...
    & isfinite(p_diffcomp_phi_emp);

merge_exponent_phi_newlinks_true_emp(valid_newlinks_true_emp) = ...
    mu_new_links_per_component_phi_emp(valid_newlinks_true_emp) ...
    .* p_diffcomp_phi_emp(valid_newlinks_true_emp);

merge_exponent_phi_newlinks_true_emp = ...
    max(merge_exponent_phi_newlinks_true_emp,0);

p_merge_phi_newlinks_true_emp = nan(1,Nb);
p_merge_phi_newlinks_true_emp(valid_newlinks_true_emp) = ...
    1-exp(-merge_exponent_phi_newlinks_true_emp(valid_newlinks_true_emp));

p_merge_phi_newlinks_true_emp = ...
    min(max(p_merge_phi_newlinks_true_emp,0),1);

valid_integration_newlinks_true_emp = ...
    isfinite(p_merge_phi_newlinks_true_emp) ...
    & isfinite(betti0_per_bin_emp) ...
    & betti0_per_bin_emp >= 0;

mass_newlinks_true_emp = ...
    sum(betti0_per_bin_emp(valid_integration_newlinks_true_emp));

if mass_newlinks_true_emp <= 0
    error('Masse empirique nulle pour la version vrais nouveaux liens.');
end

component_probability_bin_newlinks_true_emp = zeros(1,Nb);
component_probability_bin_newlinks_true_emp(valid_integration_newlinks_true_emp) = ...
    betti0_per_bin_emp(valid_integration_newlinks_true_emp) ...
    ./ mass_newlinks_true_emp;

p_merge_newlinks_true_emp = ...
    sum(p_merge_phi_newlinks_true_emp(valid_integration_newlinks_true_emp) ...
    .* component_probability_bin_newlinks_true_emp(valid_integration_newlinks_true_emp));

p_disp_fusion_newlinks_true_emp = ...
    0.5*p_merge_newlinks_true_emp;

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
plot(rad2deg(phi_vals),p_merge_phi_th_betti_emp, ...
    '--','LineWidth',2, ...
    'DisplayName','Beta_0 empirique');
plot(rad2deg(phi_vals),p_merge_phi_th_betti_eta_emp, ...
    ':','LineWidth',2, ...
    'DisplayName','Beta_0 + eta empiriques');
plot(rad2deg(phi_vals),p_merge_phi_th_betti_true_emp, ...
    '-.','LineWidth',2, ...
    'DisplayName','Beta_0 + P(C_i~=C_j|nouveau lien) empiriques');
plot(rad2deg(phi_vals),p_merge_phi_newlinks_emp_eta_th, ...
    '--','LineWidth',2, ...
    'DisplayName','Vrais nouveaux liens + \eta_{sweep}^{th}');
plot(rad2deg(phi_vals),p_merge_phi_newlinks_true_emp, ...
    '-','LineWidth',2.2, ...
    'DisplayName','Vrais nouveaux liens + vraie proba topologique');
grid on;
xlabel('Latitude \phi (deg)');
ylabel('p_{merge}^{\Delta}(\phi)');
title('Probabilite locale de fusion');
legend('Location','best');
ylim([0,max(1.05*max([p_merge_phi_th, ...
    p_merge_phi_th_betti_emp, ...
    p_merge_phi_th_betti_eta_emp, ...
    p_merge_phi_th_betti_true_emp, ...
    p_merge_phi_newlinks_emp_eta_th, ...
    p_merge_phi_newlinks_true_emp],[],'omitnan'),1e-3)]);
hold off;

figure;
plot(rad2deg(phi_vals),merge_exponent_phi,'LineWidth',1.8);
grid on;
xlabel('Latitude \phi (deg)');
ylabel('\Lambda_{merge}^{\Delta}(\phi)');
title('Exposant local du modele de fusion');


figure;
hold on;
plot(rad2deg(phi_vals),eta_sweep_phi, ...
    'LineWidth',1.8, ...
    'DisplayName','\eta_{sweep}^{th}');
plot(rad2deg(phi_vals),eta_sweep_phi_emp, ...
    '--','LineWidth',1.8, ...
    'DisplayName','\eta_{sweep}^{emp}');
plot(rad2deg(phi_vals),p_diffcomp_phi_emp, ...
    '-.','LineWidth',1.8, ...
    'DisplayName','P(C_i\neq C_j | nouveau lien)^{emp}');
grid on;
xlabel('Latitude \phi (deg)');
ylabel('Facteur de fusion');
title('Comparaison des facteurs correctifs locaux');
legend('Location','best');
ylim([0 1]);
hold off;



figure;
hold on;
plot(rad2deg(phi_vals),n_new_links_per_step_phi_th, ...
    'LineWidth',1.8, ...
    'DisplayName','Nouveaux liens theoriques');
plot(rad2deg(phi_vals),n_new_links_per_step_phi_emp, ...
    '--','LineWidth',1.8, ...
    'DisplayName','Nouveaux liens empiriques');
grid on;
xlabel('Latitude \phi (deg)');
ylabel('Nombre moyen de nouveaux liens / pas');
title('Diagnostic du flux brut de nouveaux liens');
legend('Location','best');
hold off;

figure;
plot(rad2deg(phi_vals),new_link_flux_correction_phi, ...
    'LineWidth',1.8);
grid on;
xlabel('Latitude \phi (deg)');
ylabel('Flux empirique / flux theorique');
title('Facteur correctif local du flux de nouveaux liens');

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
fprintf('p_merge avec beta0 empirique        : %.12f\n', ...
    p_merge_th_betti_emp);
fprintf('p_merge avec beta0 + eta empiriques : %.12f\n', ...
    p_merge_th_betti_eta_emp);
fprintf('p_merge beta0 + vraie proba emp     : %.12f\n', ...
    p_merge_th_betti_true_emp);
fprintf('p_merge vrais liens + eta th        : %.12f\n', ...
    p_merge_newlinks_emp_eta_th);
fprintf('p_merge vrais liens + vraie proba   : %.12f\n', ...
    p_merge_newlinks_true_emp);
fprintf('P(diff comp | nouveau lien) globale : %.12f\n', ...
    sum(n_merge_links_phi)/max(sum(n_new_links_phi),1));
fprintf('------------------------------------------------------------\n');
fprintf('Masse totale de phi_bin_probability : %.12f\n', ...
    phi_probability_mass);
fprintf('Nouveaux liens theo / pas           : %.12f\n', ...
    n_new_links_per_step_th);
fprintf('Nouveaux liens emp / pas            : %.12f\n', ...
    n_new_links_per_step_emp);
fprintf('Rapport flux theo / flux emp         : %.12f\n', ...
    new_link_flux_ratio);
fprintf('p_disp avec beta0 empirique         : %.12f\n', ...
    p_disp_fusion_th_betti_emp);
fprintf('p_disp beta0 + eta empiriques       : %.12f\n', ...
    p_disp_fusion_th_betti_eta_emp);
fprintf('p_disp beta0 + vraie proba emp      : %.12f\n', ...
    p_disp_fusion_th_betti_true_emp);
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
    'mu_new_link_per_sat_phi','N_per_bin_th', ...
    'phi_probability_mass', ...
    'n_new_links_per_step_phi_th','n_new_links_per_step_phi_emp', ...
    'n_new_links_per_step_th','n_new_links_per_step_emp', ...
    'new_link_flux_ratio','new_link_flux_correction_phi', ...
    'merge_exponent_phi','p_merge_phi_th', ...
    'component_mass_used','satellite_probability_mass_used', ...
    'p_merge_th','p_disp_fusion_th', ...
    'p_merge_th_satellite_weighting', ...
    'betti0_per_bin_emp','betti0_density_phi_emp', ...
    'component_probability_bin_emp','component_mass_emp_used', ...
    'mean_satellites_per_component_phi_emp', ...
    'eta_sweep_phi_emp', ...
    'p_diffcomp_phi_emp','n_new_links_phi','n_merge_links_phi', ...
    'merge_exponent_phi_betti_emp', ...
    'p_merge_phi_th_betti_emp', ...
    'p_merge_th_betti_emp', ...
    'p_disp_fusion_th_betti_emp', ...
    'merge_exponent_phi_betti_eta_emp', ...
    'p_merge_phi_th_betti_eta_emp', ...
    'p_merge_th_betti_eta_emp', ...
    'p_disp_fusion_th_betti_eta_emp', ...
    'merge_exponent_phi_betti_true_emp', ...
    'p_merge_phi_th_betti_true_emp', ...
    'p_merge_th_betti_true_emp', ...
    'p_disp_fusion_th_betti_true_emp', ...
    'component_probability_bin_emp_true', ...
    'component_mass_emp_true_used', ...
    'mu_new_links_per_component_phi_emp', ...
    'merge_exponent_phi_newlinks_emp_eta_th', ...
    'p_merge_phi_newlinks_emp_eta_th', ...
    'component_probability_bin_newlinks_eta_th', ...
    'p_merge_newlinks_emp_eta_th', ...
    'merge_exponent_phi_newlinks_true_emp', ...
    'p_merge_phi_newlinks_true_emp', ...
    'component_probability_bin_newlinks_true_emp', ...
    'p_merge_newlinks_true_emp', ...
    'p_disp_fusion_newlinks_true_emp', ...
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


function [p_diffcomp_phi_emp,n_new_links_phi,n_merge_links_phi] = ...
    empirical_diffcomp_probability_phi(Positions,Adjacency,phi_edges,Nb)

    Nt = min(numel(Positions),numel(Adjacency));

    n_new_links_phi = zeros(1,Nb);
    n_merge_links_phi = zeros(1,Nb);

    for t = 1:Nt-1

        pos = Positions{t};
        A0 = Adjacency{t};
        A1 = Adjacency{t+1};

        if isempty(pos) || isempty(A0) || isempty(A1)
            continue;
        end

        pos = double(pos);
        A0 = logical(A0);
        A1 = logical(A1);

        if size(pos,2) ~= 3
            error('Positions{%d} doit etre une matrice N x 3.',t);
        end

        if size(A0,1) ~= size(A1,1) || ...
                size(A0,2) ~= size(A1,2)
            error('Adjacency{%d} et Adjacency{%d} ont des tailles incompatibles.', ...
                t,t+1);
        end

        % Composantes du graphe AVANT apparition des nouveaux liens.
        G0 = graph(sparse(A0),'upper');
        component_id = conncomp(G0).';

        % Nouveaux liens entre t et t+dt.
        new_edges = triu(A1 & ~A0,1);
        [ii,jj] = find(new_edges);

        if isempty(ii)
            continue;
        end

        % Latitude des deux extremites a l'instant t.
        radius = sqrt(sum(pos.^2,2));
        phi_sat = asin(max(min(pos(:,3)./radius,1),-1));

        % On conditionne maintenant le comptage par la latitude
        % de CHAQUE extremite, comme dans le modele theorique.
        % Chaque nouveau lien apporte 1/2 dans la tranche de i
        % et 1/2 dans la tranche de j afin de compter une seule
        % fois le lien au total.
        bin_i = discretize(phi_sat(ii),phi_edges);
        bin_j = discretize(phi_sat(jj),phi_edges);

        % Un nouveau lien est topologiquement utile s'il connecte
        % deux composantes qui etaient distinctes avant sa creation.
        is_merge_link = ...
            component_id(ii) ~= component_id(jj);

        for k = 1:numel(ii)

            bi = bin_i(k);
            bj = bin_j(k);

            if ~isnan(bi)
                n_new_links_phi(bi) = ...
                    n_new_links_phi(bi) + 0.5;

                if is_merge_link(k)
                    n_merge_links_phi(bi) = ...
                        n_merge_links_phi(bi) + 0.5;
                end
            end

            if ~isnan(bj)
                n_new_links_phi(bj) = ...
                    n_new_links_phi(bj) + 0.5;

                if is_merge_link(k)
                    n_merge_links_phi(bj) = ...
                        n_merge_links_phi(bj) + 0.5;
                end
            end
        end
    end

    p_diffcomp_phi_emp = nan(1,Nb);

    valid = n_new_links_phi > 0;

    p_diffcomp_phi_emp(valid) = ...
        n_merge_links_phi(valid) ...
        ./ n_new_links_phi(valid);
end


function f = latitude_pdf_delta(phi,inc)
    % Loi de latitude pour une phase orbitale uniforme.
    %
    % f_Phi(phi) = cos(phi)/(pi*sqrt(sin(i)^2-sin(phi)^2))
    % pour |phi| < i.

    rad = sin(inc).^2 - sin(phi).^2;

    f = zeros(size(phi));

    valid = ...
        abs(phi) < inc ...
        & rad > 0;

    f(valid) = ...
        cos(phi(valid)) ...
        ./ (pi*sqrt(rad(valid)));
end

function g = geometry_factor_delta(phi,inc)
    % Facteur local de vitesse relative :
    %
    % g(phi) = sqrt(sin(i)^2-sin(phi)^2)/cos(phi).

    rad = sin(inc).^2 - sin(phi).^2;

    g = zeros(size(phi));

    valid = ...
        abs(phi) < inc ...
        & rad >= 0 ...
        & abs(cos(phi)) > 1e-12;

    g(valid) = ...
        sqrt(max(rad(valid),0)) ...
        ./ cos(phi(valid));
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
