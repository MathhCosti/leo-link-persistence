%% pbreak_phi_th.m
% Probabilite theorique locale de rupture p_break^Delta(phi)
% Walker Delta a uniformite orbitale.
%
% Formule locale :
%
%   p_break^Delta(phi)
%     = [B0(phi)-N1(phi)]/B0(phi)
%       * {1-exp[-mu_break^non-isole(phi)]}
%
% avec
%
%   mu_break^non-isole(phi)
%     = E(phi)/[B0(phi)-N1(phi)]
%       * p_break^lien(phi)
%       * p_bridge,bord(phi),
%
%   p_break^lien(phi)
%     = (2/pi) v_rel^Delta(phi) dt/dmax,
%
%   v_rel^Delta(phi)
%     = 2 v_orb sqrt(sin(i)^2-sin(phi)^2)/cos(phi),
%
% et, pour la version theorique,
%
%   p_bridge,bord(phi)
%     ~= exp[-lambda(phi) A_inter(dmax)].
%
% Deux versions corrigees sont comparees :
%   - beta0 empirique + eta_sweep empirique utilise comme p_bridge ;
%   - beta0 empirique + vraie P(pont | rupture,phi).
%
% La probabilite globale est une moyenne par composante :
%
%   p_break^Delta
%     = sum_b p_break^Delta(phi_b) B0_b
%       --------------------------------
%                  sum_b B0_b.
%
% Entrees :
%   edges_phi_results.mat
%   N1_phi_results.mat
%   betti_phi_results.mat
%   densite_phi_results.mat
%   analysis_temp_results.mat
%
% Sortie :
%   pbreak_phi_th_results.mat

clear; clc; close all;

%% 1. Parametres
dt_default = 10;
mu_default = 398600;
lambda_source = 'theory';

%% 2. Chargement
script_dir = fileparts(mfilename('fullpath'));

betti_file = fullfile(script_dir, '..', '..', '..', '..', 'Statique', 'Uniformité orbitale', 'Betti', 'betti_phi_results.mat');

density_file = fullfile(script_dir, '..', 'Paramètres', 'densite_phi_results.mat');

eta_file_candidates = {
    fullfile(script_dir, '..', 'Paramètres', 'eta_sweep_phi_results.mat')
    fullfile(script_dir, '..', 'Paramètres', 'eta_sweep_phi_results(2).mat')
    fullfile(script_dir, 'eta_sweep_phi_results.mat')
    fullfile(script_dir, 'eta_sweep_phi_results(2).mat')
};

eta_file = '';
for k = 1:numel(eta_file_candidates)
    if isfile(eta_file_candidates{k})
        eta_file = eta_file_candidates{k};
        break;
    end
end

if isempty(eta_file)
    error('Fichier eta_sweep_phi_results.mat introuvable.');
end

edges_file = fullfile(script_dir, '..', '..', '..', '..', 'Statique', 'Uniformité orbitale', 'Valeurs locales', 'edges_phi_results.mat');

N1_file = fullfile(script_dir, '..', '..', '..', '..', 'Statique', 'Uniformité orbitale', 'Betti', 'N1_phi_results.mat');

SE = load(edges_file);
S1 = load(N1_file);
SB = load(betti_file);
SD = load(density_file);
Seta = load(eta_file);

require_fields(SE,{'phi_vals','edges_density_phi_th'},edges_file);
require_fields(S1,{'phi_vals','N1_density_phi_th'},N1_file);
require_fields(SB,{'R','inc','dmax','phi_vals','betti0_density_phi_th'},betti_file);
require_fields(SD,{'phi_centers','lambda_theory_bins','lambda_empirical'},density_file);
require_fields(Seta,{'phi_vals', ...
    'eta_sweep_phi_emp_direct', ...
    'p_bridge_bord_phi_emp_true'},eta_file);

%% 3. Grille de reference et quantites topologiques
phi_vals = row_vector(SB.phi_vals);
Nb = numel(phi_vals);
R = double(SB.R);
inc = double(SB.inc);
dmax = double(SB.dmax);

betti0_density_phi = row_vector(SB.betti0_density_phi_th);

if isfield(SB,'betti0_per_bin_th')
    betti0_per_bin_th = row_vector(SB.betti0_per_bin_th);
elseif isfield(SB,'betti0_bin_th')
    betti0_per_bin_th = row_vector(SB.betti0_bin_th);
else
    phi_edges = build_edges(phi_vals,-inc,inc);
    dphi_bins = diff(phi_edges);
    betti0_per_bin_th = betti0_density_phi .* dphi_bins;
end

betti0_per_bin_th = max(betti0_per_bin_th,0);

phi_edges_source = row_vector(SE.phi_vals);
edges_density_source = row_vector(SE.edges_density_phi_th);
edges_density_phi = interp1(phi_edges_source,edges_density_source,phi_vals,'pchip','extrap');
edges_density_phi = max(edges_density_phi,0);

phi_N1 = row_vector(S1.phi_vals);
N1_density_source = row_vector(S1.N1_density_phi_th);
N1_density_phi = interp1(phi_N1,N1_density_source,phi_vals,'pchip','extrap');
N1_density_phi = max(N1_density_phi,0);
N1_density_phi = min(N1_density_phi,betti0_density_phi);

nonisolated_density_phi = max(betti0_density_phi-N1_density_phi,0);

%% 4. Recuperation de dt et vitesse orbitale
analysis_candidates = {
    fullfile(script_dir,'..','analysis_temp_results.mat')
};

analysis_file = '';
dt = dt_default;
mu = mu_default;

for k = 1:numel(analysis_candidates)
    if isfile(analysis_candidates{k})
        analysis_file = analysis_candidates{k};
        SA = load(analysis_file);
        if isfield(SA,'dt'), dt = double(SA.dt); end
        if isfield(SA,'mu')
            mu = double(SA.mu);
        elseif isfield(SA,'mu_earth')
            mu = double(SA.mu_earth);
        end
        break;
    end
end

v_orb = sqrt(mu/R);

if isempty(analysis_file)
    error(['analysis_temp_results.mat est necessaire pour reconstruire ', ...
           'beta0(phi) empiriquement.']);
end

if ~isfield(SA,'Positions') || ~isfield(SA,'Adjacency')
    error(['Le fichier %s doit contenir Positions et Adjacency ', ...
           'pour reconstruire beta0(phi) empiriquement.'],analysis_file);
end

%% 4.b Reconstruction empirique de beta0(phi)
%
% Chaque composante de taille s contribue pour 1/s dans la tranche
% de chacun de ses satellites. Sa contribution totale vaut ainsi 1.
%% ============================================================

if isfield(SB,'phi_edges')
    phi_edges_betti = row_vector(SB.phi_edges);
else
    phi_edges_betti = build_edges(phi_vals,-inc,inc);
end

dphi_betti = diff(phi_edges_betti);

betti0_per_bin_emp = empirical_betti_per_bin( ...
    SA.Positions,SA.Adjacency,phi_edges_betti,Nb);

betti0_density_phi_emp = ...
    betti0_per_bin_emp ./ dphi_betti;

% N1 reste theorique : la correction demandee porte ici sur beta0.
N1_density_phi_emp_corrected = ...
    min(N1_density_phi,betti0_density_phi_emp);

nonisolated_density_phi_emp = ...
    max(betti0_density_phi_emp-N1_density_phi_emp_corrected,0);

%% 5. Densite surfacique locale
phi_density = row_vector(SD.phi_centers);
switch lower(lambda_source)
    case 'theory'
        lambda_density_source = row_vector(SD.lambda_theory_bins);
    case 'empirical'
        lambda_density_source = row_vector(SD.lambda_empirical);
    otherwise
        error('lambda_source doit valoir ''theory'' ou ''empirical''.');
end

lambda_phi = interp1(phi_density,lambda_density_source,phi_vals,'pchip','extrap');
lambda_phi = max(lambda_phi,0);

%% 6. Vitesse relative locale et rupture d'un lien
geometry_factor = nan(1,Nb);
radicand = sin(inc)^2-sin(phi_vals).^2;
valid_geometry = abs(phi_vals) < inc & radicand >= 0 & abs(cos(phi_vals)) > 1e-12;
geometry_factor(valid_geometry) = sqrt(max(radicand(valid_geometry),0)) ./ cos(phi_vals(valid_geometry));

v_rel_phi = nan(1,Nb);
v_rel_phi(valid_geometry) = 2*v_orb.*geometry_factor(valid_geometry);

p_break_link_phi = nan(1,Nb);
p_break_link_phi(valid_geometry) = (2/pi).*v_rel_phi(valid_geometry).*dt/dmax;
p_break_link_phi = min(max(p_break_link_phi,0),1-eps);

%% 7. Probabilite locale qu'un lien de bord soit un pont
A_inter_at_dmax = (2*pi/3-sqrt(3)/2)*dmax^2;
p_bridge_bord_phi = exp(-lambda_phi.*A_inter_at_dmax);
p_bridge_bord_phi = min(max(p_bridge_bord_phi,0),1);

% Deux corrections empiriques sont comparees :
%
% 1) approximation empirique issue de eta_sweep :
%
%    p_bridge,bord^eta(phi) := eta_sweep^emp(phi),
%
%    comme dans le code corrections_empiriques ;
%
% 2) vraie probabilite empirique :
%
%    p_bridge,bord^vrai(phi)
%      = P(pont dans G_t | lien rompu entre t et t+dt,phi).
%
phi_emp = row_vector(Seta.phi_vals);

%% 7.a Approximation empirique par eta_sweep
p_bridge_bord_phi_eta_source = ...
    row_vector(Seta.eta_sweep_phi_emp_direct);

valid_eta_source = ...
    isfinite(phi_emp) ...
    & isfinite(p_bridge_bord_phi_eta_source);

if nnz(valid_eta_source) < 2
    error(['Pas assez de tranches valides dans ', ...
           'eta_sweep_phi_emp_direct pour interpoler.']);
end

p_bridge_bord_phi_emp_eta = interp1( ...
    phi_emp(valid_eta_source), ...
    p_bridge_bord_phi_eta_source(valid_eta_source), ...
    phi_vals,'pchip','extrap');

p_bridge_bord_phi_emp_eta = ...
    min(max(p_bridge_bord_phi_emp_eta,0),1);

%% 7.b Vraie probabilite empirique de pont
p_bridge_bord_phi_true_source = ...
    row_vector(Seta.p_bridge_bord_phi_emp_true);

valid_true_source = ...
    isfinite(phi_emp) ...
    & isfinite(p_bridge_bord_phi_true_source);

if nnz(valid_true_source) < 2
    error(['Pas assez de tranches valides dans ', ...
           'p_bridge_bord_phi_emp_true pour interpoler.']);
end

p_bridge_bord_phi_emp_true = interp1( ...
    phi_emp(valid_true_source), ...
    p_bridge_bord_phi_true_source(valid_true_source), ...
    phi_vals,'linear','extrap');

p_bridge_bord_phi_emp_true = ...
    min(max(p_bridge_bord_phi_emp_true,0),1);

%% Facteurs correctifs locaux
correction_p_bridge_eta_phi = nan(1,Nb);
correction_p_bridge_true_phi = nan(1,Nb);

valid_p_bridge_correction = p_bridge_bord_phi > 0;

correction_p_bridge_eta_phi(valid_p_bridge_correction) = ...
    p_bridge_bord_phi_emp_eta(valid_p_bridge_correction) ...
    ./ p_bridge_bord_phi(valid_p_bridge_correction);

correction_p_bridge_true_phi(valid_p_bridge_correction) = ...
    p_bridge_bord_phi_emp_true(valid_p_bridge_correction) ...
    ./ p_bridge_bord_phi(valid_p_bridge_correction);

%% 8. Passage a l'echelle d'une composante non isolee
mean_links_per_nonisolated_component_phi = nan(1,Nb);
valid_nonisolated = isfinite(edges_density_phi) & isfinite(nonisolated_density_phi) & nonisolated_density_phi > 0;
mean_links_per_nonisolated_component_phi(valid_nonisolated) = edges_density_phi(valid_nonisolated) ./ nonisolated_density_phi(valid_nonisolated);

mu_break_nonisolated_phi = nan(1,Nb);
valid_local = valid_nonisolated & valid_geometry & isfinite(p_break_link_phi) & isfinite(p_bridge_bord_phi);
mu_break_nonisolated_phi(valid_local) = mean_links_per_nonisolated_component_phi(valid_local) .* p_break_link_phi(valid_local) .* p_bridge_bord_phi(valid_local);
mu_break_nonisolated_phi = max(mu_break_nonisolated_phi,0);

p_break_nonisolated_phi = nan(1,Nb);
p_break_nonisolated_phi(valid_local) = 1-exp(-mu_break_nonisolated_phi(valid_local));
p_break_nonisolated_phi = min(max(p_break_nonisolated_phi,0),1);

%% 9. Deconditionnement local
fraction_nonisolated_phi = zeros(1,Nb);
valid_betti = isfinite(betti0_density_phi) & betti0_density_phi > 0;
fraction_nonisolated_phi(valid_betti) = nonisolated_density_phi(valid_betti) ./ betti0_density_phi(valid_betti);
fraction_nonisolated_phi = min(max(fraction_nonisolated_phi,0),1);

p_break_phi_th = nan(1,Nb);
p_break_phi_th(valid_local) = fraction_nonisolated_phi(valid_local) .* p_break_nonisolated_phi(valid_local);
p_break_phi_th = min(max(p_break_phi_th,0),1);

p_break_phi_th_linear = nan(1,Nb);
p_break_phi_th_linear(valid_local) = fraction_nonisolated_phi(valid_local) .* mu_break_nonisolated_phi(valid_local);
p_break_phi_th_linear = min(max(p_break_phi_th_linear,0),1);

p_break_phi_th_linear_simplified = nan(1,Nb);
valid_linear_simplified = valid_betti & valid_geometry & isfinite(edges_density_phi) & isfinite(p_break_link_phi) & isfinite(p_bridge_bord_phi);
p_break_phi_th_linear_simplified(valid_linear_simplified) = edges_density_phi(valid_linear_simplified) ./ betti0_density_phi(valid_linear_simplified) .* p_break_link_phi(valid_linear_simplified) .* p_bridge_bord_phi(valid_linear_simplified);
p_break_phi_th_linear_simplified = min(max(p_break_phi_th_linear_simplified,0),1);

%% 9.b Versions corrigees avec beta0 empirique
%
% La densite d'aretes et N1 restent theoriques. On compare :
%
%   - correction eta :
%       B0_th(phi) -> B0_emp(phi)
%       p_bridge_th(phi) -> eta_sweep_emp(phi)
%
%   - correction vraie :
%       B0_th(phi) -> B0_emp(phi)
%       p_bridge_th(phi) -> P(pont | rupture,phi).
%% ============================================================

mean_links_per_nonisolated_component_phi_corrected = nan(1,Nb);

valid_nonisolated_corrected = ...
    isfinite(edges_density_phi) ...
    & isfinite(nonisolated_density_phi_emp) ...
    & nonisolated_density_phi_emp > 0;

mean_links_per_nonisolated_component_phi_corrected( ...
        valid_nonisolated_corrected) = ...
    edges_density_phi(valid_nonisolated_corrected) ...
    ./ nonisolated_density_phi_emp(valid_nonisolated_corrected);

fraction_nonisolated_phi_corrected = zeros(1,Nb);

valid_betti_emp = ...
    isfinite(betti0_density_phi_emp) ...
    & betti0_density_phi_emp > 0;

fraction_nonisolated_phi_corrected(valid_betti_emp) = ...
    nonisolated_density_phi_emp(valid_betti_emp) ...
    ./ betti0_density_phi_emp(valid_betti_emp);

fraction_nonisolated_phi_corrected = ...
    min(max(fraction_nonisolated_phi_corrected,0),1);

% ---------- Correction avec eta_sweep empirique ----------
[mu_break_nonisolated_phi_corrected_eta, ...
 p_break_nonisolated_phi_corrected_eta, ...
 p_break_phi_th_corrected_eta, ...
 p_break_phi_th_linear_corrected_eta, ...
 valid_local_corrected_eta] = corrected_local_model( ...
    mean_links_per_nonisolated_component_phi_corrected, ...
    fraction_nonisolated_phi_corrected, ...
    p_break_link_phi, ...
    p_bridge_bord_phi_emp_eta, ...
    valid_nonisolated_corrected, ...
    valid_geometry);

% ---------- Correction avec la vraie probabilite de pont ----------
[mu_break_nonisolated_phi_corrected_true, ...
 p_break_nonisolated_phi_corrected_true, ...
 p_break_phi_th_corrected_true, ...
 p_break_phi_th_linear_corrected_true, ...
 valid_local_corrected_true] = corrected_local_model( ...
    mean_links_per_nonisolated_component_phi_corrected, ...
    fraction_nonisolated_phi_corrected, ...
    p_break_link_phi, ...
    p_bridge_bord_phi_emp_true, ...
    valid_nonisolated_corrected, ...
    valid_geometry);

% Facteur local lie a beta0 dans l'exposant conditionnel.
correction_betti_phi = nan(1,Nb);
valid_betti_correction = ...
    nonisolated_density_phi > 0 ...
    & nonisolated_density_phi_emp > 0;

correction_betti_phi(valid_betti_correction) = ...
    nonisolated_density_phi(valid_betti_correction) ...
    ./ nonisolated_density_phi_emp(valid_betti_correction);

%% 10. Integration globale par composante
valid_integration = isfinite(p_break_phi_th) & isfinite(betti0_per_bin_th) & betti0_per_bin_th >= 0;
component_mass_used = sum(betti0_per_bin_th(valid_integration));
if component_mass_used <= 0
    error('La masse theorique totale de composantes est nulle.');
end

component_probability_bin = zeros(1,Nb);
component_probability_bin(valid_integration) = betti0_per_bin_th(valid_integration) ./ component_mass_used;

p_break_th = sum(p_break_phi_th(valid_integration) .* component_probability_bin(valid_integration));

valid_linear_integration = isfinite(p_break_phi_th_linear) & isfinite(betti0_per_bin_th) & betti0_per_bin_th >= 0;
linear_mass = sum(betti0_per_bin_th(valid_linear_integration));
p_break_th_linear = sum(p_break_phi_th_linear(valid_linear_integration) .* betti0_per_bin_th(valid_linear_integration)) / linear_mass;
p_break_th_linear_simplified = sum(p_break_phi_th_linear_simplified(valid_linear_integration) .* betti0_per_bin_th(valid_linear_integration)) / linear_mass;

% Integration corrigee selon la loi empirique des composantes.
valid_integration_corrected_eta = ...
    isfinite(p_break_phi_th_corrected_eta) ...
    & isfinite(betti0_per_bin_emp) ...
    & betti0_per_bin_emp >= 0;

valid_integration_corrected_true = ...
    isfinite(p_break_phi_th_corrected_true) ...
    & isfinite(betti0_per_bin_emp) ...
    & betti0_per_bin_emp >= 0;

valid_integration_corrected = ...
    valid_integration_corrected_eta ...
    & valid_integration_corrected_true;

component_mass_emp_used = ...
    sum(betti0_per_bin_emp(valid_integration_corrected));

if component_mass_emp_used <= 0
    error('La masse empirique totale de composantes est nulle.');
end

component_probability_bin_emp = zeros(1,Nb);
component_probability_bin_emp(valid_integration_corrected) = ...
    betti0_per_bin_emp(valid_integration_corrected) ...
    ./ component_mass_emp_used;

p_break_th_corrected_eta = ...
    sum(p_break_phi_th_corrected_eta(valid_integration_corrected) ...
    .* component_probability_bin_emp(valid_integration_corrected));

p_break_th_linear_corrected_eta = ...
    sum(p_break_phi_th_linear_corrected_eta(valid_integration_corrected) ...
    .* component_probability_bin_emp(valid_integration_corrected));

p_break_th_corrected_true = ...
    sum(p_break_phi_th_corrected_true(valid_integration_corrected) ...
    .* component_probability_bin_emp(valid_integration_corrected));

p_break_th_linear_corrected_true = ...
    sum(p_break_phi_th_linear_corrected_true(valid_integration_corrected) ...
    .* component_probability_bin_emp(valid_integration_corrected));

% Alias conserves pour compatibilite avec pbreak_temp.m :
% la version "corrected" designe maintenant la correction vraie.
p_break_phi_th_corrected = p_break_phi_th_corrected_true;
p_break_phi_th_linear_corrected = ...
    p_break_phi_th_linear_corrected_true;
p_break_th_corrected = p_break_th_corrected_true;
p_break_th_linear_corrected = ...
    p_break_th_linear_corrected_true;

%% 11. Figures
figure;
hold on;
plot(rad2deg(phi_vals),p_break_phi_th,'LineWidth',2,'DisplayName','Forme probabiliste');
plot(rad2deg(phi_vals),p_break_phi_th_linear,'--','LineWidth',1.8,'DisplayName','Approximation lineaire');
plot(rad2deg(phi_vals),p_break_phi_th_corrected_eta, ...
    ':','LineWidth',2, ...
    'DisplayName','Correction \beta_0 + \eta_{sweep}^{emp}');
plot(rad2deg(phi_vals),p_break_phi_th_corrected_true, ...
    '-.','LineWidth',2, ...
    'DisplayName','Correction \beta_0 + vraie P(pont|rupture)');
grid on;
xlabel('Latitude \phi (deg)');
ylabel('p_{break}^{\Delta}(\phi)');
title('Probabilite theorique locale de rupture');
legend('Location','best');
hold off;

figure;
plot(rad2deg(phi_vals),p_break_link_phi,'LineWidth',2);
grid on;
xlabel('Latitude \phi (deg)');
ylabel('p_{break}^{lien,\Delta}(\phi)');
title('Probabilite locale de rupture d''un lien');

figure;
hold on;
plot(rad2deg(phi_vals),p_bridge_bord_phi,'LineWidth',2,'DisplayName','p_{bridge,bord}^{\Delta}(\phi)');
plot(rad2deg(phi_vals),p_bridge_bord_phi_emp_eta, ...
    ':','LineWidth',2, ...
    'DisplayName','Approximation \eta_{sweep}^{emp}(\phi)');
plot(rad2deg(phi_vals),p_bridge_bord_phi_emp_true, ...
    '-.','LineWidth',2, ...
    'DisplayName','Vraie P(pont | rupture,\phi)');
plot(rad2deg(phi_vals),fraction_nonisolated_phi,'--','LineWidth',1.8,'DisplayName','Fraction non isolee theorique');
plot(rad2deg(phi_vals),fraction_nonisolated_phi_corrected, ...
    '-.','LineWidth',1.8, ...
    'DisplayName','Fraction non isolee corrigee');
grid on;
xlabel('Latitude \phi (deg)');
ylabel('Probabilite');
title('Facteurs topologiques locaux');
legend('Location','best');
hold off;

figure;
hold on;
plot(rad2deg(phi_vals),component_probability_bin, ...
    'LineWidth',1.8,'DisplayName','Betti theorique');
plot(rad2deg(phi_vals),component_probability_bin_emp, ...
    '--','LineWidth',1.8,'DisplayName','Betti empirique');
grid on;
xlabel('Latitude \phi (deg)');
ylabel('P(\Phi_C dans la tranche)');
title('Loi de latitude des composantes');
legend('Location','best');
hold off;

%% 12. Affichage console
fprintf('\n');
fprintf('============================================================\n');
fprintf(' p_break(phi) THEORIQUE - WALKER DELTA ORBITAL\n');
fprintf('============================================================\n');
fprintf('Rayon orbital R                     : %.8f km\n',R);
fprintf('Inclinaison                         : %.8f deg\n',rad2deg(inc));
fprintf('dmax                                : %.8f km\n',dmax);
fprintf('dt                                  : %.8f s\n',dt);
fprintf('v_orb                               : %.8f km/s\n',v_orb);
fprintf('Source de lambda                    : %s\n',lambda_source);
fprintf('Nombre theorique de composantes     : %.12f\n',component_mass_used);
fprintf('Somme des poids de composantes      : %.12f\n',sum(component_probability_bin));
fprintf('------------------------------------------------------------\n');
fprintf('p_break global probabiliste         : %.12f\n',p_break_th);
fprintf('p_break global lineaire             : %.12f\n',p_break_th_linear);
fprintf('p_break lineaire simplifie          : %.12f\n',p_break_th_linear_simplified);
fprintf('------------------------------------------------------------\n');
fprintf('Nombre empirique moyen composantes  : %.12f\n',component_mass_emp_used);

if isfield(Seta,'p_bridge_bord_emp_true_global')
    fprintf('Vraie P(pont | rupture) globale     : %.12f\n', ...
        double(Seta.p_bridge_bord_emp_true_global));
end

fprintf('p_break corrige beta0 + eta emp     : %.12f\n', ...
    p_break_th_corrected_eta);
fprintf('p_break corrige beta0 + vrai pbridge: %.12f\n', ...
    p_break_th_corrected_true);
fprintf('p_break lineaire corrige eta        : %.12f\n', ...
    p_break_th_linear_corrected_eta);
fprintf('p_break lineaire corrige vrai pont  : %.12f\n', ...
    p_break_th_linear_corrected_true);
fprintf('============================================================\n');

%% 13. Sauvegarde
output_file = fullfile(script_dir,'pbreak_phi_th_results.mat');
save(output_file, ...
    'R','inc','dmax','dt','mu','v_orb', ...
    'phi_vals', ...
    'edges_density_phi', ...
    'N1_density_phi', ...
    'betti0_density_phi','betti0_per_bin_th', ...
    'nonisolated_density_phi', ...
    'component_probability_bin','component_mass_used', ...
    'lambda_phi','lambda_source', ...
    'geometry_factor','v_rel_phi', ...
    'p_break_link_phi', ...
    'A_inter_at_dmax','p_bridge_bord_phi', ...
    'mean_links_per_nonisolated_component_phi', ...
    'mu_break_nonisolated_phi', ...
    'p_break_nonisolated_phi', ...
    'fraction_nonisolated_phi', ...
    'p_break_phi_th', ...
    'p_break_phi_th_linear', ...
    'p_break_phi_th_linear_simplified', ...
    'p_break_th', ...
    'p_break_th_linear', ...
    'p_break_th_linear_simplified', ...
    'betti0_per_bin_emp','betti0_density_phi_emp', ...
    'nonisolated_density_phi_emp', ...
    'component_probability_bin_emp','component_mass_emp_used', ...
    'p_bridge_bord_phi_eta_source','valid_eta_source', ...
    'p_bridge_bord_phi_true_source','valid_true_source', ...
    'p_bridge_bord_phi_emp_eta','p_bridge_bord_phi_emp_true', ...
    'correction_p_bridge_eta_phi', ...
    'correction_p_bridge_true_phi', ...
    'mean_links_per_nonisolated_component_phi_corrected', ...
    'fraction_nonisolated_phi_corrected', ...
    'mu_break_nonisolated_phi_corrected_eta', ...
    'p_break_nonisolated_phi_corrected_eta', ...
    'p_break_phi_th_corrected_eta', ...
    'p_break_phi_th_linear_corrected_eta', ...
    'valid_local_corrected_eta', ...
    'mu_break_nonisolated_phi_corrected_true', ...
    'p_break_nonisolated_phi_corrected_true', ...
    'p_break_phi_th_corrected_true', ...
    'p_break_phi_th_linear_corrected_true', ...
    'valid_local_corrected_true', ...
    'correction_betti_phi', ...
    'p_break_th_corrected_eta', ...
    'p_break_th_linear_corrected_eta', ...
    'p_break_th_corrected_true', ...
    'p_break_th_linear_corrected_true', ...
    'p_break_phi_th_corrected', ...
    'p_break_phi_th_linear_corrected', ...
    'p_break_th_corrected', ...
    'p_break_th_linear_corrected', ...
    'edges_file','N1_file','betti_file','density_file', ...
    'eta_file','analysis_file');

fprintf('Resultats sauvegardes dans %s\n',output_file);

%% Fonctions locales
function path_out = find_result_file(script_dir,file_names)
    search_dirs = {
        script_dir
        fullfile(script_dir,'..')
        fullfile(script_dir,'..','Valeurs locales')
        fullfile(script_dir,'..','Paramètres')
        fullfile(script_dir,'..','Betti')
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
    error('Aucun fichier trouve parmi : %s.',strjoin(file_names,', '));
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

        radius = sqrt(sum(pos.^2,2));
        phi_sat = asin(max(min(pos(:,3)./radius,1),-1));

        G = graph(sparse(A),'upper');
        component_id = conncomp(G).';

        bin_id = discretize(phi_sat,phi_edges);
        local_mass = zeros(1,Nb);

        for c = unique(component_id).'
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

function [mu_phi,p_nonisolated_phi,p_global_phi, ...
        p_linear_phi,valid_local] = corrected_local_model( ...
        mean_links_phi,fraction_nonisolated_phi, ...
        p_break_link_phi,p_bridge_phi, ...
        valid_nonisolated,valid_geometry)

    valid_local = ...
        valid_nonisolated ...
        & valid_geometry ...
        & isfinite(p_break_link_phi) ...
        & isfinite(p_bridge_phi);

    mu_phi = nan(size(mean_links_phi));
    mu_phi(valid_local) = ...
        mean_links_phi(valid_local) ...
        .* p_break_link_phi(valid_local) ...
        .* p_bridge_phi(valid_local);

    mu_phi = max(mu_phi,0);

    p_nonisolated_phi = nan(size(mean_links_phi));
    p_nonisolated_phi(valid_local) = ...
        1-exp(-mu_phi(valid_local));

    p_nonisolated_phi = ...
        min(max(p_nonisolated_phi,0),1);

    p_global_phi = nan(size(mean_links_phi));
    p_global_phi(valid_local) = ...
        fraction_nonisolated_phi(valid_local) ...
        .* p_nonisolated_phi(valid_local);

    p_global_phi = min(max(p_global_phi,0),1);

    p_linear_phi = nan(size(mean_links_phi));
    p_linear_phi(valid_local) = ...
        fraction_nonisolated_phi(valid_local) ...
        .* mu_phi(valid_local);

    p_linear_phi = min(max(p_linear_phi,0),1);
end

function require_fields(S,names,file_name)
    for k = 1:numel(names)
        if ~isfield(S,names{k})
            error('Le fichier %s doit contenir %s.',file_name,names{k});
        end
    end
end

function x = row_vector(x)
    x = double(x(:).');
end

function edges = build_edges(centers,left_edge,right_edge)
    centers = row_vector(centers);
    edges = zeros(1,numel(centers)+1);
    edges(1) = left_edge;
    edges(end) = right_edge;
    edges(2:end-1) = 0.5*(centers(1:end-1)+centers(2:end));
end
