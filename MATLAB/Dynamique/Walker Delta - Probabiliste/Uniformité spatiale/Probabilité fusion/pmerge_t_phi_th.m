%% pmerge_t_phi_script.m
% Calcul de p_merge(t,phi) a partir des resultats locaux deja produits.
%
% Formule locale :
%
%   p_merge(t,phi)
%     = 1-exp[-2*dmax*DeltaT*lambda(t,phi)
%              * N_C(t,phi)*v_rel(t,phi)*eta_sweep(t,phi)]
%
% avec, sur les tranches de latitude,
%
%   N_C(t,b) = N_b(t)/beta_0(t,b).
%
% La moyenne globale est ponderee par le nombre local de composantes :
%
%   p_merge_global(t)
%     = sum_b p_merge(t,b) beta_0(t,b)
%       / sum_b beta_0(t,b).
%
% Le script calcule aussi :
%
%   p_disp_fusion(t,phi) = p_merge(t,phi)/2,
%
% sous l'hypothese de fusions binaires.
%
% Fichiers requis :
%   N_t_phi_results.mat
%   lambda_t_phi_results.mat
%   eta_sweep_t_phi_results.mat
%   betti_t_phi_results.mat
%   vrel_t_phi_results.mat
%
% Sortie :
%   pmerge_t_phi_results.mat

clear; clc; close all;

%% ============================================================
%  1. Parametres du script
%% ============================================================

% Laisser vide pour utiliser median(diff(time_values)).
DeltaT = [];

make_figures = true;
save_results = true;

% Dossier contenant les fichiers de resultats.
script_dir = fileparts(mfilename('fullpath'));

if isempty(script_dir)
    script_dir = pwd;
end

search_dir = script_dir;

%% ============================================================
%  2. Recherche et chargement des fichiers
%% ============================================================

file_N = fullfile(script_dir, '..', 'Paramètres', 'N_t_phi_results.mat');
file_lambda = fullfile(script_dir, '..', 'Paramètres', 'lambda_t_phi_results.mat');
file_eta = fullfile(script_dir, '..', 'Paramètres', 'eta_sweep_t_phi_results.mat');
file_betti = fullfile(script_dir, '..', 'Paramètres', 'Betti', 'betti_t_phi_results.mat');
file_vrel = fullfile(script_dir, '..', 'Paramètres', 'Vitesse relative', 'vrel_t_phi_results.mat');

% Nouvelle theorie de beta0(t,phi) par expansion PPP en tailles de
% composantes, produite par beta0_t_phi_th_delta_spatial.m.
beta0_new_candidates = { ...
    fullfile(script_dir,'beta0_t_phi_results.mat'), ...
    fullfile(script_dir,'..','beta0_t_phi_results.mat'), ...
    fullfile(script_dir,'..','Paramètres','Betti','beta0_t_phi_results.mat'), ...
    fullfile(script_dir,'..','..','Paramètres','Betti','beta0_t_phi_results.mat')};

file_beta0_new = '';
for k = 1:numel(beta0_new_candidates)
    if isfile(beta0_new_candidates{k})
        file_beta0_new = beta0_new_candidates{k};
        break;
    end
end

if isempty(file_beta0_new)
    error(['beta0_t_phi_results.mat introuvable. ', ...
           'Executer auparavant le code beta0(t,phi) PPP.']);
end

SN = load(file_N);
SL = load(file_lambda);
SE = load(file_eta);
SB = load(file_betti);
SV = load(file_vrel);
SBnew = load(file_beta0_new);

check_fields(SN,{ ...
    'time_values','phi_vals_emp','phi_edges_emp', ...
    'satellite_count_th'},file_N);

check_fields(SL,{ ...
    'time_values','phi_vals_emp','lambda_bin_th'},file_lambda);

check_fields(SE,{ ...
    'time_values','phi_vals_emp','eta_sweep_bin_th','eta_sweep_emp'},file_eta);

check_fields(SB,{ ...
    'time_values','phi_vals','beta0_th','beta0_emp_true'},file_betti);

check_fields(SV,{ ...
    'time_values','phi_vals_emp','v_rel_th_on_emp'},file_vrel);

check_fields(SBnew,{ ...
    'time_values','phi_vals','beta0_bin_th'},file_beta0_new);

%% ============================================================
%  3. Grille commune
%% ============================================================

time_values = double(SN.time_values(:));
phi_vals = double(SN.phi_vals_emp(:).');
phi_edges = double(SN.phi_edges_emp(:).');

Nt = numel(time_values);
Nb = numel(phi_vals);

assert_same_grid(time_values,phi_vals, ...
    double(SL.time_values(:)),double(SL.phi_vals_emp(:).'),file_lambda);

assert_same_grid(time_values,phi_vals, ...
    double(SE.time_values(:)),double(SE.phi_vals_emp(:).'),file_eta);

assert_same_grid(time_values,phi_vals, ...
    double(SB.time_values(:)),double(SB.phi_vals(:).'),file_betti);

assert_same_grid(time_values,phi_vals, ...
    double(SV.time_values(:)),double(SV.phi_vals_emp(:).'),file_vrel);

assert_same_grid(time_values,phi_vals, ...
    double(SBnew.time_values(:)),double(SBnew.phi_vals(:).'),file_beta0_new);

if isempty(DeltaT)
    if Nt < 2
        error('Impossible de deduire DeltaT avec un seul instant.');
    end

    DeltaT = median(diff(time_values));
end

%% ============================================================
%  4. Recuperation des grandeurs locales
%% ============================================================

N_sat_th = double(SN.satellite_count_th);
lambda_th = double(SL.lambda_bin_th);
eta_sweep_th = double(SE.eta_sweep_bin_th);
eta_sweep_emp = double(SE.eta_sweep_emp);
beta0_th = double(SB.beta0_th);
beta0_emp = double(SB.beta0_emp_true);

% Nouvelle theorie PPP de beta0(t,phi), deja exprimee en nombre moyen
% de composantes PAR TRANCHE, donc directement compatible avec beta0_th.
beta0_ppp = double(SBnew.beta0_bin_th);

v_rel_th = double(SV.v_rel_th_on_emp);

check_size(N_sat_th,Nt,Nb,'satellite_count_th');
check_size(lambda_th,Nt,Nb,'lambda_bin_th');
check_size(eta_sweep_th,Nt,Nb,'eta_sweep_bin_th');
check_size(eta_sweep_emp,Nt,Nb,'eta_sweep_emp');
check_size(beta0_th,Nt,Nb,'beta0_th');
check_size(beta0_emp,Nt,Nb,'beta0_emp_true');
check_size(beta0_ppp,Nt,Nb,'beta0_bin_th (nouveau beta0 PPP)');
check_size(v_rel_th,Nt,Nb,'v_rel_th_on_emp');

% dmax est prioritairement lu depuis lambda/eta/vrel/betti.
dmax = get_first_scalar_field( ...
    {SL,SE,SV,SB,SN},'dmax');

%% ============================================================
%  4.b Donnees empiriques pour le flux de nouveaux liens
%% ============================================================

analysis_candidates = { ...
    fullfile(script_dir,'analysis_temp_results.mat'), ...
    fullfile(script_dir,'..','analysis_temp_results.mat'), ...
    fullfile(script_dir,'..','..','analysis_temp_results.mat')};

analysis_file = '';
for k = 1:numel(analysis_candidates)
    if isfile(analysis_candidates{k})
        analysis_file = analysis_candidates{k};
        break;
    end
end

if isempty(analysis_file)
    error('analysis_temp_results.mat introuvable.');
end

SA = load(analysis_file,'Positions','Adjacency','time_values');

if ~isfield(SA,'Positions') || ~isfield(SA,'Adjacency')
    error('analysis_temp_results.mat doit contenir Positions et Adjacency.');
end

if isfield(SA,'time_values')
    time_analysis = double(SA.time_values(:));
    if numel(time_analysis) ~= Nt || max(abs(time_analysis-time_values)) > 1e-9
        error('Grille temporelle incompatible dans analysis_temp_results.mat.');
    end
end

%% ============================================================
%  5. Taille moyenne locale d'une composante
%% ============================================================

mean_component_size_th = nan(Nt,Nb);

valid_beta = ...
    isfinite(beta0_th) ...
    & beta0_th > 0 ...
    & isfinite(N_sat_th) ...
    & N_sat_th >= 0;

mean_component_size_th(valid_beta) = ...
    N_sat_th(valid_beta)./beta0_th(valid_beta);

%% ============================================================
%  5.a bis Nouvelle variante THEORIQUE avec beta0 PPP
%
% On remplace uniquement l'ancien beta0(t,phi) par beta0_ppp(t,phi).
% N_sat, lambda, v_rel et eta_sweep restent strictement identiques.
%% ============================================================

mean_component_size_ppp = nan(Nt,Nb);

valid_beta_ppp = ...
    isfinite(beta0_ppp) ...
    & beta0_ppp > 0 ...
    & isfinite(N_sat_th) ...
    & N_sat_th >= 0;

mean_component_size_ppp(valid_beta_ppp) = ...
    N_sat_th(valid_beta_ppp)./beta0_ppp(valid_beta_ppp);

%% ============================================================
%  5.b Correction avec beta0 empirique
%
% On conserve N_sat(t,phi), lambda(t,phi), v_rel(t,phi) et
% eta_sweep(t,phi) theoriques, mais on remplace uniquement
%
%   beta0_th(t,phi) -> beta0_emp_true(t,phi).
%
% Ainsi :
%
%   N_C^corr(t,phi) = N_sat^th(t,phi)/beta0_emp(t,phi).
%% ============================================================

mean_component_size_beta0_emp = nan(Nt,Nb);

valid_beta_emp = ...
    isfinite(beta0_emp) ...
    & beta0_emp > 0 ...
    & isfinite(N_sat_th) ...
    & N_sat_th >= 0;

mean_component_size_beta0_emp(valid_beta_emp) = ...
    N_sat_th(valid_beta_emp)./beta0_emp(valid_beta_emp);

%% ============================================================
%  6. Exposant local et probabilite de fusion
%% ============================================================

mu_merge_th = nan(Nt,Nb);

valid_formula = ...
    valid_beta ...
    & isfinite(lambda_th) ...
    & lambda_th >= 0 ...
    & isfinite(v_rel_th) ...
    & v_rel_th >= 0 ...
    & isfinite(eta_sweep_th) ...
    & eta_sweep_th >= 0;

mu_merge_th(valid_formula) = ...
    2*dmax*DeltaT ...
    .* lambda_th(valid_formula) ...
    .* mean_component_size_th(valid_formula) ...
    .* v_rel_th(valid_formula) ...
    .* eta_sweep_th(valid_formula);

mu_merge_th = max(mu_merge_th,0);

p_merge_th = 1-exp(-mu_merge_th);
p_merge_th = min(max(p_merge_th,0),1);

% Disparition d'une barre H0 sous l'hypothese de fusions binaires.
p_disp_fusion_th = 0.5*p_merge_th;

%% ============================================================
%  6.a bis Nouvelle theorie avec beta0 PPP
%% ============================================================

mu_merge_ppp = nan(Nt,Nb);

valid_formula_ppp = ...
    valid_beta_ppp ...
    & isfinite(lambda_th) ...
    & lambda_th >= 0 ...
    & isfinite(v_rel_th) ...
    & v_rel_th >= 0 ...
    & isfinite(eta_sweep_th) ...
    & eta_sweep_th >= 0;

mu_merge_ppp(valid_formula_ppp) = ...
    2*dmax*DeltaT ...
    .* lambda_th(valid_formula_ppp) ...
    .* mean_component_size_ppp(valid_formula_ppp) ...
    .* v_rel_th(valid_formula_ppp) ...
    .* eta_sweep_th(valid_formula_ppp);

mu_merge_ppp = max(mu_merge_ppp,0);

p_merge_ppp = 1-exp(-mu_merge_ppp);
p_merge_ppp = min(max(p_merge_ppp,0),1);

p_disp_fusion_ppp = 0.5*p_merge_ppp;

%% ============================================================
%  6.b Version corrigee avec beta0 empirique
%% ============================================================

mu_merge_beta0_emp = nan(Nt,Nb);

valid_formula_beta0_emp = ...
    valid_beta_emp ...
    & isfinite(lambda_th) ...
    & lambda_th >= 0 ...
    & isfinite(v_rel_th) ...
    & v_rel_th >= 0 ...
    & isfinite(eta_sweep_th) ...
    & eta_sweep_th >= 0;

mu_merge_beta0_emp(valid_formula_beta0_emp) = ...
    2*dmax*DeltaT ...
    .* lambda_th(valid_formula_beta0_emp) ...
    .* mean_component_size_beta0_emp(valid_formula_beta0_emp) ...
    .* v_rel_th(valid_formula_beta0_emp) ...
    .* eta_sweep_th(valid_formula_beta0_emp);

mu_merge_beta0_emp = max(mu_merge_beta0_emp,0);

p_merge_beta0_emp = 1-exp(-mu_merge_beta0_emp);
p_merge_beta0_emp = min(max(p_merge_beta0_emp,0),1);

p_disp_fusion_beta0_emp = 0.5*p_merge_beta0_emp;

%% ============================================================
%  6.c Version beta0 empirique + eta_sweep empirique
%% ============================================================

mu_merge_beta0_eta_emp = nan(Nt,Nb);

valid_formula_beta0_eta_emp = ...
    valid_beta_emp ...
    & isfinite(lambda_th) & lambda_th >= 0 ...
    & isfinite(v_rel_th) & v_rel_th >= 0 ...
    & isfinite(eta_sweep_emp) & eta_sweep_emp >= 0;

mu_merge_beta0_eta_emp(valid_formula_beta0_eta_emp) = ...
    2*dmax*DeltaT ...
    .* lambda_th(valid_formula_beta0_eta_emp) ...
    .* mean_component_size_beta0_emp(valid_formula_beta0_eta_emp) ...
    .* v_rel_th(valid_formula_beta0_eta_emp) ...
    .* eta_sweep_emp(valid_formula_beta0_eta_emp);

mu_merge_beta0_eta_emp = max(mu_merge_beta0_eta_emp,0);
p_merge_beta0_eta_emp = 1-exp(-mu_merge_beta0_eta_emp);
p_merge_beta0_eta_emp = min(max(p_merge_beta0_eta_emp,0),1);
p_disp_fusion_beta0_eta_emp = 0.5*p_merge_beta0_eta_emp;

%% ============================================================
%  6.d Correction par le flux empirique de nouveaux liens
%
% Chaque nouveau lien A_ij(t)=0, A_ij(t+DeltaT)=1 est compte une
% seule fois, avec 1/2 dans la tranche de latitude de chaque extremite.
%
% Le taux incident par composante vaut :
%
%   mu_new,comp^emp(t,b) = 2*N_new^emp(t,b)/beta0^emp(t,b).
%% ============================================================

[n_new_links_emp_t_phi,n_valid_transitions_flux] = ...
    empirical_new_links_t_phi( ...
        SA.Positions,SA.Adjacency,phi_edges,Nt,Nb);

% Flux theorique brut, avant eta_sweep.
n_new_links_th_t_phi = nan(Nt,Nb);

valid_flux_th = ...
    isfinite(N_sat_th) & N_sat_th >= 0 ...
    & isfinite(lambda_th) & lambda_th >= 0 ...
    & isfinite(v_rel_th) & v_rel_th >= 0;

n_new_links_th_t_phi(valid_flux_th) = ...
    0.5 .* N_sat_th(valid_flux_th) ...
    .* 2*dmax*DeltaT ...
    .* lambda_th(valid_flux_th) ...
    .* v_rel_th(valid_flux_th);

% Pas de transition apres le dernier instant.
n_new_links_th_t_phi(end,:) = NaN;

new_link_flux_correction_t_phi = nan(Nt,Nb);
valid_flux_ratio = ...
    isfinite(n_new_links_th_t_phi) ...
    & n_new_links_th_t_phi > 0 ...
    & isfinite(n_new_links_emp_t_phi);

new_link_flux_correction_t_phi(valid_flux_ratio) = ...
    n_new_links_emp_t_phi(valid_flux_ratio) ...
    ./ n_new_links_th_t_phi(valid_flux_ratio);

mu_new_per_component_emp = nan(Nt,Nb);
valid_emp_flux_component = ...
    isfinite(n_new_links_emp_t_phi) ...
    & isfinite(beta0_emp) ...
    & beta0_emp > 0;

mu_new_per_component_emp(valid_emp_flux_component) = ...
    2 .* n_new_links_emp_t_phi(valid_emp_flux_component) ...
    ./ beta0_emp(valid_emp_flux_component);

% Vrai flux + eta theorique
mu_merge_flux_emp_eta_th = nan(Nt,Nb);
valid_flux_eta_th = ...
    valid_emp_flux_component ...
    & isfinite(eta_sweep_th) ...
    & eta_sweep_th >= 0;

mu_merge_flux_emp_eta_th(valid_flux_eta_th) = ...
    mu_new_per_component_emp(valid_flux_eta_th) ...
    .* eta_sweep_th(valid_flux_eta_th);

p_merge_flux_emp_eta_th = ...
    1-exp(-mu_merge_flux_emp_eta_th);
p_merge_flux_emp_eta_th = ...
    min(max(p_merge_flux_emp_eta_th,0),1);

% Vrai flux + eta empirique
mu_merge_flux_eta_emp = nan(Nt,Nb);
valid_flux_eta_emp = ...
    valid_emp_flux_component ...
    & isfinite(eta_sweep_emp) ...
    & eta_sweep_emp >= 0;

mu_merge_flux_eta_emp(valid_flux_eta_emp) = ...
    mu_new_per_component_emp(valid_flux_eta_emp) ...
    .* eta_sweep_emp(valid_flux_eta_emp);

p_merge_flux_eta_emp = ...
    1-exp(-mu_merge_flux_eta_emp);
p_merge_flux_eta_emp = ...
    min(max(p_merge_flux_eta_emp,0),1);

%% ============================================================
%  7. Moyennes globales par composante
%% ============================================================

component_weight = beta0_th;
component_weight(~valid_formula) = 0;

weight_sum = sum(component_weight,2,'omitnan');

p_merge_global_th = nan(Nt,1);
p_disp_fusion_global_th = nan(Nt,1);
mu_merge_global_component_weighted = nan(Nt,1);

valid_time = weight_sum > 0;

p_merge_global_th(valid_time) = ...
    sum(p_merge_th.*component_weight,2,'omitnan') ...
    ./ weight_sum;

p_disp_fusion_global_th(valid_time) = ...
    sum(p_disp_fusion_th.*component_weight,2,'omitnan') ...
    ./ weight_sum;

mu_merge_global_component_weighted(valid_time) = ...
    sum(mu_merge_th.*component_weight,2,'omitnan') ...
    ./ weight_sum;

% Diagnostic : ne pas confondre E[1-exp(-mu)] avec
% 1-exp(-E[mu]).
p_merge_from_mean_mu = ...
    1-exp(-mu_merge_global_component_weighted);

nonlinearity_gap = ...
    p_merge_global_th-p_merge_from_mean_mu;

% ------------------------------------------------------------
% Nouvelle theorie beta0 PPP : ponderation par sa propre masse de
% composantes.
% ------------------------------------------------------------
component_weight_ppp = beta0_ppp;
component_weight_ppp(~valid_formula_ppp) = 0;

weight_sum_ppp = sum(component_weight_ppp,2,'omitnan');

p_merge_global_ppp = nan(Nt,1);
p_disp_fusion_global_ppp = nan(Nt,1);

valid_time_ppp = weight_sum_ppp > 0;

num_merge_ppp = ...
    sum(p_merge_ppp.*component_weight_ppp,2,'omitnan');

p_merge_global_ppp(valid_time_ppp) = ...
    num_merge_ppp(valid_time_ppp) ...
    ./ weight_sum_ppp(valid_time_ppp);

num_disp_ppp = ...
    sum(p_disp_fusion_ppp.*component_weight_ppp,2,'omitnan');

p_disp_fusion_global_ppp(valid_time_ppp) = ...
    num_disp_ppp(valid_time_ppp) ...
    ./ weight_sum_ppp(valid_time_ppp);

% Version corrigee : la ponderation globale utilise elle aussi
% la masse empirique de composantes.
component_weight_beta0_emp = beta0_emp;
component_weight_beta0_emp(~valid_formula_beta0_emp) = 0;

weight_sum_beta0_emp = ...
    sum(component_weight_beta0_emp,2,'omitnan');

p_merge_global_beta0_emp = nan(Nt,1);
p_disp_fusion_global_beta0_emp = nan(Nt,1);

valid_time_beta0_emp = weight_sum_beta0_emp > 0;

p_merge_global_beta0_emp(valid_time_beta0_emp) = ...
    sum(p_merge_beta0_emp.*component_weight_beta0_emp,2,'omitnan') ...
    ./ weight_sum_beta0_emp(valid_time_beta0_emp);

p_disp_fusion_global_beta0_emp(valid_time_beta0_emp) = ...
    sum(p_disp_fusion_beta0_emp.*component_weight_beta0_emp,2,'omitnan') ...
    ./ weight_sum_beta0_emp(valid_time_beta0_emp);

% beta0 + eta empiriques
component_weight_beta0_eta_emp = beta0_emp;
component_weight_beta0_eta_emp(~valid_formula_beta0_eta_emp) = 0;
weight_sum_beta0_eta_emp = sum(component_weight_beta0_eta_emp,2,'omitnan');

p_merge_global_beta0_eta_emp = nan(Nt,1);
valid_time_beta0_eta_emp = weight_sum_beta0_eta_emp > 0;
num_beta0_eta_emp = ...
    sum(p_merge_beta0_eta_emp.*component_weight_beta0_eta_emp,2,'omitnan');

p_merge_global_beta0_eta_emp(valid_time_beta0_eta_emp) = ...
    num_beta0_eta_emp(valid_time_beta0_eta_emp) ...
    ./ weight_sum_beta0_eta_emp(valid_time_beta0_eta_emp);

% beta0 + vrai flux + eta theorique
component_weight_flux_eta_th = beta0_emp;
component_weight_flux_eta_th(~valid_flux_eta_th) = 0;
weight_sum_flux_eta_th = sum(component_weight_flux_eta_th,2,'omitnan');

p_merge_global_flux_emp_eta_th = nan(Nt,1);
valid_time_flux_eta_th = weight_sum_flux_eta_th > 0;
num_flux_eta_th = ...
    sum(p_merge_flux_emp_eta_th.*component_weight_flux_eta_th,2,'omitnan');

p_merge_global_flux_emp_eta_th(valid_time_flux_eta_th) = ...
    num_flux_eta_th(valid_time_flux_eta_th) ...
    ./ weight_sum_flux_eta_th(valid_time_flux_eta_th);

% beta0 + vrai flux + eta empirique
component_weight_flux_eta_emp = beta0_emp;
component_weight_flux_eta_emp(~valid_flux_eta_emp) = 0;
weight_sum_flux_eta_emp = sum(component_weight_flux_eta_emp,2,'omitnan');

p_merge_global_flux_eta_emp = nan(Nt,1);
valid_time_flux_eta_emp = weight_sum_flux_eta_emp > 0;
num_flux_eta_emp = ...
    sum(p_merge_flux_eta_emp.*component_weight_flux_eta_emp,2,'omitnan');

p_merge_global_flux_eta_emp(valid_time_flux_eta_emp) = ...
    num_flux_eta_emp(valid_time_flux_eta_emp) ...
    ./ weight_sum_flux_eta_emp(valid_time_flux_eta_emp);

%% ============================================================
%  8. Figures
%% ============================================================

if make_figures

    figure;
    imagesc(time_values,rad2deg(phi_vals),p_merge_th.');
    axis xy;
    colorbar;
    caxis([0 1]);
    xlabel('Temps (s)');
    ylabel('Latitude \phi (deg)');
    title('p_{merge}^{th}(t,\phi)');

    figure;
    imagesc(time_values,rad2deg(phi_vals),p_disp_fusion_th.');
    axis xy;
    colorbar;
    caxis([0 0.5]);
    xlabel('Temps (s)');
    ylabel('Latitude \phi (deg)');
    title('p_{disp,fusion}^{th}(t,\phi)=p_{merge}^{th}/2');

    n_selected_times = 5;
    selected_indices = unique(round( ...
        linspace(1,Nt,n_selected_times)));

    figure;
    tiledlayout(numel(selected_indices),1, ...
        'TileSpacing','compact','Padding','compact');

    for k = 1:numel(selected_indices)

        it = selected_indices(k);
        nexttile;
        hold on;

        plot(rad2deg(phi_vals),p_merge_th(it,:), ...
            'LineWidth',1.8, ...
            'DisplayName','p_{merge}^{th}, ancien \beta_0');

        plot(rad2deg(phi_vals),p_merge_ppp(it,:), ...
            'LineWidth',2.0, ...
            'DisplayName','p_{merge}^{th}, nouveau \beta_0 PPP');

        plot(rad2deg(phi_vals),p_disp_fusion_th(it,:), ...
            '--','LineWidth',1.5, ...
            'DisplayName','p_{disp,fusion}^{th}');

        grid on;
        ylim([0 1]);
        ylabel('Probabilite');
        title(sprintf('t = %.1f s',time_values(it)));

        if k == 1
            legend('Location','best');
        end

        if k == numel(selected_indices)
            xlabel('Latitude \phi (deg)');
        end

        hold off;
    end

    figure;
    hold on;

    plot(time_values,p_merge_global_th, ...
        'LineWidth',2, ...
        'DisplayName','p_{merge}^{th}, ancien \beta_0');

    plot(time_values,p_merge_global_ppp, ...
        'LineWidth',2.2, ...
        'DisplayName','p_{merge}^{th}, nouveau \beta_0 PPP');

    plot(time_values,p_merge_global_beta0_emp, ...
        '--','LineWidth',2, ...
        'DisplayName','p_{merge}, beta_0 empirique');

    plot(time_values,p_merge_global_beta0_eta_emp, ...
        ':','LineWidth',2, ...
        'DisplayName','p_{merge}, beta_0 + eta empiriques');

    plot(time_values,p_merge_global_flux_emp_eta_th, ...
        '-.','LineWidth',1.8, ...
        'DisplayName','p_{merge}, vrai flux + eta theorique');

    plot(time_values,p_merge_global_flux_eta_emp, ...
        'LineWidth',2.2, ...
        'DisplayName','p_{merge}, vrai flux + eta empirique');

    plot(time_values,p_disp_fusion_global_th, ...
        '--','LineWidth',1.8, ...
        'DisplayName','p_{disp,fusion}^{th} global');

    plot(time_values,p_merge_from_mean_mu, ...
        ':','LineWidth',1.4, ...
        'DisplayName','1-exp(-moyenne composante de \mu)');

    grid on;
    ylim([0 1]);
    xlabel('Temps (s)');
    ylabel('Probabilite');
    title('Probabilites globales ponderees par les composantes');
    legend('Location','best');
    hold off;

    figure;
    imagesc(time_values,rad2deg(phi_vals), ...
        mean_component_size_th.');
    axis xy;
    colorbar;
    xlabel('Temps (s)');
    ylabel('Latitude \phi (deg)');
    title('Taille moyenne locale N_{sat}(t,\phi)/\beta_0(t,\phi)');
end

%% ============================================================
%  9. Affichage console
%% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' P_MERGE(t,phi) - MODELE LOCAL TEMPOREL ET EN LATITUDE\n');
fprintf('============================================================\n');
fprintf('Fichier N                            : %s\n',file_N);
fprintf('Fichier lambda                       : %s\n',file_lambda);
fprintf('Fichier eta_sweep                    : %s\n',file_eta);
fprintf('Fichier beta0 ancien                 : %s\n',file_betti);
fprintf('Fichier beta0 PPP                    : %s\n',file_beta0_new);
fprintf('Fichier v_rel                        : %s\n',file_vrel);
fprintf('dmax                                 : %.10f km\n',dmax);
fprintf('DeltaT                               : %.10f s\n',DeltaT);
fprintf('Nombre d''instants                   : %d\n',Nt);
fprintf('Nombre de tranches                   : %d\n',Nb);
fprintf('------------------------------------------------------------\n');
fprintf('Moyenne p_merge ancien beta0         : %.10f\n', ...
    mean(p_merge_global_th,'omitnan'));
fprintf('Moyenne p_merge nouveau beta0 PPP    : %.10f\n', ...
    mean(p_merge_global_ppp,'omitnan'));
fprintf('Moyenne p_merge avec beta0 empirique : %.10f\n', ...
    mean(p_merge_global_beta0_emp,'omitnan'));
fprintf('Moyenne beta0 + eta empiriques       : %.10f\n', ...
    mean(p_merge_global_beta0_eta_emp,'omitnan'));
fprintf('Moyenne vrai flux + eta theorique    : %.10f\n', ...
    mean(p_merge_global_flux_emp_eta_th,'omitnan'));
fprintf('Moyenne vrai flux + eta empirique    : %.10f\n', ...
    mean(p_merge_global_flux_eta_emp,'omitnan'));

flux_th_global_mean = mean(sum(n_new_links_th_t_phi,2,'omitnan'),'omitnan');
flux_emp_global_mean = mean(sum(n_new_links_emp_t_phi,2,'omitnan'),'omitnan');

fprintf('Nouveaux liens theo / pas            : %.10f\n', ...
    flux_th_global_mean);
fprintf('Nouveaux liens emp / pas             : %.10f\n', ...
    flux_emp_global_mean);
fprintf('Rapport flux theo / emp              : %.10f\n', ...
    flux_th_global_mean/max(flux_emp_global_mean,eps));
fprintf('Moyenne temporelle p_disp fusion     : %.10f\n', ...
    mean(p_disp_fusion_global_th,'omitnan'));
fprintf('Ecart non-linearite moyen            : %.10e\n', ...
    mean(nonlinearity_gap,'omitnan'));
fprintf('============================================================\n');

%% ============================================================
%  10. Sauvegarde
%% ============================================================

output_file = fullfile(search_dir,'pmerge_t_phi_th_results.mat');

if save_results
    save(output_file, ...
        'file_N','file_lambda','file_eta','file_betti','file_vrel', ...
        'dmax','DeltaT', ...
        'time_values','phi_vals','phi_edges','Nt','Nb', ...
        'N_sat_th','lambda_th','eta_sweep_th','eta_sweep_emp', ...
        'beta0_th','beta0_emp','beta0_ppp','v_rel_th', ...
        'file_beta0_new', ...
        'mean_component_size_th','mean_component_size_ppp', ...
        'mean_component_size_beta0_emp', ...
        'mu_merge_th','p_merge_th','p_disp_fusion_th', ...
        'mu_merge_ppp','p_merge_ppp','p_disp_fusion_ppp', ...
        'component_weight_ppp','weight_sum_ppp', ...
        'p_merge_global_ppp','p_disp_fusion_global_ppp', ...
        'component_weight','weight_sum', ...
        'p_merge_global_th','p_disp_fusion_global_th', ...
        'mu_merge_beta0_emp','p_merge_beta0_emp','p_disp_fusion_beta0_emp', ...
        'mu_merge_beta0_eta_emp','p_merge_beta0_eta_emp', ...
        'p_disp_fusion_beta0_eta_emp', ...
        'n_new_links_emp_t_phi','n_new_links_th_t_phi', ...
        'new_link_flux_correction_t_phi','n_valid_transitions_flux', ...
        'mu_new_per_component_emp', ...
        'mu_merge_flux_emp_eta_th','p_merge_flux_emp_eta_th', ...
        'mu_merge_flux_eta_emp','p_merge_flux_eta_emp', ...
        'component_weight_beta0_emp','weight_sum_beta0_emp', ...
        'p_merge_global_beta0_emp','p_disp_fusion_global_beta0_emp', ...
        'component_weight_beta0_eta_emp','weight_sum_beta0_eta_emp', ...
        'p_merge_global_beta0_eta_emp', ...
        'component_weight_flux_eta_th','weight_sum_flux_eta_th', ...
        'p_merge_global_flux_emp_eta_th', ...
        'component_weight_flux_eta_emp','weight_sum_flux_eta_emp', ...
        'p_merge_global_flux_eta_emp', ...
        'mu_merge_global_component_weighted', ...
        'p_merge_from_mean_mu','nonlinearity_gap', ...
        'analysis_file','flux_th_global_mean','flux_emp_global_mean', ...
        '-v7.3');

    fprintf('Resultats sauvegardes dans %s\n',output_file);
end
%% ============================================================
%  Fonctions utilitaires
%% ============================================================

function [n_new_links_t_phi,n_valid_transitions] = ...
    empirical_new_links_t_phi(Positions,Adjacency,phi_edges,Nt,Nb)

    n_new_links_t_phi = nan(Nt,Nb);
    n_new_links_t_phi(1:Nt-1,:) = 0;
    n_valid_transitions = 0;

    Nt_data = min(numel(Positions),numel(Adjacency));
    Nt_loop = min(Nt-1,Nt_data-1);

    for t = 1:Nt_loop
        pos = Positions{t};
        A0 = Adjacency{t};
        A1 = Adjacency{t+1};

        if isempty(pos) || isempty(A0) || isempty(A1)
            n_new_links_t_phi(t,:) = NaN;
            continue;
        end

        pos = double(pos);
        A0 = logical(A0);
        A1 = logical(A1);

        new_edges = triu(A1 & ~A0,1);
        [ii,jj] = find(new_edges);

        n_valid_transitions = n_valid_transitions + 1;

        if isempty(ii)
            continue;
        end

        radius = sqrt(sum(pos.^2,2));
        phi_sat = asin(max(min(pos(:,3)./radius,1),-1));

        bin_i = discretize(phi_sat(ii),phi_edges);
        bin_j = discretize(phi_sat(jj),phi_edges);

        for k = 1:numel(ii)
            bi = bin_i(k);
            bj = bin_j(k);

            if ~isnan(bi)
                n_new_links_t_phi(t,bi) = ...
                    n_new_links_t_phi(t,bi) + 0.5;
            end
            if ~isnan(bj)
                n_new_links_t_phi(t,bj) = ...
                    n_new_links_t_phi(t,bj) + 0.5;
            end
        end
    end
end

function file_path = find_result_file(search_dir,file_name)

    candidates = { ...
        fullfile(search_dir,file_name), ...
        fullfile(search_dir,'..',file_name)};

    file_path = '';

    for k = 1:numel(candidates)
        if isfile(candidates{k})
            file_path = candidates{k};
            break;
        end
    end

    if isempty(file_path)
        error('Fichier introuvable : %s',file_name);
    end
end

function check_fields(S,fields,file_name)

    for k = 1:numel(fields)
        if ~isfield(S,fields{k})
            error('Le fichier %s doit contenir %s.', ...
                file_name,fields{k});
        end
    end
end

function assert_same_grid(time_ref,phi_ref,time_other,phi_other,file_name)

    if numel(time_other) ~= numel(time_ref) ...
            || max(abs(time_other-time_ref)) > 1e-9
        error('Grille temporelle incompatible dans %s.',file_name);
    end

    if numel(phi_other) ~= numel(phi_ref) ...
            || max(abs(phi_other-phi_ref)) > 1e-12
        error('Grille de latitude incompatible dans %s.',file_name);
    end
end

function check_size(X,Nt,Nb,name)

    if ~isequal(size(X),[Nt,Nb])
        error('%s doit etre de taille %d x %d.',name,Nt,Nb);
    end
end

function value = get_first_scalar_field(structures,field_name)

    value = [];

    for k = 1:numel(structures)

        S = structures{k};

        if isfield(S,field_name)

            candidate = double(S.(field_name));

            if isscalar(candidate) ...
                    && isfinite(candidate) ...
                    && candidate > 0

                value = candidate;
                break;
            end
        end
    end

    if isempty(value)
        error('Impossible de recuperer une valeur valide de %s.', ...
            field_name);
    end
end
