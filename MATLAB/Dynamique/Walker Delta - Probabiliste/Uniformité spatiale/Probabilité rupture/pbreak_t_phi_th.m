%% pbreak_t_phi_script.m
% Calcul de p_break(t,phi) à partir des résultats locaux.
%
% ATTENTION :
% la formule traitée ici est celle de p_break, et non de p_merge.
%
% Formule locale conditionnelle à une composante non isolée :
%
%   p_break^noniso(t,phi)
%     = 1-exp[-mu_break(t,phi)]
%
% avec
%
%   mu_break(t,phi)
%     = E(t,phi)/(beta0(t,phi)-N1(t,phi))
%       * p_break^lien(t,phi)
%       * p_bridge,bord(t,phi)
%
% et
%
%   p_break^lien(t,phi)
%     ~= (2/pi) * v_rel^lien(t,phi) * DeltaT / dmax.
%
% Déconditionnement :
%
%   p_break(t,phi)
%     = (beta0(t,phi)-N1(t,phi))/beta0(t,phi)
%       * p_break^noniso(t,phi).
%
% Ici :
%
%   p_bridge,bord(t,phi) = eta_sweep(t,phi).
%
% Fichiers requis :
%   liens_t_phi_results.mat
%   betti_t_phi_results.mat
%   N1_t_phi_results.mat
%   eta_sweep_t_phi_results.mat
%   vrel_t_phi_results.mat
%
% Sortie :
%   pbreak_t_phi_results.mat

clear; clc; close all;

%% ============================================================
%  1. Paramètres du script
%% ============================================================

% Laisser vide pour utiliser median(diff(time_values)).
DeltaT = [];

make_figures = true;
save_results = true;

% Dossier contenant les fichiers de résultats.
script_dir = fileparts(mfilename('fullpath'));

if isempty(script_dir)
    script_dir = pwd;
end

search_dir = script_dir;

%% ============================================================
%  2. Chargement des fichiers
%% ============================================================

file_edges = fullfile(script_dir, '..', 'Paramètres', 'Nombre liens', 'liens_t_phi_results.mat');
file_betti = fullfile(script_dir, '..', 'Paramètres', 'Betti', 'betti_t_phi_results.mat');
file_N1 = fullfile(script_dir, '..', 'Paramètres', 'Betti', 'N1_t_phi_results.mat');
file_eta = fullfile(script_dir, '..', 'Paramètres', 'eta_sweep_t_phi_results.mat');
file_vrel = fullfile(script_dir, '..', 'Paramètres', 'Vitesse relative', 'vrel_t_phi_results.mat');

SE = load(file_edges);
SB = load(file_betti);
SN1 = load(file_N1);
Seta = load(file_eta);
SV = load(file_vrel);

check_fields(SE,{ ...
    'time_theory','phi_vals','edges_per_bin_th'},file_edges);

check_fields(SB,{ ...
    'time_values','phi_vals','beta0_th','beta0_emp_true'},file_betti);

check_fields(SN1,{ ...
    'time_values','phi_vals_emp','isolated_count_th'},file_N1);

check_fields(Seta,{ ...
    'time_values','phi_vals_emp','eta_sweep_bin_th','eta_sweep_emp'},file_eta);

check_fields(SV,{ ...
    'time_values','phi_vals_emp','v_rel_th_on_emp'},file_vrel);

%% ============================================================
%  3. Grille commune
%% ============================================================

time_values = double(SB.time_values(:));
phi_vals = double(SB.phi_vals(:).');

Nt = numel(time_values);
Nb = numel(phi_vals);

assert_same_grid( ...
    time_values,phi_vals, ...
    double(SE.time_theory(:)),double(SE.phi_vals(:).'), ...
    file_edges);

assert_same_grid( ...
    time_values,phi_vals, ...
    double(SN1.time_values(:)),double(SN1.phi_vals_emp(:).'), ...
    file_N1);

assert_same_grid( ...
    time_values,phi_vals, ...
    double(Seta.time_values(:)),double(Seta.phi_vals_emp(:).'), ...
    file_eta);

assert_same_grid( ...
    time_values,phi_vals, ...
    double(SV.time_values(:)),double(SV.phi_vals_emp(:).'), ...
    file_vrel);

if isempty(DeltaT)
    if Nt < 2
        error('Impossible de déduire DeltaT avec un seul instant.');
    end

    DeltaT = median(diff(time_values));
end

%% ============================================================
%  4. Grandeurs locales théoriques
%% ============================================================

edges_th = double(SE.edges_per_bin_th);
beta0_th = double(SB.beta0_th);
beta0_emp = double(SB.beta0_emp_true);
N1_th = double(SN1.isolated_count_th);
p_bridge_bord_th = double(Seta.eta_sweep_bin_th);

% Dans ce modele spatial, eta_sweep joue le role de p_bridge,bord.
% Sa mesure empirique est donc utilisee comme correction empirique
% de p_bridge,bord.
p_bridge_bord_emp = double(Seta.eta_sweep_emp);

v_rel_link_th = double(SV.v_rel_th_on_emp);

check_size(edges_th,Nt,Nb,'edges_per_bin_th');
check_size(beta0_th,Nt,Nb,'beta0_th');
check_size(beta0_emp,Nt,Nb,'beta0_emp_true');
check_size(N1_th,Nt,Nb,'isolated_count_th');
check_size(p_bridge_bord_th,Nt,Nb,'eta_sweep_bin_th');
check_size(p_bridge_bord_emp,Nt,Nb,'eta_sweep_emp');
check_size(v_rel_link_th,Nt,Nb,'v_rel_th_on_emp');

dmax = get_first_scalar_field( ...
    {SN1,Seta,SV,SB,SE},'dmax');

%% ============================================================
%  4.b Donnees empiriques pour le flux de ruptures
%
% Un lien rompu verifie :
%   A_ij(t)=1 et A_ij(t+DeltaT)=0.
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
%  5. Probabilité de rupture d'un lien
%% ============================================================

p_break_link_th = ...
    (2/pi) ...
    .* v_rel_link_th ...
    .* DeltaT ...
    ./ dmax;

% La formule est une approximation de premier ordre.
% On borne néanmoins la valeur pour éviter une probabilité > 1.
p_break_link_th = min(max(p_break_link_th,0),1);

%% ============================================================
%  6. Passage à l'échelle d'une composante non isolée
%% ============================================================

nonisolated_component_count_th = ...
    beta0_th-N1_th;

valid_nonisolated = ...
    isfinite(nonisolated_component_count_th) ...
    & nonisolated_component_count_th > 0;

mean_edges_per_nonisolated_component_th = nan(Nt,Nb);

mean_edges_per_nonisolated_component_th(valid_nonisolated) = ...
    edges_th(valid_nonisolated) ...
    ./ nonisolated_component_count_th(valid_nonisolated);

mu_break_nonisolated_th = nan(Nt,Nb);

valid_formula = ...
    valid_nonisolated ...
    & isfinite(mean_edges_per_nonisolated_component_th) ...
    & mean_edges_per_nonisolated_component_th >= 0 ...
    & isfinite(p_break_link_th) ...
    & p_break_link_th >= 0 ...
    & isfinite(p_bridge_bord_th) ...
    & p_bridge_bord_th >= 0;

mu_break_nonisolated_th(valid_formula) = ...
    mean_edges_per_nonisolated_component_th(valid_formula) ...
    .* p_break_link_th(valid_formula) ...
    .* p_bridge_bord_th(valid_formula);

mu_break_nonisolated_th = max(mu_break_nonisolated_th,0);

p_break_nonisolated_th = ...
    1-exp(-mu_break_nonisolated_th);

p_break_nonisolated_th = min(max( ...
    p_break_nonisolated_th,0),1);

%% ============================================================
%  7. Déconditionnement par rapport aux composantes isolées
%% ============================================================

fraction_nonisolated_th = zeros(Nt,Nb);

valid_beta0 = ...
    isfinite(beta0_th) ...
    & beta0_th > 0;

fraction_nonisolated_th(valid_beta0) = ...
    nonisolated_component_count_th(valid_beta0) ...
    ./ beta0_th(valid_beta0);

fraction_nonisolated_th = min(max( ...
    fraction_nonisolated_th,0),1);

p_break_th = ...
    fraction_nonisolated_th ...
    .* p_break_nonisolated_th;

% Version linéarisée, correspondant au nombre moyen de nouvelles
% composantes par composante initiale lorsque mu << 1.
p_break_linear_th = ...
    fraction_nonisolated_th ...
    .* mu_break_nonisolated_th;

%% ============================================================
%  7.b Version corrigee avec beta0 empirique
%
% Seul beta0(t,phi) est remplace par beta0_emp_true(t,phi).
% N1(t,phi), E(t,phi), p_break^lien et p_bridge restent theoriques.
%% ============================================================

nonisolated_component_count_beta0_emp = ...
    beta0_emp-N1_th;

valid_nonisolated_beta0_emp = ...
    isfinite(nonisolated_component_count_beta0_emp) ...
    & nonisolated_component_count_beta0_emp > 0;

mean_edges_per_nonisolated_component_beta0_emp = ...
    nan(Nt,Nb);

mean_edges_per_nonisolated_component_beta0_emp( ...
    valid_nonisolated_beta0_emp) = ...
    edges_th(valid_nonisolated_beta0_emp) ...
    ./ nonisolated_component_count_beta0_emp( ...
        valid_nonisolated_beta0_emp);

mu_break_nonisolated_beta0_emp = nan(Nt,Nb);

valid_formula_beta0_emp = ...
    valid_nonisolated_beta0_emp ...
    & isfinite(mean_edges_per_nonisolated_component_beta0_emp) ...
    & mean_edges_per_nonisolated_component_beta0_emp >= 0 ...
    & isfinite(p_break_link_th) ...
    & p_break_link_th >= 0 ...
    & isfinite(p_bridge_bord_th) ...
    & p_bridge_bord_th >= 0;

mu_break_nonisolated_beta0_emp(valid_formula_beta0_emp) = ...
    mean_edges_per_nonisolated_component_beta0_emp( ...
        valid_formula_beta0_emp) ...
    .* p_break_link_th(valid_formula_beta0_emp) ...
    .* p_bridge_bord_th(valid_formula_beta0_emp);

mu_break_nonisolated_beta0_emp = ...
    max(mu_break_nonisolated_beta0_emp,0);

p_break_nonisolated_beta0_emp = ...
    1-exp(-mu_break_nonisolated_beta0_emp);

p_break_nonisolated_beta0_emp = ...
    min(max(p_break_nonisolated_beta0_emp,0),1);

fraction_nonisolated_beta0_emp = zeros(Nt,Nb);

valid_beta0_emp = ...
    isfinite(beta0_emp) ...
    & beta0_emp > 0;

fraction_nonisolated_beta0_emp(valid_beta0_emp) = ...
    nonisolated_component_count_beta0_emp(valid_beta0_emp) ...
    ./ beta0_emp(valid_beta0_emp);

fraction_nonisolated_beta0_emp = ...
    min(max(fraction_nonisolated_beta0_emp,0),1);

p_break_beta0_emp = ...
    fraction_nonisolated_beta0_emp ...
    .* p_break_nonisolated_beta0_emp;

p_break_linear_beta0_emp = ...
    fraction_nonisolated_beta0_emp ...
    .* mu_break_nonisolated_beta0_emp;

%% ============================================================
%  7.c Version beta0 empirique + p_bridge,bord empirique
%
% On conserve E(t,phi), N1(t,phi) et p_break^lien(t,phi)
% theoriques, mais on remplace :
%
%   p_bridge,bord^th(t,phi) -> p_bridge,bord^emp(t,phi).
%% ============================================================

mu_break_nonisolated_beta0_pbridge_emp = nan(Nt,Nb);

valid_formula_beta0_pbridge_emp = ...
    valid_nonisolated_beta0_emp ...
    & isfinite(mean_edges_per_nonisolated_component_beta0_emp) ...
    & mean_edges_per_nonisolated_component_beta0_emp >= 0 ...
    & isfinite(p_break_link_th) ...
    & p_break_link_th >= 0 ...
    & isfinite(p_bridge_bord_emp) ...
    & p_bridge_bord_emp >= 0;

mu_break_nonisolated_beta0_pbridge_emp( ...
    valid_formula_beta0_pbridge_emp) = ...
    mean_edges_per_nonisolated_component_beta0_emp( ...
        valid_formula_beta0_pbridge_emp) ...
    .* p_break_link_th(valid_formula_beta0_pbridge_emp) ...
    .* p_bridge_bord_emp(valid_formula_beta0_pbridge_emp);

mu_break_nonisolated_beta0_pbridge_emp = ...
    max(mu_break_nonisolated_beta0_pbridge_emp,0);

p_break_nonisolated_beta0_pbridge_emp = ...
    1-exp(-mu_break_nonisolated_beta0_pbridge_emp);

p_break_nonisolated_beta0_pbridge_emp = ...
    min(max(p_break_nonisolated_beta0_pbridge_emp,0),1);

p_break_beta0_pbridge_emp = ...
    fraction_nonisolated_beta0_emp ...
    .* p_break_nonisolated_beta0_pbridge_emp;

%% ============================================================
%  7.d Correction par le flux empirique de ruptures
%
% Le flux theorique brut de ruptures dans une tranche est :
%
%   N_break^th(t,b)
%      = E^th(t,b) * p_break^lien(t,b).
%
% Le flux empirique est mesure directement dans les graphes par :
%
%   A_ij(t)=1, A_ij(t+DeltaT)=0.
%
% Chaque lien rompu est compte une seule fois. Il n'y a donc PAS
% de facteur 2 dans le passage au taux par composante.
%% ============================================================

[n_broken_links_emp_t_phi,n_valid_transitions_flux] = ...
    empirical_broken_links_t_phi( ...
        SA.Positions,SA.Adjacency,phi_vals,Nt,Nb);

n_broken_links_th_t_phi = ...
    edges_th .* p_break_link_th;

% Le dernier instant ne possede pas de transition sortante.
n_broken_links_th_t_phi(end,:) = NaN;

break_flux_correction_t_phi = nan(Nt,Nb);

valid_break_flux_ratio = ...
    isfinite(n_broken_links_th_t_phi) ...
    & n_broken_links_th_t_phi > 0 ...
    & isfinite(n_broken_links_emp_t_phi);

break_flux_correction_t_phi(valid_break_flux_ratio) = ...
    n_broken_links_emp_t_phi(valid_break_flux_ratio) ...
    ./ n_broken_links_th_t_phi(valid_break_flux_ratio);

% Nombre moyen empirique de ruptures de liens par composante non isolee.
mean_broken_links_per_nonisolated_component_emp = nan(Nt,Nb);

valid_emp_break_component = ...
    isfinite(n_broken_links_emp_t_phi) ...
    & isfinite(nonisolated_component_count_beta0_emp) ...
    & nonisolated_component_count_beta0_emp > 0;

mean_broken_links_per_nonisolated_component_emp( ...
    valid_emp_break_component) = ...
    n_broken_links_emp_t_phi(valid_emp_break_component) ...
    ./ nonisolated_component_count_beta0_emp(valid_emp_break_component);

% ------------------------------------------------------------
% 7.d.1 beta0 empirique + vrai flux + p_bridge theorique
% ------------------------------------------------------------

mu_break_nonisolated_flux_emp_pbridge_th = nan(Nt,Nb);

valid_flux_pbridge_th = ...
    valid_emp_break_component ...
    & isfinite(p_bridge_bord_th) ...
    & p_bridge_bord_th >= 0;

mu_break_nonisolated_flux_emp_pbridge_th(valid_flux_pbridge_th) = ...
    mean_broken_links_per_nonisolated_component_emp(valid_flux_pbridge_th) ...
    .* p_bridge_bord_th(valid_flux_pbridge_th);

p_break_nonisolated_flux_emp_pbridge_th = ...
    1-exp(-mu_break_nonisolated_flux_emp_pbridge_th);

p_break_nonisolated_flux_emp_pbridge_th = ...
    min(max(p_break_nonisolated_flux_emp_pbridge_th,0),1);

p_break_flux_emp_pbridge_th = ...
    fraction_nonisolated_beta0_emp ...
    .* p_break_nonisolated_flux_emp_pbridge_th;

% ------------------------------------------------------------
% 7.d.2 beta0 empirique + vrai flux + p_bridge empirique
% ------------------------------------------------------------

mu_break_nonisolated_flux_pbridge_emp = nan(Nt,Nb);

valid_flux_pbridge_emp = ...
    valid_emp_break_component ...
    & isfinite(p_bridge_bord_emp) ...
    & p_bridge_bord_emp >= 0;

mu_break_nonisolated_flux_pbridge_emp(valid_flux_pbridge_emp) = ...
    mean_broken_links_per_nonisolated_component_emp(valid_flux_pbridge_emp) ...
    .* p_bridge_bord_emp(valid_flux_pbridge_emp);

p_break_nonisolated_flux_pbridge_emp = ...
    1-exp(-mu_break_nonisolated_flux_pbridge_emp);

p_break_nonisolated_flux_pbridge_emp = ...
    min(max(p_break_nonisolated_flux_pbridge_emp,0),1);

p_break_flux_pbridge_emp = ...
    fraction_nonisolated_beta0_emp ...
    .* p_break_nonisolated_flux_pbridge_emp;

%% ============================================================
%  8. Moyennes globales pondérées par les composantes
%% ============================================================

component_weight = beta0_th;
component_weight(~valid_beta0) = 0;

weight_sum = sum(component_weight,2,'omitnan');

p_break_global_th = nan(Nt,1);
p_break_nonisolated_global_th = nan(Nt,1);
p_break_linear_global_th = nan(Nt,1);

valid_time = weight_sum > 0;

p_break_global_th(valid_time) = ...
    sum(p_break_th.*component_weight,2,'omitnan') ...
    ./ weight_sum(valid_time);

% Pour la probabilité conditionnelle, on pondère seulement par les
% composantes non isolées.
nonisolated_weight = nonisolated_component_count_th;
nonisolated_weight(~valid_nonisolated) = 0;

nonisolated_weight_sum = sum( ...
    nonisolated_weight,2,'omitnan');

valid_nonisolated_time = nonisolated_weight_sum > 0;

p_break_nonisolated_global_th(valid_nonisolated_time) = ...
    sum( ...
        p_break_nonisolated_th.*nonisolated_weight, ...
        2,'omitnan') ...
    ./ nonisolated_weight_sum(valid_nonisolated_time);

p_break_linear_global_th(valid_time) = ...
    sum(p_break_linear_th.*component_weight,2,'omitnan') ...
    ./ weight_sum(valid_time);

% Integration globale de la version beta0 empirique.
component_weight_beta0_emp = beta0_emp;
component_weight_beta0_emp(~valid_beta0_emp) = 0;

weight_sum_beta0_emp = ...
    sum(component_weight_beta0_emp,2,'omitnan');

p_break_global_beta0_emp = nan(Nt,1);
p_break_nonisolated_global_beta0_emp = nan(Nt,1);
p_break_linear_global_beta0_emp = nan(Nt,1);

valid_time_beta0_emp = weight_sum_beta0_emp > 0;

p_break_global_beta0_emp(valid_time_beta0_emp) = ...
    sum(p_break_beta0_emp.*component_weight_beta0_emp,2,'omitnan') ...
    ./ weight_sum_beta0_emp(valid_time_beta0_emp);

nonisolated_weight_beta0_emp = ...
    nonisolated_component_count_beta0_emp;
nonisolated_weight_beta0_emp(~valid_nonisolated_beta0_emp) = 0;

nonisolated_weight_sum_beta0_emp = ...
    sum(nonisolated_weight_beta0_emp,2,'omitnan');

valid_nonisolated_time_beta0_emp = ...
    nonisolated_weight_sum_beta0_emp > 0;

p_break_nonisolated_global_beta0_emp( ...
    valid_nonisolated_time_beta0_emp) = ...
    sum(p_break_nonisolated_beta0_emp ...
        .* nonisolated_weight_beta0_emp,2,'omitnan') ...
    ./ nonisolated_weight_sum_beta0_emp( ...
        valid_nonisolated_time_beta0_emp);

p_break_linear_global_beta0_emp(valid_time_beta0_emp) = ...
    sum(p_break_linear_beta0_emp ...
        .* component_weight_beta0_emp,2,'omitnan') ...
    ./ weight_sum_beta0_emp(valid_time_beta0_emp);

% beta0 + p_bridge,bord empirique
component_weight_beta0_pbridge_emp = beta0_emp;
component_weight_beta0_pbridge_emp(~valid_formula_beta0_pbridge_emp) = 0;

weight_sum_beta0_pbridge_emp = ...
    sum(component_weight_beta0_pbridge_emp,2,'omitnan');

p_break_global_beta0_pbridge_emp = nan(Nt,1);
valid_time_beta0_pbridge_emp = weight_sum_beta0_pbridge_emp > 0;

num_beta0_pbridge_emp = ...
    sum(p_break_beta0_pbridge_emp ...
        .* component_weight_beta0_pbridge_emp,2,'omitnan');

p_break_global_beta0_pbridge_emp(valid_time_beta0_pbridge_emp) = ...
    num_beta0_pbridge_emp(valid_time_beta0_pbridge_emp) ...
    ./ weight_sum_beta0_pbridge_emp(valid_time_beta0_pbridge_emp);

% beta0 + vrai flux + p_bridge theorique
component_weight_flux_pbridge_th = beta0_emp;
component_weight_flux_pbridge_th(~valid_flux_pbridge_th) = 0;

weight_sum_flux_pbridge_th = ...
    sum(component_weight_flux_pbridge_th,2,'omitnan');

p_break_global_flux_emp_pbridge_th = nan(Nt,1);
valid_time_flux_pbridge_th = weight_sum_flux_pbridge_th > 0;

num_flux_pbridge_th = ...
    sum(p_break_flux_emp_pbridge_th ...
        .* component_weight_flux_pbridge_th,2,'omitnan');

p_break_global_flux_emp_pbridge_th(valid_time_flux_pbridge_th) = ...
    num_flux_pbridge_th(valid_time_flux_pbridge_th) ...
    ./ weight_sum_flux_pbridge_th(valid_time_flux_pbridge_th);

% beta0 + vrai flux + p_bridge empirique
component_weight_flux_pbridge_emp = beta0_emp;
component_weight_flux_pbridge_emp(~valid_flux_pbridge_emp) = 0;

weight_sum_flux_pbridge_emp = ...
    sum(component_weight_flux_pbridge_emp,2,'omitnan');

p_break_global_flux_pbridge_emp = nan(Nt,1);
valid_time_flux_pbridge_emp = weight_sum_flux_pbridge_emp > 0;

num_flux_pbridge_emp = ...
    sum(p_break_flux_pbridge_emp ...
        .* component_weight_flux_pbridge_emp,2,'omitnan');

p_break_global_flux_pbridge_emp(valid_time_flux_pbridge_emp) = ...
    num_flux_pbridge_emp(valid_time_flux_pbridge_emp) ...
    ./ weight_sum_flux_pbridge_emp(valid_time_flux_pbridge_emp);

%% ============================================================
%  9. Figures
%% ============================================================

if make_figures

    figure;
    imagesc(time_values,rad2deg(phi_vals),p_break_th.');
    axis xy;
    colorbar;
    caxis([0 1]);
    xlabel('Temps (s)');
    ylabel('Latitude \phi (deg)');
    title('p_{break}^{th}(t,\phi) déconditionnée');

    figure;
    imagesc( ...
        time_values, ...
        rad2deg(phi_vals), ...
        p_break_nonisolated_th.');
    axis xy;
    colorbar;
    caxis([0 1]);
    xlabel('Temps (s)');
    ylabel('Latitude \phi (deg)');
    title('p_{break}^{th}(t,\phi \mid non isolée)');

    n_selected_times = 5;
    selected_indices = unique(round( ...
        linspace(1,Nt,n_selected_times)));

    figure;
    tiledlayout(numel(selected_indices),1, ...
        'TileSpacing','compact', ...
        'Padding','compact');

    for k = 1:numel(selected_indices)

        it = selected_indices(k);
        nexttile;
        hold on;

        plot(rad2deg(phi_vals), ...
            p_break_nonisolated_th(it,:), ...
            '--','LineWidth',1.5, ...
            'DisplayName','Conditionnelle non isolée');

        plot(rad2deg(phi_vals), ...
            p_break_th(it,:), ...
            'LineWidth',1.8, ...
            'DisplayName','Déconditionnée');

        plot(rad2deg(phi_vals), ...
            p_break_linear_th(it,:), ...
            ':','LineWidth',1.4, ...
            'DisplayName','Version linéaire');

        grid on;
        ylabel('Probabilité');
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

    plot(time_values,p_break_nonisolated_global_th, ...
        '--','LineWidth',1.7, ...
        'DisplayName','Conditionnelle non isolée');

    plot(time_values,p_break_global_th, ...
        'LineWidth',2, ...
        'DisplayName','Déconditionnée');

    plot(time_values,p_break_global_beta0_emp, ...
        '--','LineWidth',2, ...
        'DisplayName','Déconditionnée, beta_0 empirique');

    plot(time_values,p_break_global_beta0_pbridge_emp, ...
        ':','LineWidth',2, ...
        'DisplayName','beta_0 + p_{bridge,bord} empiriques');

    plot(time_values,p_break_global_flux_emp_pbridge_th, ...
        '-.','LineWidth',1.8, ...
        'DisplayName','vrai flux + p_{bridge,bord}^{th}');

    plot(time_values,p_break_global_flux_pbridge_emp, ...
        'LineWidth',2.2, ...
        'DisplayName','vrai flux + p_{bridge,bord}^{emp}');

    plot(time_values,p_break_linear_global_th, ...
        ':','LineWidth',1.5, ...
        'DisplayName','Version linéaire');

    grid on;
    xlabel('Temps (s)');
    ylabel('Probabilité');
    title('Probabilité globale de rupture');
    legend('Location','best');
    hold off;
end

%% ============================================================
%  10. Affichage console
%% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' P_BREAK(t,phi) - MODELE LOCAL TEMPOREL ET EN LATITUDE\n');
fprintf('============================================================\n');
fprintf('Fichier liens                     : %s\n',file_edges);
fprintf('Fichier beta0                     : %s\n',file_betti);
fprintf('Fichier N1                        : %s\n',file_N1);
fprintf('Fichier eta / p_bridge_bord       : %s\n',file_eta);
fprintf('Fichier vitesse relative          : %s\n',file_vrel);
fprintf('dmax                              : %.10f km\n',dmax);
fprintf('DeltaT                            : %.10f s\n',DeltaT);
fprintf('Nombre d''instants                : %d\n',Nt);
fprintf('Nombre de tranches                : %d\n',Nb);
fprintf('------------------------------------------------------------\n');
fprintf('Moyenne p_break conditionnelle    : %.10f\n', ...
    mean(p_break_nonisolated_global_th,'omitnan'));
fprintf('Moyenne p_break déconditionnée    : %.10f\n', ...
    mean(p_break_global_th,'omitnan'));
fprintf('Moyenne p_break beta0 empirique   : %.10f\n', ...
    mean(p_break_global_beta0_emp,'omitnan'));
fprintf('Moyenne beta0 + pbridge empirique : %.10f\n', ...
    mean(p_break_global_beta0_pbridge_emp,'omitnan'));
fprintf('Moyenne vrai flux + pbridge th    : %.10f\n', ...
    mean(p_break_global_flux_emp_pbridge_th,'omitnan'));
fprintf('Moyenne vrai flux + pbridge emp   : %.10f\n', ...
    mean(p_break_global_flux_pbridge_emp,'omitnan'));

flux_break_th_global_mean = ...
    mean(sum(n_broken_links_th_t_phi,2,'omitnan'),'omitnan');
flux_break_emp_global_mean = ...
    mean(sum(n_broken_links_emp_t_phi,2,'omitnan'),'omitnan');

fprintf('Liens rompus theo / pas           : %.10f\n', ...
    flux_break_th_global_mean);
fprintf('Liens rompus emp / pas            : %.10f\n', ...
    flux_break_emp_global_mean);
fprintf('Rapport flux theo / emp           : %.10f\n', ...
    flux_break_th_global_mean/max(flux_break_emp_global_mean,eps));

fprintf('Moyenne p_break linéarisée        : %.10f\n', ...
    mean(p_break_linear_global_th,'omitnan'));
fprintf('============================================================\n');

%% ============================================================
%  11. Sauvegarde
%% ============================================================

output_file = fullfile(search_dir,'pbreak_t_phi_th_results.mat');

if save_results
    save(output_file, ...
        'file_edges','file_betti','file_N1','file_eta','file_vrel', ...
        'dmax','DeltaT', ...
        'time_values','phi_vals','Nt','Nb', ...
        'edges_th','beta0_th','beta0_emp','N1_th', ...
        'p_bridge_bord_th','p_bridge_bord_emp','v_rel_link_th', ...
        'p_break_link_th', ...
        'nonisolated_component_count_th', ...
        'mean_edges_per_nonisolated_component_th', ...
        'mu_break_nonisolated_th', ...
        'p_break_nonisolated_th', ...
        'fraction_nonisolated_th', ...
        'p_break_th','p_break_linear_th', ...
        'component_weight','weight_sum', ...
        'nonisolated_weight','nonisolated_weight_sum', ...
        'p_break_global_th', ...
        'p_break_nonisolated_global_th', ...
        'p_break_linear_global_th', ...
        'nonisolated_component_count_beta0_emp', ...
        'mean_edges_per_nonisolated_component_beta0_emp', ...
        'mu_break_nonisolated_beta0_emp', ...
        'p_break_nonisolated_beta0_emp', ...
        'fraction_nonisolated_beta0_emp', ...
        'p_break_beta0_emp','p_break_linear_beta0_emp', ...
        'component_weight_beta0_emp','weight_sum_beta0_emp', ...
        'nonisolated_weight_beta0_emp','nonisolated_weight_sum_beta0_emp', ...
        'p_break_global_beta0_emp', ...
        'p_break_nonisolated_global_beta0_emp', ...
        'p_break_linear_global_beta0_emp', ...
        'mu_break_nonisolated_beta0_pbridge_emp', ...
        'p_break_nonisolated_beta0_pbridge_emp', ...
        'p_break_beta0_pbridge_emp', ...
        'n_broken_links_emp_t_phi','n_broken_links_th_t_phi', ...
        'break_flux_correction_t_phi','n_valid_transitions_flux', ...
        'mean_broken_links_per_nonisolated_component_emp', ...
        'mu_break_nonisolated_flux_emp_pbridge_th', ...
        'p_break_nonisolated_flux_emp_pbridge_th', ...
        'p_break_flux_emp_pbridge_th', ...
        'mu_break_nonisolated_flux_pbridge_emp', ...
        'p_break_nonisolated_flux_pbridge_emp', ...
        'p_break_flux_pbridge_emp', ...
        'component_weight_beta0_pbridge_emp', ...
        'weight_sum_beta0_pbridge_emp', ...
        'p_break_global_beta0_pbridge_emp', ...
        'component_weight_flux_pbridge_th', ...
        'weight_sum_flux_pbridge_th', ...
        'p_break_global_flux_emp_pbridge_th', ...
        'component_weight_flux_pbridge_emp', ...
        'weight_sum_flux_pbridge_emp', ...
        'p_break_global_flux_pbridge_emp', ...
        'analysis_file', ...
        'flux_break_th_global_mean','flux_break_emp_global_mean', ...
        '-v7.3');

    fprintf('Résultats sauvegardés dans %s\n',output_file);
end
%% ============================================================
%  Fonctions utilitaires
%% ============================================================

function [n_broken_links_t_phi,n_valid_transitions] = ...
    empirical_broken_links_t_phi(Positions,Adjacency,phi_vals,Nt,Nb)

    n_broken_links_t_phi = nan(Nt,Nb);
    n_broken_links_t_phi(1:Nt-1,:) = 0;
    n_valid_transitions = 0;

    % Construction des bords de tranches a partir des centres phi_vals.
    phi_vals = double(phi_vals(:).');
    phi_edges = zeros(1,Nb+1);

    if Nb == 1
        phi_edges = [-pi/2,pi/2];
    else
        phi_edges(2:end-1) = ...
            0.5*(phi_vals(1:end-1)+phi_vals(2:end));

        dphi_left = phi_vals(2)-phi_vals(1);
        dphi_right = phi_vals(end)-phi_vals(end-1);

        phi_edges(1) = phi_vals(1)-0.5*dphi_left;
        phi_edges(end) = phi_vals(end)+0.5*dphi_right;
    end

    Nt_data = min(numel(Positions),numel(Adjacency));
    Nt_loop = min(Nt-1,Nt_data-1);

    for t = 1:Nt_loop

        pos = Positions{t};
        A0 = Adjacency{t};
        A1 = Adjacency{t+1};

        if isempty(pos) || isempty(A0) || isempty(A1)
            n_broken_links_t_phi(t,:) = NaN;
            continue;
        end

        pos = double(pos);
        A0 = logical(A0);
        A1 = logical(A1);

        broken_edges = triu(A0 & ~A1,1);
        [ii,jj] = find(broken_edges);

        n_valid_transitions = n_valid_transitions + 1;

        if isempty(ii)
            continue;
        end

        radius = sqrt(sum(pos.^2,2));
        phi_sat = asin(max(min(pos(:,3)./radius,1),-1));

        bin_i = discretize(phi_sat(ii),phi_edges);
        bin_j = discretize(phi_sat(jj),phi_edges);

        % Chaque lien non oriente compte une seule fois :
        % 1/2 sur chacune de ses extremites.
        for k = 1:numel(ii)

            bi = bin_i(k);
            bj = bin_j(k);

            if ~isnan(bi)
                n_broken_links_t_phi(t,bi) = ...
                    n_broken_links_t_phi(t,bi) + 0.5;
            end

            if ~isnan(bj)
                n_broken_links_t_phi(t,bj) = ...
                    n_broken_links_t_phi(t,bj) + 0.5;
            end
        end
    end
end

function file_path = find_result_file(search_dir,file_name)

    candidates = { ...
        fullfile(search_dir,file_name), ...
        fullfile(search_dir,'..',file_name), ...
        fullfile(search_dir,'Paramètres',file_name), ...
        fullfile(search_dir,'..','Paramètres',file_name), ...
        fullfile(search_dir,'Paramètres','Betti',file_name), ...
        fullfile(search_dir,'..','Paramètres','Betti',file_name), ...
        fullfile(search_dir,'Paramètres','Vitesse relative',file_name), ...
        fullfile(search_dir,'..','Paramètres','Vitesse relative',file_name)};

    file_path = '';

    for k = 1:numel(candidates)
        if isfile(candidates{k})
            file_path = candidates{k};
            return;
        end
    end

    error('Fichier introuvable : %s',file_name);
end

function check_fields(S,fields,file_name)

    for k = 1:numel(fields)
        if ~isfield(S,fields{k})
            error('Le fichier %s doit contenir %s.', ...
                file_name,fields{k});
        end
    end
end

function assert_same_grid( ...
    time_ref,phi_ref,time_other,phi_other,file_name)

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
        error('%s doit être de taille %d x %d.',name,Nt,Nb);
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
                return;
            end
        end
    end

    error('Impossible de récupérer une valeur valide de %s.', ...
        field_name);
end
