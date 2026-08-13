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

SN = load(file_N);
SL = load(file_lambda);
SE = load(file_eta);
SB = load(file_betti);
SV = load(file_vrel);

check_fields(SN,{ ...
    'time_values','phi_vals_emp','phi_edges_emp', ...
    'satellite_count_th'},file_N);

check_fields(SL,{ ...
    'time_values','phi_vals_emp','lambda_bin_th'},file_lambda);

check_fields(SE,{ ...
    'time_values','phi_vals_emp','eta_sweep_bin_th'},file_eta);

check_fields(SB,{ ...
    'time_values','phi_vals','beta0_th'},file_betti);

check_fields(SV,{ ...
    'time_values','phi_vals_emp','v_rel_th_on_emp'},file_vrel);

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
beta0_th = double(SB.beta0_th);
v_rel_th = double(SV.v_rel_th_on_emp);

check_size(N_sat_th,Nt,Nb,'satellite_count_th');
check_size(lambda_th,Nt,Nb,'lambda_bin_th');
check_size(eta_sweep_th,Nt,Nb,'eta_sweep_bin_th');
check_size(beta0_th,Nt,Nb,'beta0_th');
check_size(v_rel_th,Nt,Nb,'v_rel_th_on_emp');

% dmax est prioritairement lu depuis lambda/eta/vrel/betti.
dmax = get_first_scalar_field( ...
    {SL,SE,SV,SB,SN},'dmax');

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
            'DisplayName','p_{merge}^{th}');

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
        'DisplayName','p_{merge}^{th} global');

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
fprintf('Fichier beta0                        : %s\n',file_betti);
fprintf('Fichier v_rel                        : %s\n',file_vrel);
fprintf('dmax                                 : %.10f km\n',dmax);
fprintf('DeltaT                               : %.10f s\n',DeltaT);
fprintf('Nombre d''instants                   : %d\n',Nt);
fprintf('Nombre de tranches                   : %d\n',Nb);
fprintf('------------------------------------------------------------\n');
fprintf('Moyenne temporelle p_merge global    : %.10f\n', ...
    mean(p_merge_global_th,'omitnan'));
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
        'N_sat_th','lambda_th','eta_sweep_th','beta0_th','v_rel_th', ...
        'mean_component_size_th', ...
        'mu_merge_th','p_merge_th','p_disp_fusion_th', ...
        'component_weight','weight_sum', ...
        'p_merge_global_th','p_disp_fusion_global_th', ...
        'mu_merge_global_component_weighted', ...
        'p_merge_from_mean_mu','nonlinearity_gap', ...
        '-v7.3');

    fprintf('Resultats sauvegardes dans %s\n',output_file);
end
%% ============================================================
%  Fonctions utilitaires
%% ============================================================

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
