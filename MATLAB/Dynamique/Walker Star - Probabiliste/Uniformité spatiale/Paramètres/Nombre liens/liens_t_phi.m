%% liens_t_phi_depuis_plink_results.m
% Calcul et comparaison du nombre local de liens E(t,phi)
% directement a partir de plink_t_phi_results.mat.
%
% Aucun recalcul de p_link(t,phi) n'est effectue.
%
% Theorie :
%   N_b^th(t) = N f_Phi(t,phi_b) Delta phi_b
%
%   E_b^th(t)
%     = 1/2 N_b^th(t) (N-1) p_link^th(t,phi_b)
%
% Empirique, pour chaque simulation r :
%   E_b^{emp,(r)}(t)
%     = 1/2 sum_{i dans b} deg_i^{(r)}(t)
%
% puis :
%   E_b^emp(t) = moyenne_r E_b^{emp,(r)}(t).
%
% Ce choix attribue la moitie de chaque lien a chacune de ses
% extremites. Ainsi,
%
%   sum_b E_b(t) = E_total(t).
%
% Entree :
%   plink_t_phi_results.mat
%
% Sortie :
%   liens_t_phi_results.mat

clear; clc; close all;

%% ============================================================
%  1. Chargement
%% ============================================================

script_dir = fileparts(mfilename('fullpath'));

input_candidates = {
    fullfile(script_dir,'plink_t_phi_results.mat')
    fullfile(script_dir,'..','plink_t_phi_results.mat')
};

input_file = '';

for k = 1:numel(input_candidates)
    if isfile(input_candidates{k})
        input_file = input_candidates{k};
        break;
    end
end

if isempty(input_file)
    error('Fichier plink_t_phi_results.mat introuvable.');
end

S = load(input_file);

required_fields = { ...
    'N', ...
    'time_theory', ...
    'phi_vals', ...
    'dphi', ...
    'f_phi_th', ...
    'p_link_th', ...
    'degree_sum_emp_iterations', ...
    'satellite_count_emp_iterations'};

for k = 1:numel(required_fields)
    if ~isfield(S,required_fields{k})
        error('Le fichier %s doit contenir %s.', ...
            input_file,required_fields{k});
    end
end

N = double(S.N);
time_theory = double(S.time_theory(:));
phi_vals = double(S.phi_vals(:).');
dphi = double(S.dphi(:).');

p_link_th = double(S.p_link_th);
f_phi_th = double(S.f_phi_th);

degree_sum_emp_iterations = ...
    double(S.degree_sum_emp_iterations);

satellite_count_emp_iterations = ...
    double(S.satellite_count_emp_iterations);

[Nt,Nb] = size(p_link_th);

if size(f_phi_th,1) ~= Nt || size(f_phi_th,2) ~= Nb
    error('Dimensions incompatibles entre f_phi_th et p_link_th.');
end

% Les tableaux empiriques doivent etre de taille :
%
%   n_iterations x Nt x Nb.
if ndims(degree_sum_emp_iterations) ~= 3
    error('degree_sum_emp_iterations doit etre un tableau 3D.');
end

if ndims(satellite_count_emp_iterations) ~= 3
    error('satellite_count_emp_iterations doit etre un tableau 3D.');
end

n_iterations = size(degree_sum_emp_iterations,1);

if size(degree_sum_emp_iterations,2) ~= Nt || ...
        size(degree_sum_emp_iterations,3) ~= Nb
    error(['Dimensions incompatibles entre ', ...
           'degree_sum_emp_iterations et p_link_th.']);
end

if ~isequal(size(satellite_count_emp_iterations), ...
        size(degree_sum_emp_iterations))
    error(['satellite_count_emp_iterations et ', ...
           'degree_sum_emp_iterations doivent avoir ', ...
           'les memes dimensions.']);
end

if numel(phi_vals) ~= Nb || numel(dphi) ~= Nb
    error('Dimensions incompatibles pour phi_vals ou dphi.');
end

if numel(time_theory) ~= Nt
    error('Dimension incompatible pour time_theory.');
end

%% ============================================================
%  2. Nombre theorique de satellites par tranche
%% ============================================================

satellite_count_th = ...
    N .* f_phi_th .* dphi;

%% ============================================================
%  3. Nombre local de liens
%
% Chaque lien est partage entre les tranches de ses deux extremites.
%
% Pour chaque realisation r :
%
%   E_b^{emp,(r)}(t)
%     = 1/2 sum_{i dans b} deg_i^{(r)}(t).
%
% La courbe empirique affichee est ensuite la moyenne de cette
% quantite sur les realisations, et non la somme des realisations.
%% ============================================================

edges_per_bin_th = ...
    0.5 .* satellite_count_th ...
    .* (N-1) ...
    .* p_link_th;

% Tableau n_iterations x Nt x Nb.
edges_per_bin_emp_iterations = ...
    0.5 .* degree_sum_emp_iterations;

% Moyenne et dispersion entre les simulations.
edges_per_bin_emp = squeeze(mean( ...
    edges_per_bin_emp_iterations,1,'omitnan'));

edges_per_bin_emp_std = squeeze(std( ...
    edges_per_bin_emp_iterations,0,1,'omitnan'));

n_valid_iterations_per_bin = squeeze(sum( ...
    isfinite(edges_per_bin_emp_iterations),1));

edges_per_bin_emp_sem = ...
    edges_per_bin_emp_std ...
    ./ sqrt(max(n_valid_iterations_per_bin,1));

% Nombre moyen de satellites par tranche, utile comme diagnostic.
satellite_count_emp = squeeze(mean( ...
    satellite_count_emp_iterations,1,'omitnan'));

satellite_count_emp_std = squeeze(std( ...
    satellite_count_emp_iterations,0,1,'omitnan'));

% p_link empirique moyen reconstruit par regroupement des realisations.
degree_sum_emp_total = squeeze(sum( ...
    degree_sum_emp_iterations,1));

satellite_count_emp_total = squeeze(sum( ...
    satellite_count_emp_iterations,1));

p_link_emp = nan(Nt,Nb);

valid_satellite_count = satellite_count_emp_total > 0;

p_link_emp(valid_satellite_count) = ...
    degree_sum_emp_total(valid_satellite_count) ...
    ./ ((N-1)*satellite_count_emp_total(valid_satellite_count));

% Densite de liens par radian de latitude.
edges_density_th = edges_per_bin_th ./ dphi;

edges_density_emp_iterations = ...
    edges_per_bin_emp_iterations ...
    ./ reshape(dphi,1,1,Nb);

edges_density_emp = ...
    edges_per_bin_emp ./ dphi;

edges_density_emp_std = ...
    edges_per_bin_emp_std ./ dphi;

edges_density_emp_sem = ...
    edges_per_bin_emp_sem ./ dphi;

%% ============================================================
%  4. Nombre total de liens
%% ============================================================

edges_total_th = sum(edges_per_bin_th,2,'omitnan');

% Total par realisation, puis moyenne entre simulations.
edges_total_emp_iterations = squeeze(sum( ...
    edges_per_bin_emp_iterations,3,'omitnan'));

edges_total_emp = mean( ...
    edges_total_emp_iterations,1,'omitnan').';

edges_total_emp_std = std( ...
    edges_total_emp_iterations,0,1,'omitnan').';

edges_total_emp_sem = ...
    edges_total_emp_std/sqrt(n_iterations);

% Verification independante depuis p_link global, si disponible.
if isfield(S,'p_link_th_global')
    p_link_th_global = double(S.p_link_th_global(:));
    edges_total_th_from_global = ...
        N*(N-1)/2 .* p_link_th_global;
else
    p_link_th_global = nan(Nt,1);
    edges_total_th_from_global = nan(Nt,1);
end

if isfield(S,'p_link_emp_global_iterations')
    p_link_emp_global_iterations = ...
        double(S.p_link_emp_global_iterations);

    p_link_emp_global = mean( ...
        p_link_emp_global_iterations,1,'omitnan').';

    edges_total_emp_from_global = ...
        N*(N-1)/2 .* p_link_emp_global;
elseif isfield(S,'p_link_emp_global_mean')
    p_link_emp_global = ...
        double(S.p_link_emp_global_mean(:));

    edges_total_emp_from_global = ...
        N*(N-1)/2 .* p_link_emp_global;
else
    p_link_emp_global_iterations = nan(n_iterations,Nt);
    p_link_emp_global = nan(Nt,1);
    edges_total_emp_from_global = nan(Nt,1);
end

%% ============================================================
%  5. Diagnostics
%% ============================================================

valid_compare = ...
    isfinite(edges_per_bin_th) ...
    & isfinite(edges_per_bin_emp);

edges_difference = nan(Nt,Nb);
edges_difference(valid_compare) = ...
    edges_per_bin_emp(valid_compare) ...
    - edges_per_bin_th(valid_compare);

rmse_grid = sqrt(mean( ...
    edges_difference(valid_compare).^2));

mae_grid = mean(abs( ...
    edges_difference(valid_compare)));

bias_grid = mean( ...
    edges_difference(valid_compare));

rmse_by_phi = nan(1,Nb);
mae_by_phi = nan(1,Nb);
bias_by_phi = nan(1,Nb);

for b = 1:Nb
    valid_b = valid_compare(:,b);

    if ~any(valid_b)
        continue;
    end

    err_b = edges_difference(valid_b,b);

    rmse_by_phi(b) = sqrt(mean(err_b.^2));
    mae_by_phi(b) = mean(abs(err_b));
    bias_by_phi(b) = mean(err_b);
end

if all(isfinite(edges_total_th_from_global))
    consistency_th = max(abs( ...
        edges_total_th-edges_total_th_from_global));
else
    consistency_th = NaN;
end

if all(isfinite(edges_total_emp_from_global))
    consistency_emp = max(abs( ...
        edges_total_emp-edges_total_emp_from_global));
else
    consistency_emp = NaN;
end

relative_bias_total = ...
    (mean(edges_total_th,'omitnan') ...
    -mean(edges_total_emp,'omitnan')) ...
    / mean(edges_total_emp,'omitnan');

%% ============================================================
%  6. Figures
%% ============================================================

figure;
imagesc( ...
    time_theory, ...
    rad2deg(phi_vals), ...
    edges_per_bin_th.');
axis xy;
colorbar;
xlabel('Temps (s)');
ylabel('Latitude \phi (deg)');
title('Nombre theorique de liens par tranche E_b^{th}(t)');

figure;
imagesc( ...
    time_theory, ...
    rad2deg(phi_vals), ...
    edges_per_bin_emp.');
axis xy;
colorbar;
xlabel('Temps (s)');
ylabel('Latitude \phi (deg)');
title(sprintf(['Nombre empirique moyen de liens par tranche ', ...
    'sur %d simulations'],n_iterations));

figure;
imagesc( ...
    time_theory, ...
    rad2deg(phi_vals), ...
    edges_difference.');
axis xy;
colorbar;
xlabel('Temps (s)');
ylabel('Latitude \phi (deg)');
title('Ecart E_b^{emp}(t)-E_b^{th}(t)');

% Coupes selon la latitude.
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
        edges_per_bin_th(it,:), ...
        'LineWidth',2, ...
        'DisplayName','Theorie');

    errorbar(rad2deg(phi_vals), ...
        edges_per_bin_emp(it,:), ...
        edges_per_bin_emp_sem(it,:), ...
        'o-', ...
        'LineWidth',1.2, ...
        'DisplayName','Empirique moyen \pm SEM');

    grid on;
    ylabel('Nombre de liens');
    title(sprintf('t = %.1f s',time_theory(it)));

    if k == 1
        legend('Location','best');
    end

    if k == numel(selected_indices)
        xlabel('Latitude \phi (deg)');
    end

    hold off;
end

% Nombre total de liens.
figure;
hold on;

plot(time_theory,edges_total_th, ...
    'LineWidth',2, ...
    'DisplayName','Theorie depuis E_b(t)');

plot(time_theory,edges_total_emp, ...
    'LineWidth',1.5, ...
    'DisplayName','Empirique moyen');

plot(time_theory, ...
    edges_total_emp+edges_total_emp_sem, ...
    ':','LineWidth',1, ...
    'DisplayName','Empirique moyen + SEM');

plot(time_theory, ...
    edges_total_emp-edges_total_emp_sem, ...
    ':','LineWidth',1, ...
    'DisplayName','Empirique moyen - SEM');

if all(isfinite(edges_total_th_from_global))
    plot(time_theory,edges_total_th_from_global, ...
        '--','LineWidth',1.5, ...
        'DisplayName','Theorie depuis p_{link}(t)');
end

grid on;
xlabel('Temps (s)');
ylabel('Nombre total de liens');
title('Nombre total de liens');
legend('Location','best');
hold off;

% Erreurs selon la latitude.
figure;
hold on;

plot(rad2deg(phi_vals),rmse_by_phi, ...
    'LineWidth',2, ...
    'DisplayName','RMSE');

plot(rad2deg(phi_vals),mae_by_phi, ...
    '--','LineWidth',1.8, ...
    'DisplayName','MAE');

plot(rad2deg(phi_vals),bias_by_phi, ...
    ':','LineWidth',1.8, ...
    'DisplayName','Biais emp-th');

grid on;
xlabel('Latitude \phi (deg)');
ylabel('Erreur en nombre de liens');
title('Erreur du nombre local de liens');
legend('Location','best');
hold off;

%% ============================================================
%  7. Affichage console
%% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' NOMBRE DE LIENS E(t,phi) DEPUIS p_link(t,phi)\n');
fprintf('============================================================\n');
fprintf('Fichier charge                     : %s\n',input_file);
fprintf('N                                  : %d\n',round(N));
fprintf('Nombre de simulations empiriques   : %d\n',n_iterations);
fprintf('Nombre d''instants                 : %d\n',Nt);
fprintf('Nombre de tranches                 : %d\n',Nb);
fprintf('------------------------------------------------------------\n');
fprintf('RMSE locale                        : %.10f liens\n', ...
    rmse_grid);
fprintf('MAE locale                         : %.10f liens\n', ...
    mae_grid);
fprintf('Biais local moyen emp-th           : %.10f liens\n', ...
    bias_grid);
fprintf('------------------------------------------------------------\n');
fprintf('Nombre total moyen empirique       : %.10f\n', ...
    mean(edges_total_emp,'omitnan'));
fprintf('Nombre total moyen theorique       : %.10f\n', ...
    mean(edges_total_th,'omitnan'));
fprintf('Biais relatif total theorie/emp    : %.6f %%\n', ...
    100*relative_bias_total);

if isfinite(consistency_th)
    fprintf('Coherence somme locale/theorie glob: %.3e liens\n', ...
        consistency_th);
end

if isfinite(consistency_emp)
    fprintf('Coherence somme locale/emp global  : %.3e liens\n', ...
        consistency_emp);
end

fprintf('============================================================\n');

%% ============================================================
%  8. Sauvegarde
%% ============================================================

output_file = fullfile(script_dir,'liens_t_phi_results.mat');

save(output_file, ...
    'input_file', ...
    'N','Nt','Nb', ...
    'time_theory','phi_vals','dphi', ...
    'p_link_th','p_link_emp', ...
    'f_phi_th', ...
    'n_iterations', ...
    'satellite_count_th', ...
    'satellite_count_emp_iterations', ...
    'satellite_count_emp', ...
    'satellite_count_emp_std', ...
    'satellite_count_emp_total', ...
    'degree_sum_emp_iterations', ...
    'degree_sum_emp_total', ...
    'edges_per_bin_th', ...
    'edges_per_bin_emp_iterations', ...
    'edges_per_bin_emp', ...
    'edges_per_bin_emp_std', ...
    'edges_per_bin_emp_sem', ...
    'n_valid_iterations_per_bin', ...
    'edges_density_th', ...
    'edges_density_emp_iterations', ...
    'edges_density_emp', ...
    'edges_density_emp_std', ...
    'edges_density_emp_sem', ...
    'edges_total_th', ...
    'edges_total_emp_iterations', ...
    'edges_total_emp', ...
    'edges_total_emp_std', ...
    'edges_total_emp_sem', ...
    'p_link_th_global', ...
    'p_link_emp_global_iterations', ...
    'p_link_emp_global', ...
    'edges_total_th_from_global', ...
    'edges_total_emp_from_global', ...
    'edges_difference','valid_compare', ...
    'rmse_grid','mae_grid','bias_grid', ...
    'rmse_by_phi','mae_by_phi','bias_by_phi', ...
    'consistency_th','consistency_emp', ...
    'relative_bias_total');

fprintf('Resultats sauvegardes dans %s\n',output_file);
