%% lambda_t_phi_compare.m
% Comparaison theorie / empirique de la densite satellitaire locale
%
%   lambda(t,phi)   [satellites / km^2]
%
% a partir de plink_t_phi_results.mat produit par plink_t_phi.m.
%
% Theorie continue :
%
%   lambda^th(t,phi)
%      = N f_Phi(t,phi) / (2*pi*R^2*cos(phi)).
%
% Comparaison par tranche b :
%
%   lambda_b^th(t)
%      = N P(phi dans b a l'instant t) / A_b,
%
% avec
%
%   A_b = 2*pi*R^2 [sin(phi_b^+) - sin(phi_b^-)].
%
% Empirique, pour chaque simulation r :
%
%   lambda_b^{emp,(r)}(t)
%      = N_b^{(r)}(t) / A_b.
%
% La courbe empirique affichee est la moyenne sur les simulations,
% accompagnee de l'erreur standard de la moyenne.
%
% Entree :
%   plink_t_phi_results.mat
%
% Sortie :
%   lambda_t_phi_results.mat

clear; clc; close all;

%% ============================================================
%  1. Dossier du script et chargement
%% ============================================================

script_dir = fileparts(mfilename('fullpath'));

if isempty(script_dir)
    script_dir = pwd;
end

input_file = fullfile(script_dir,'plink_t_phi_results.mat');

if ~isfile(input_file)
    error('Fichier introuvable : %s',input_file);
end

S = load(input_file);

required_fields = { ...
    'N','R','time_theory', ...
    'phi_edges_emp','phi_vals_emp','dphi_emp', ...
    'phi_vals_th','dphi_th', ...
    'f_phi_th_fine', ...
    'f_phi_mass_on_emp', ...
    'satellite_count_emp_iterations'};

for k = 1:numel(required_fields)
    if ~isfield(S,required_fields{k})
        error('Le fichier %s doit contenir %s.', ...
            input_file,required_fields{k});
    end
end

N = double(S.N);
R = double(S.R);

time_values = double(S.time_theory(:));
phi_edges_emp = double(S.phi_edges_emp(:).');
phi_vals_emp = double(S.phi_vals_emp(:).');
dphi_emp = double(S.dphi_emp(:).');

phi_vals_th = double(S.phi_vals_th(:).');
dphi_th = double(S.dphi_th(:).');

f_phi_th_fine = double(S.f_phi_th_fine);
f_phi_mass_on_emp = double(S.f_phi_mass_on_emp);

satellite_count_emp_iterations = ...
    double(S.satellite_count_emp_iterations);

[Nt,Nb] = size(f_phi_mass_on_emp);
n_iterations = size(satellite_count_emp_iterations,1);

if numel(time_values) ~= Nt
    error('Dimension incompatible pour time_theory.');
end

if numel(phi_vals_emp) ~= Nb || numel(dphi_emp) ~= Nb
    error('Dimensions incompatibles pour la grille empirique.');
end

if size(satellite_count_emp_iterations,2) ~= Nt || ...
        size(satellite_count_emp_iterations,3) ~= Nb
    error(['satellite_count_emp_iterations doit etre de taille ', ...
           'n_iterations x Nt x Nb.']);
end

%% ============================================================
%  2. Aire exacte des tranches de latitude
%% ============================================================

area_bin = 2*pi*R^2 .* ...
    (sin(phi_edges_emp(2:end)) ...
    - sin(phi_edges_emp(1:end-1)));

if any(area_bin <= 0)
    error('Les aires de tranche doivent etre strictement positives.');
end

%% ============================================================
%  3. Densite theorique par tranche
%
% On utilise la masse theorique exacte deja calculee :
%
%   P_b^th(t) = int_{tranche b} f_Phi(t,phi) dphi.
%% ============================================================

satellite_count_th = N .* f_phi_mass_on_emp;

lambda_bin_th = ...
    satellite_count_th ./ area_bin;

%% ============================================================
%  4. Densite empirique par simulation
%% ============================================================

lambda_bin_emp_iterations = ...
    satellite_count_emp_iterations ...
    ./ reshape(area_bin,1,1,Nb);

lambda_bin_emp = squeeze(mean( ...
    lambda_bin_emp_iterations,1,'omitnan'));

lambda_bin_emp_std = squeeze(std( ...
    lambda_bin_emp_iterations,0,1,'omitnan'));

n_valid_iterations = squeeze(sum( ...
    isfinite(lambda_bin_emp_iterations),1));

lambda_bin_emp_sem = ...
    lambda_bin_emp_std ...
    ./ sqrt(max(n_valid_iterations,1));

%% ============================================================
%  5. Densite theorique continue sur la grille fine
%
% L'element de surface integre sur la longitude est
%
%   dA = 2*pi*R^2*cos(phi) dphi.
%% ============================================================

cos_phi_th = cos(phi_vals_th);

lambda_th_fine = ...
    N .* f_phi_th_fine ...
    ./ (2*pi*R^2 .* cos_phi_th);

lambda_th_fine(~isfinite(lambda_th_fine)) = NaN;

%% ============================================================
%  6. Verifications globales
%% ============================================================

% Le nombre total de satellites doit etre retrouve en sommant
% densite x aire.
N_th_reconstructed = sum( ...
    lambda_bin_th .* area_bin,2,'omitnan');

N_emp_reconstructed_iterations = squeeze(sum( ...
    lambda_bin_emp_iterations ...
    .* reshape(area_bin,1,1,Nb), ...
    3,'omitnan'));

N_emp_reconstructed = mean( ...
    N_emp_reconstructed_iterations,1,'omitnan').';

max_error_N_th = max(abs(N_th_reconstructed-N));
max_error_N_emp = max(abs(N_emp_reconstructed-N));

%% ============================================================
%  7. Diagnostics theorie / empirique
%% ============================================================

valid_compare = ...
    isfinite(lambda_bin_th) ...
    & isfinite(lambda_bin_emp);

difference = nan(Nt,Nb);
difference(valid_compare) = ...
    lambda_bin_emp(valid_compare) ...
    - lambda_bin_th(valid_compare);

rmse_grid = sqrt(mean( ...
    difference(valid_compare).^2));

mae_grid = mean(abs( ...
    difference(valid_compare)));

bias_grid = mean( ...
    difference(valid_compare));

% Erreurs par latitude.
rmse_by_phi = nan(1,Nb);
mae_by_phi = nan(1,Nb);
bias_by_phi = nan(1,Nb);

for b = 1:Nb
    valid_b = valid_compare(:,b);

    if ~any(valid_b)
        continue;
    end

    err_b = difference(valid_b,b);

    rmse_by_phi(b) = sqrt(mean(err_b.^2));
    mae_by_phi(b) = mean(abs(err_b));
    bias_by_phi(b) = mean(err_b);
end

%% ============================================================
%  8. Figures
%% ============================================================

% Carte theorique par tranche.
figure;
imagesc( ...
    time_values, ...
    rad2deg(phi_vals_emp), ...
    lambda_bin_th.');
axis xy;
colorbar;
xlabel('Temps (s)');
ylabel('Latitude \phi (deg)');
title('\lambda^{th}(t,\phi) par tranche');

% Carte empirique moyenne.
figure;
imagesc( ...
    time_values, ...
    rad2deg(phi_vals_emp), ...
    lambda_bin_emp.');
axis xy;
colorbar;
xlabel('Temps (s)');
ylabel('Latitude \phi (deg)');
title(sprintf( ...
    '\\lambda^{emp}(t,\\phi) moyen sur %d simulations', ...
    n_iterations));

% Carte de l'ecart.
figure;
imagesc( ...
    time_values, ...
    rad2deg(phi_vals_emp), ...
    difference.');
axis xy;
colorbar;
xlabel('Temps (s)');
ylabel('Latitude \phi (deg)');
title('\lambda^{emp}(t,\phi)-\lambda^{th}(t,\phi)');

% Coupes a plusieurs instants.
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

    plot(rad2deg(phi_vals_th), ...
        lambda_th_fine(it,:), ...
        'LineWidth',2, ...
        'DisplayName','Theorie fine');

    plot(rad2deg(phi_vals_emp), ...
        lambda_bin_th(it,:), ...
        's', ...
        'LineWidth',1.1, ...
        'DisplayName','Theorie moyenne par tranche');

    errorbar(rad2deg(phi_vals_emp), ...
        lambda_bin_emp(it,:), ...
        lambda_bin_emp_sem(it,:), ...
        'o-', ...
        'LineWidth',1.2, ...
        'DisplayName','Empirique moyen \pm SEM');

    grid on;
    ylabel('\lambda (sat./km^2)');
    title(sprintf('t = %.1f s',time_values(it)));

    if k == 1
        legend('Location','best');
    end

    if k == numel(selected_indices)
        xlabel('Latitude \phi (deg)');
    end

    hold off;
end

% Moyenne temporelle selon la latitude.
figure;
hold on;

plot(rad2deg(phi_vals_emp), ...
    mean(lambda_bin_th,1,'omitnan'), ...
    'LineWidth',2, ...
    'DisplayName','Theorie');

errorbar(rad2deg(phi_vals_emp), ...
    mean(lambda_bin_emp,1,'omitnan'), ...
    std(lambda_bin_emp,0,1,'omitnan') ...
        / sqrt(Nt), ...
    'o-', ...
    'LineWidth',1.2, ...
    'DisplayName','Empirique moyen temporel');

grid on;
xlabel('Latitude \phi (deg)');
ylabel('\lambda moyenne (sat./km^2)');
title('Densite moyenne selon la latitude');
legend('Location','best');
hold off;

% Erreur selon la latitude.
figure;
hold on;

plot(rad2deg(phi_vals_emp),rmse_by_phi, ...
    'LineWidth',2, ...
    'DisplayName','RMSE');

plot(rad2deg(phi_vals_emp),mae_by_phi, ...
    '--','LineWidth',1.8, ...
    'DisplayName','MAE');

plot(rad2deg(phi_vals_emp),bias_by_phi, ...
    ':','LineWidth',1.8, ...
    'DisplayName','Biais emp-th');

grid on;
xlabel('Latitude \phi (deg)');
ylabel('Erreur de densite (sat./km^2)');
title('Erreur de \lambda(t,\phi) selon la latitude');
legend('Location','best');
hold off;

%% ============================================================
%  9. Affichage console
%% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' DENSITE LAMBDA(t,phi) - THEORIE / EMPIRIQUE\n');
fprintf('============================================================\n');
fprintf('Fichier charge                    : %s\n',input_file);
fprintf('N                                 : %d\n',round(N));
fprintf('R                                 : %.8f km\n',R);
fprintf('Nombre de simulations             : %d\n',n_iterations);
fprintf('Nombre d''instants                : %d\n',Nt);
fprintf('Nombre de tranches                : %d\n',Nb);
fprintf('------------------------------------------------------------\n');
fprintf('RMSE                              : %.10e sat./km^2\n', ...
    rmse_grid);
fprintf('MAE                               : %.10e sat./km^2\n', ...
    mae_grid);
fprintf('Biais moyen emp-theorie           : %.10e sat./km^2\n', ...
    bias_grid);
fprintf('------------------------------------------------------------\n');
fprintf('Erreur max reconstruction N th    : %.10e satellites\n', ...
    max_error_N_th);
fprintf('Erreur max reconstruction N emp   : %.10e satellites\n', ...
    max_error_N_emp);
fprintf('============================================================\n');

%% ============================================================
%  10. Sauvegarde
%% ============================================================

output_file = fullfile(script_dir,'lambda_t_phi_results.mat');

save(output_file, ...
    'input_file', ...
    'N','R','Nt','Nb','n_iterations', ...
    'time_values', ...
    'phi_edges_emp','phi_vals_emp','dphi_emp', ...
    'phi_vals_th','dphi_th', ...
    'area_bin', ...
    'f_phi_th_fine','f_phi_mass_on_emp', ...
    'satellite_count_th', ...
    'satellite_count_emp_iterations', ...
    'lambda_th_fine', ...
    'lambda_bin_th', ...
    'lambda_bin_emp_iterations', ...
    'lambda_bin_emp', ...
    'lambda_bin_emp_std', ...
    'lambda_bin_emp_sem', ...
    'n_valid_iterations', ...
    'N_th_reconstructed', ...
    'N_emp_reconstructed_iterations', ...
    'N_emp_reconstructed', ...
    'max_error_N_th','max_error_N_emp', ...
    'difference','valid_compare', ...
    'rmse_grid','mae_grid','bias_grid', ...
    'rmse_by_phi','mae_by_phi','bias_by_phi', ...
    '-v7.3');

fprintf('Resultats sauvegardes dans %s\n',output_file);
