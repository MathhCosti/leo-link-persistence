%% nombre_satellites_t_phi_compare.m
% Comparaison theorie / empirique du nombre local de satellites
%
%   N_b(t)
%
% dans chaque tranche de latitude b, en fonction du temps.
%
% Theorie continue :
%
%   n_phi^th(t,phi) = N f_Phi(t,phi)
%
% avec n_phi exprime en satellites par radian de latitude.
%
% Pour une tranche b = [phi_b^-,phi_b^+] :
%
%   N_b^th(t)
%     = N int_{phi_b^-}^{phi_b^+} f_Phi(t,phi) dphi.
%
% Empirique, pour chaque simulation r :
%
%   N_b^{emp,(r)}(t)
%     = nombre de satellites presents dans la tranche b.
%
% La comparaison locale est effectuee sur les memes tranches de
% latitude que celles utilisees dans plink_t_phi.m.
%
% Entree :
%   plink_t_phi_results.mat
%
% Sortie :
%   nombre_satellites_t_phi_results.mat

clear; clc; close all;

%% ============================================================
%  1. Chargement
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
    'N', ...
    'time_theory', ...
    'phi_edges_emp', ...
    'phi_vals_emp', ...
    'dphi_emp', ...
    'phi_vals_th', ...
    'dphi_th', ...
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
%  2. Nombre theorique de satellites par tranche
%
% f_phi_mass_on_emp contient deja :
%
%   int_{tranche b} f_Phi(t,phi) dphi.
%% ============================================================

satellite_count_th = ...
    N .* f_phi_mass_on_emp;

%% ============================================================
%  3. Nombre empirique moyen sur les simulations
%% ============================================================

satellite_count_emp = squeeze(mean( ...
    satellite_count_emp_iterations,1,'omitnan'));

satellite_count_emp_std = squeeze(std( ...
    satellite_count_emp_iterations,0,1,'omitnan'));

n_valid_iterations = squeeze(sum( ...
    isfinite(satellite_count_emp_iterations),1));

satellite_count_emp_sem = ...
    satellite_count_emp_std ...
    ./ sqrt(max(n_valid_iterations,1));

%% ============================================================
%  4. Densite en latitude : satellites par radian
%
% Cette quantite est utile pour visualiser la forme continue :
%
%   n_phi^th(t,phi) = N f_Phi(t,phi).
%% ============================================================

satellite_density_phi_th_fine = ...
    N .* f_phi_th_fine;

satellite_density_phi_th_bin = ...
    satellite_count_th ./ dphi_emp;

satellite_density_phi_emp_iterations = ...
    satellite_count_emp_iterations ...
    ./ reshape(dphi_emp,1,1,Nb);

satellite_density_phi_emp = squeeze(mean( ...
    satellite_density_phi_emp_iterations,1,'omitnan'));

satellite_density_phi_emp_std = squeeze(std( ...
    satellite_density_phi_emp_iterations,0,1,'omitnan'));

satellite_density_phi_emp_sem = ...
    satellite_density_phi_emp_std ...
    ./ sqrt(max(n_valid_iterations,1));

%% ============================================================
%  5. Verifications de conservation du nombre total
%% ============================================================

N_total_th = sum( ...
    satellite_count_th,2,'omitnan');

N_total_emp_iterations = squeeze(sum( ...
    satellite_count_emp_iterations,3,'omitnan'));

N_total_emp = mean( ...
    N_total_emp_iterations,1,'omitnan').';

N_total_emp_std = std( ...
    N_total_emp_iterations,0,1,'omitnan').';

N_total_emp_sem = ...
    N_total_emp_std/sqrt(n_iterations);

max_error_total_th = max(abs(N_total_th-N));
max_error_total_emp = max(abs(N_total_emp-N));

%% ============================================================
%  6. Diagnostics theorie / empirique
%% ============================================================

valid_compare = ...
    isfinite(satellite_count_th) ...
    & isfinite(satellite_count_emp);

difference = nan(Nt,Nb);

difference(valid_compare) = ...
    satellite_count_emp(valid_compare) ...
    - satellite_count_th(valid_compare);

rmse_grid = sqrt(mean( ...
    difference(valid_compare).^2));

mae_grid = mean(abs( ...
    difference(valid_compare)));

bias_grid = mean( ...
    difference(valid_compare));

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
%  7. Figures
%% ============================================================

% Carte theorique.
figure;
imagesc( ...
    time_values, ...
    rad2deg(phi_vals_emp), ...
    satellite_count_th.');
axis xy;
colorbar;
xlabel('Temps (s)');
ylabel('Latitude \phi (deg)');
title('Nombre theorique de satellites par tranche N_b^{th}(t)');

% Carte empirique moyenne.
figure;
imagesc( ...
    time_values, ...
    rad2deg(phi_vals_emp), ...
    satellite_count_emp.');
axis xy;
colorbar;
xlabel('Temps (s)');
ylabel('Latitude \phi (deg)');
title(sprintf( ...
    'Nombre empirique moyen de satellites sur %d simulations', ...
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
title('Ecart N_b^{emp}(t)-N_b^{th}(t)');

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

    % Courbe continue exprimee en nombre attendu dans une tranche
    % empirique de largeur dphi locale.
    continuous_count_equivalent = ...
        satellite_density_phi_th_fine(it,:) ...
        .* mean(dphi_emp);

    plot(rad2deg(phi_vals_th), ...
        continuous_count_equivalent, ...
        'LineWidth',2, ...
        'DisplayName','Theorie fine (equivalent tranche)');

    plot(rad2deg(phi_vals_emp), ...
        satellite_count_th(it,:), ...
        's', ...
        'LineWidth',1.1, ...
        'DisplayName','Theorie integree par tranche');

    errorbar(rad2deg(phi_vals_emp), ...
        satellite_count_emp(it,:), ...
        satellite_count_emp_sem(it,:), ...
        'o-', ...
        'LineWidth',1.2, ...
        'DisplayName','Empirique moyen \pm SEM');

    grid on;
    ylabel('Nombre de satellites');
    title(sprintf('t = %.1f s',time_values(it)));

    if k == 1
        legend('Location','best');
    end

    if k == numel(selected_indices)
        xlabel('Latitude \phi (deg)');
    end

    hold off;
end

% Nombre total de satellites.
figure;
hold on;

plot(time_values,N_total_th, ...
    'LineWidth',2, ...
    'DisplayName','Theorie');

plot(time_values,N_total_emp, ...
    'LineWidth',1.5, ...
    'DisplayName','Empirique moyen');

plot(time_values,N_total_emp+N_total_emp_sem, ...
    ':','LineWidth',1, ...
    'DisplayName','Empirique + SEM');

plot(time_values,N_total_emp-N_total_emp_sem, ...
    ':','LineWidth',1, ...
    'DisplayName','Empirique - SEM');

yline(N,'--','N impose', ...
    'HandleVisibility','off');

grid on;
xlabel('Temps (s)');
ylabel('Nombre total de satellites');
title('Verification de la conservation de N');
legend('Location','best');
hold off;

% Moyenne temporelle selon la latitude.
figure;
hold on;

plot(rad2deg(phi_vals_emp), ...
    mean(satellite_count_th,1,'omitnan'), ...
    'LineWidth',2, ...
    'DisplayName','Theorie');

errorbar(rad2deg(phi_vals_emp), ...
    mean(satellite_count_emp,1,'omitnan'), ...
    std(satellite_count_emp,0,1,'omitnan') ...
        / sqrt(Nt), ...
    'o-', ...
    'LineWidth',1.2, ...
    'DisplayName','Empirique moyen temporel');

grid on;
xlabel('Latitude \phi (deg)');
ylabel('Nombre moyen de satellites');
title('Nombre moyen de satellites par tranche');
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
ylabel('Erreur en nombre de satellites');
title('Erreur de N_b(t) selon la latitude');
legend('Location','best');
hold off;

%% ============================================================
%  8. Affichage console
%% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' NOMBRE DE SATELLITES N_b(t,phi) - THEORIE / EMPIRIQUE\n');
fprintf('============================================================\n');
fprintf('Fichier charge                    : %s\n',input_file);
fprintf('N impose                          : %d\n',round(N));
fprintf('Nombre de simulations             : %d\n',n_iterations);
fprintf('Nombre d''instants                : %d\n',Nt);
fprintf('Nombre de tranches                : %d\n',Nb);
fprintf('------------------------------------------------------------\n');
fprintf('RMSE locale                       : %.10e satellites\n', ...
    rmse_grid);
fprintf('MAE locale                        : %.10e satellites\n', ...
    mae_grid);
fprintf('Biais moyen emp-theorie           : %.10e satellites\n', ...
    bias_grid);
fprintf('------------------------------------------------------------\n');
fprintf('Erreur max conservation N th      : %.10e satellites\n', ...
    max_error_total_th);
fprintf('Erreur max conservation N emp     : %.10e satellites\n', ...
    max_error_total_emp);
fprintf('============================================================\n');

%% ============================================================
%  9. Sauvegarde
%% ============================================================

output_file = fullfile( ...
    script_dir, ...
    'N_t_phi_results.mat');

save(output_file, ...
    'input_file', ...
    'N','Nt','Nb','n_iterations', ...
    'time_values', ...
    'phi_edges_emp','phi_vals_emp','dphi_emp', ...
    'phi_vals_th','dphi_th', ...
    'f_phi_th_fine','f_phi_mass_on_emp', ...
    'satellite_count_th', ...
    'satellite_count_emp_iterations', ...
    'satellite_count_emp', ...
    'satellite_count_emp_std', ...
    'satellite_count_emp_sem', ...
    'n_valid_iterations', ...
    'satellite_density_phi_th_fine', ...
    'satellite_density_phi_th_bin', ...
    'satellite_density_phi_emp_iterations', ...
    'satellite_density_phi_emp', ...
    'satellite_density_phi_emp_std', ...
    'satellite_density_phi_emp_sem', ...
    'N_total_th', ...
    'N_total_emp_iterations', ...
    'N_total_emp', ...
    'N_total_emp_std', ...
    'N_total_emp_sem', ...
    'max_error_total_th', ...
    'max_error_total_emp', ...
    'difference','valid_compare', ...
    'rmse_grid','mae_grid','bias_grid', ...
    'rmse_by_phi','mae_by_phi','bias_by_phi', ...
    '-v7.3');

fprintf('Resultats sauvegardes dans %s\n',output_file);
