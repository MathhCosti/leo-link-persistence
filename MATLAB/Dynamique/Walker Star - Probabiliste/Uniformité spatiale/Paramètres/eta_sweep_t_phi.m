%% eta_sweep_t_phi.m
% Calcul de eta_sweep(t,phi) a partir de lambda_t_phi_results.mat
%
% Modele local :
%
%   eta_sweep(t,phi)
%      = exp[-lambda(t,phi) A_inter(dmax)],
%
% avec
%
%   A_inter(dmax)
%      = (2*pi/3-sqrt(3)/2) dmax^2.
%
% Le script calcule :
%   - eta_sweep theorique sur la grille fine ;
%   - eta_sweep theorique moyenne par tranche ;
%   - eta_sweep empirique pour chaque simulation ;
%   - la moyenne empirique et sa SEM.
%
% IMPORTANT :
% comme l'exponentielle est non lineaire, la valeur empirique principale
% est calculee en moyennant
%
%   exp[-A_inter lambda_emp^(r)(t,phi)]
%
% sur les simulations, et non en appliquant directement l'exponentielle
% a la densite empirique moyenne.
%
% Entree :
%   lambda_t_phi_results.mat
%
% Sortie :
%   eta_sweep_t_phi_results.mat

clear; clc; close all;

%% ============================================================
%  1. Chargement
%% ============================================================

script_dir = fileparts(mfilename('fullpath'));

if isempty(script_dir)
    script_dir = pwd;
end

lambda_file = fullfile(script_dir,'lambda_t_phi_results.mat');

if ~isfile(lambda_file)
    error(['Fichier introuvable : %s\n', ...
           'Execute d''abord lambda_t_phi.m.'],lambda_file);
end

S = load(lambda_file);

required_fields = { ...
    'N','R','Nt','Nb','n_iterations', ...
    'time_values', ...
    'phi_edges_emp','phi_vals_emp','dphi_emp', ...
    'phi_vals_th','dphi_th', ...
    'lambda_th_fine', ...
    'lambda_bin_th', ...
    'lambda_bin_emp_iterations', ...
    'lambda_bin_emp', ...
    'lambda_bin_emp_sem'};

for k = 1:numel(required_fields)
    if ~isfield(S,required_fields{k})
        error('Le fichier %s doit contenir %s.', ...
            lambda_file,required_fields{k});
    end
end

N = double(S.N);
R = double(S.R);

Nt = double(S.Nt);
Nb = double(S.Nb);
n_iterations = double(S.n_iterations);

time_values = double(S.time_values(:));

phi_edges_emp = double(S.phi_edges_emp(:).');
phi_vals_emp = double(S.phi_vals_emp(:).');
dphi_emp = double(S.dphi_emp(:).');

phi_vals_th = double(S.phi_vals_th(:).');
dphi_th = double(S.dphi_th(:).');

lambda_th_fine = double(S.lambda_th_fine);
lambda_bin_th = double(S.lambda_bin_th);

lambda_bin_emp_iterations = ...
    double(S.lambda_bin_emp_iterations);

lambda_bin_emp = double(S.lambda_bin_emp);
lambda_bin_emp_sem = double(S.lambda_bin_emp_sem);

%% ============================================================
%  2. Recuperation de dmax
%% ============================================================

% lambda_t_phi_results.mat ne contient pas necessairement dmax.
% On le cherche d'abord dans le fichier, puis dans le fichier plink
% utilise pour produire lambda(t,phi).

if isfield(S,'dmax')
    dmax = double(S.dmax);

elseif isfield(S,'input_file') && isfile(S.input_file)
    P = load(S.input_file,'dmax');

    if ~isfield(P,'dmax')
        error('Le fichier %s ne contient pas dmax.',S.input_file);
    end

    dmax = double(P.dmax);

else
    plink_file = fullfile(script_dir,'plink_t_phi_results.mat');

    if ~isfile(plink_file)
        error(['Impossible de recuperer dmax. ', ...
               'Le fichier plink_t_phi_results.mat est absent.']);
    end

    P = load(plink_file,'dmax');

    if ~isfield(P,'dmax')
        error('Le fichier %s ne contient pas dmax.',plink_file);
    end

    dmax = double(P.dmax);
end

%% ============================================================
%  3. Aire d'intersection des deux disques de liaison
%% ============================================================

A_intersection = ...
    (2*pi/3-sqrt(3)/2)*dmax^2;

%% ============================================================
%  4. eta_sweep theorique
%% ============================================================

% Valeur ponctuelle sur la grille fine.
eta_sweep_th_fine = exp( ...
    -A_intersection*lambda_th_fine);

% Valeur correspondant a la densite moyenne de chaque tranche.
eta_sweep_bin_th = exp( ...
    -A_intersection*lambda_bin_th);

eta_sweep_th_fine = min(max(eta_sweep_th_fine,0),1);
eta_sweep_bin_th = min(max(eta_sweep_bin_th,0),1);

%% ============================================================
%  5. eta_sweep empirique par simulation
%% ============================================================

eta_sweep_emp_iterations = exp( ...
    -A_intersection*lambda_bin_emp_iterations);

eta_sweep_emp_iterations = min(max( ...
    eta_sweep_emp_iterations,0),1);

eta_sweep_emp = squeeze(mean( ...
    eta_sweep_emp_iterations,1,'omitnan'));

eta_sweep_emp_std = squeeze(std( ...
    eta_sweep_emp_iterations,0,1,'omitnan'));

n_valid_iterations = squeeze(sum( ...
    isfinite(eta_sweep_emp_iterations),1));

eta_sweep_emp_sem = ...
    eta_sweep_emp_std ...
    ./ sqrt(max(n_valid_iterations,1));

% Diagnostic : application de la formule a la densite empirique moyenne.
% Cette quantite n'est pas identique a la moyenne des eta en raison de
% la non-linearite de l'exponentielle.
eta_sweep_from_mean_lambda_emp = exp( ...
    -A_intersection*lambda_bin_emp);

eta_sweep_from_mean_lambda_emp = min(max( ...
    eta_sweep_from_mean_lambda_emp,0),1);

jensen_difference = ...
    eta_sweep_emp-eta_sweep_from_mean_lambda_emp;

%% ============================================================
%  6. Diagnostics theorie / empirique
%% ============================================================

valid_compare = ...
    isfinite(eta_sweep_bin_th) ...
    & isfinite(eta_sweep_emp);

difference = nan(Nt,Nb);

difference(valid_compare) = ...
    eta_sweep_emp(valid_compare) ...
    - eta_sweep_bin_th(valid_compare);

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
    eta_sweep_bin_th.');
axis xy;
colorbar;
caxis([0 1]);
xlabel('Temps (s)');
ylabel('Latitude \phi (deg)');
title('\eta_{sweep}^{th}(t,\phi)');

% Carte empirique moyenne.
figure;
imagesc( ...
    time_values, ...
    rad2deg(phi_vals_emp), ...
    eta_sweep_emp.');
axis xy;
colorbar;
caxis([0 1]);
xlabel('Temps (s)');
ylabel('Latitude \phi (deg)');
title(sprintf( ...
    '\\eta_{sweep}^{emp}(t,\\phi) moyen sur %d simulations', ...
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
title('\eta_{sweep}^{emp}(t,\phi)-\eta_{sweep}^{th}(t,\phi)');

% Coupes temporelles.
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
        eta_sweep_th_fine(it,:), ...
        'LineWidth',2, ...
        'DisplayName','Theorie fine');

    plot(rad2deg(phi_vals_emp), ...
        eta_sweep_bin_th(it,:), ...
        's', ...
        'LineWidth',1.1, ...
        'DisplayName','Theorie par tranche');

    errorbar(rad2deg(phi_vals_emp), ...
        eta_sweep_emp(it,:), ...
        eta_sweep_emp_sem(it,:), ...
        'o-', ...
        'LineWidth',1.2, ...
        'DisplayName','Empirique moyen \pm SEM');

    grid on;
    ylim([0 1]);
    ylabel('\eta_{sweep}');
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
    mean(eta_sweep_bin_th,1,'omitnan'), ...
    'LineWidth',2, ...
    'DisplayName','Theorie');

errorbar(rad2deg(phi_vals_emp), ...
    mean(eta_sweep_emp,1,'omitnan'), ...
    std(eta_sweep_emp,0,1,'omitnan')/sqrt(Nt), ...
    'o-', ...
    'LineWidth',1.2, ...
    'DisplayName','Empirique moyen temporel');

grid on;
ylim([0 1]);
xlabel('Latitude \phi (deg)');
ylabel('\eta_{sweep} moyen');
title('Facteur moyen de redondance spatiale');
legend('Location','best');
hold off;

% Effet de la non-linearite.
figure;
imagesc( ...
    time_values, ...
    rad2deg(phi_vals_emp), ...
    jensen_difference.');
axis xy;
colorbar;
xlabel('Temps (s)');
ylabel('Latitude \phi (deg)');
title(['Moyenne de exp(-A\lambda) ', ...
       '- exp(-A moyenne(\lambda))']);

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
ylabel('Erreur sur \eta_{sweep}');
title('Erreur de \eta_{sweep}(t,\phi) selon la latitude');
legend('Location','best');
hold off;

%% ============================================================
%  8. Affichage console
%% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' ETA_SWEEP(t,phi) DEPUIS LAMBDA(t,phi)\n');
fprintf('============================================================\n');
fprintf('Fichier charge                    : %s\n',lambda_file);
fprintf('N                                 : %d\n',round(N));
fprintf('R                                 : %.8f km\n',R);
fprintf('dmax                              : %.8f km\n',dmax);
fprintf('A_intersection                    : %.10f km^2\n', ...
    A_intersection);
fprintf('Nombre de simulations             : %d\n',n_iterations);
fprintf('Nombre d''instants                : %d\n',Nt);
fprintf('Nombre de tranches                : %d\n',Nb);
fprintf('------------------------------------------------------------\n');
fprintf('RMSE                              : %.10e\n',rmse_grid);
fprintf('MAE                               : %.10e\n',mae_grid);
fprintf('Biais moyen emp-theorie           : %.10e\n',bias_grid);
fprintf('------------------------------------------------------------\n');
fprintf('Moyenne eta theorique             : %.10f\n', ...
    mean(eta_sweep_bin_th,'all','omitnan'));
fprintf('Moyenne eta empirique             : %.10f\n', ...
    mean(eta_sweep_emp,'all','omitnan'));
fprintf('Effet Jensen moyen                : %.10e\n', ...
    mean(jensen_difference,'all','omitnan'));
fprintf('============================================================\n');

%% ============================================================
%  9. Sauvegarde
%% ============================================================

output_file = fullfile( ...
    script_dir, ...
    'eta_sweep_t_phi_results.mat');

save(output_file, ...
    'lambda_file', ...
    'N','R','dmax','A_intersection', ...
    'Nt','Nb','n_iterations', ...
    'time_values', ...
    'phi_edges_emp','phi_vals_emp','dphi_emp', ...
    'phi_vals_th','dphi_th', ...
    'lambda_th_fine','lambda_bin_th', ...
    'lambda_bin_emp_iterations', ...
    'lambda_bin_emp','lambda_bin_emp_sem', ...
    'eta_sweep_th_fine', ...
    'eta_sweep_bin_th', ...
    'eta_sweep_emp_iterations', ...
    'eta_sweep_emp', ...
    'eta_sweep_emp_std', ...
    'eta_sweep_emp_sem', ...
    'eta_sweep_from_mean_lambda_emp', ...
    'jensen_difference', ...
    'difference','valid_compare', ...
    'rmse_grid','mae_grid','bias_grid', ...
    'rmse_by_phi','mae_by_phi','bias_by_phi', ...
    '-v7.3');

fprintf('Resultats sauvegardes dans %s\n',output_file);
