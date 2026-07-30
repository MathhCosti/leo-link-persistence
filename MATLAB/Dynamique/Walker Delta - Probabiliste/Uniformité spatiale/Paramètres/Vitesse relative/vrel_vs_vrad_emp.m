clear; clc; close all;

%% ============================================================
%  SCRIPT DE DIAGNOSTIC : ||v_rel|| vs composante radiale
%
%  Ce script utilise la fonction plot_vrel_vs_vrad.m.
%
%  Il charge les positions temporelles depuis leo_zigzag_analysis_results.mat,
%  reconstruit les vitesses par differences finies, puis trace :
%    - ||v_j - v_i|| moyen
%    - |(v_j - v_i) . r_hat_ij| moyen
%    - [- (v_j - v_i) . r_hat_ij]_+ moyen, utile pour la fusion
%
%  Les paires proches de la frontiere dmax sont les plus pertinentes pour
%  diagnostiquer les creations/ruptures de liens.
%% ============================================================

%% Fichier contenant Positions, time_values et dmax
script_dir = fileparts(mfilename('fullpath'));
results_file = fullfile(script_dir, '..', '..', 'analysis_temp_results.mat');

if ~isfile(results_file)
    error(['Fichier %s introuvable. Lance d''abord analysis_temp.m, ' ...
           'ou place leo_zigzag_analysis_results.mat dans le dossier courant.'], results_file);
end

S = load(results_file);

if ~isfield(S, 'Positions')
    error('Le fichier ne contient pas la variable Positions.');
end

Positions = S.Positions;

if isfield(S, 'time_values')
    time_values = S.time_values(:);
else
    error('Le fichier ne contient pas time_values.');
end

if isfield(S, 'dmax')
    dmax = S.dmax;
else
    error('Le fichier ne contient pas dmax.');
end

fprintf('\n=== Chargement donnees ===\n');
fprintf('Fichier : %s\n', results_file);
fprintf('Nombre de pas temporels : %d\n', numel(time_values));
fprintf('dmax = %.3f km\n', dmax);

%% ============================================================
%  Reconstruction des vitesses par differences finies
%% ============================================================

Nt = numel(Positions);
Velocities = cell(Nt,1);

for k = 1:Nt
    if k == 1
        dt_loc = time_values(k+1) - time_values(k);
        Velocities{k} = (Positions{k+1} - Positions{k}) ./ dt_loc;
    elseif k == Nt
        dt_loc = time_values(k) - time_values(k-1);
        Velocities{k} = (Positions{k} - Positions{k-1}) ./ dt_loc;
    else
        dt_loc = time_values(k+1) - time_values(k-1);
        Velocities{k} = (Positions{k+1} - Positions{k-1}) ./ dt_loc;
    end
end

fprintf('Vitesses reconstruites par differences finies.\n');

%% ============================================================
%  1. Diagnostic sur les liens existants : d_ij <= dmax
%% ============================================================

[t_linked, vrel_linked, vrad_abs_linked, ratio_abs_linked] = vrel_vs_vrad( ...
    Positions, Velocities, time_values, dmax, ...
    'Pairs', 'linked', ...
    'UseAbsRadial', true, ...
    'MakeFigure', true);

%% ============================================================
%  2. Diagnostic pres de la frontiere : dmax +/- 5% dmax
%     C'est le plus pertinent pour les ruptures et creations de liens.
%% ============================================================

boundary_width = 0.05 * dmax;

[t_bound, vrel_bound, vrad_abs_bound, ratio_abs_bound] = vrel_vs_vrad( ...
    Positions, Velocities, time_values, dmax, ...
    'Pairs', 'near_boundary', ...
    'Width', boundary_width, ...
    'UseAbsRadial', true, ...
    'MakeFigure', true);

%% ============================================================
%  3. Diagnostic fusion : seulement la partie qui rapproche
%     [-v_rel . r_hat]_+
%% ============================================================

[t_merge, vrel_merge, vrad_merge, ratio_merge] = vrel_vs_vrad( ...
    Positions, Velocities, time_values, dmax, ...
    'Pairs', 'near_boundary', ...
    'Width', boundary_width, ...
    'UseAbsRadial', false, ...
    'MakeFigure', true);

%% ============================================================
%  Figure de synthese plus lisible
%% ============================================================

figure;
plot(time_values, vrel_bound, 'LineWidth', 1.5); hold on;
plot(time_values, vrad_abs_bound, 'LineWidth', 1.5);
plot(time_values, vrad_merge, 'LineWidth', 1.5);
grid on;
xlabel('Temps (s)');
ylabel('Vitesse moyenne (km/s)');
title('Diagnostic vitesse relative : paires proches de d_{max}');
legend('||v_{rel}|| moyen', ...
       '|v_{rel}\cdot\hat r| moyen', ...
       '[-v_{rel}\cdot\hat r]_+ moyen', ...
       'Location', 'best');

figure;
plot(time_values, ratio_abs_bound, 'LineWidth', 1.5); hold on;
plot(time_values, ratio_merge, 'LineWidth', 1.5);
grid on;
xlabel('Temps (s)');
ylabel('Ratio radial / total');
title('Fraction utile de la vitesse relative pres de d_{max}');
legend('|v_{rel}\cdot\hat r| / ||v_{rel}||', ...
       '[-v_{rel}\cdot\hat r]_+ / ||v_{rel}||', ...
       'Location', 'best');

%% ============================================================
%  Sauvegarde des resultats
%% ============================================================

save('vrel_vs_vrad_emp_results.mat', ...
    'time_values', 'dmax', 'boundary_width', ...
    'vrel_linked', 'vrad_abs_linked', 'ratio_abs_linked', ...
    'vrel_bound', 'vrad_abs_bound', 'ratio_abs_bound', ...
    'vrel_merge', 'vrad_merge', 'ratio_merge');

fprintf('\n=== Resume diagnostic frontiere ===\n');
fprintf('Moyenne ||v_rel|| pres de dmax              : %.6g km/s\n', mean(vrel_bound, 'omitnan'));
fprintf('Moyenne |v_rel . r_hat| pres de dmax        : %.6g km/s\n', mean(vrad_abs_bound, 'omitnan'));
fprintf('Moyenne [-v_rel . r_hat]_+ pres de dmax     : %.6g km/s\n', mean(vrad_merge, 'omitnan'));
fprintf('Ratio moyen radial absolu / total           : %.6g\n', mean(ratio_abs_bound, 'omitnan'));
fprintf('Ratio moyen radial fusion / total           : %.6g\n', mean(ratio_merge, 'omitnan'));
fprintf('\nResultats sauvegardes dans diagnostic_vrel_vrad_results.mat\n');
