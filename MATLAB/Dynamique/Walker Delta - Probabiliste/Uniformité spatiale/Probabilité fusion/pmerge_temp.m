%% pmerge_temp.m
% Comparaison temporelle de p_merge theorique et empirique
% pour le Walker-Delta a uniformite spatiale.
%
% Entrees :
%   - pmerge_th_results.mat
%   - pmerge_emp_results.mat
%
% Sortie :
%   - figure de comparaison avec les valeurs moyennes
%   - pmerge_temp_results.mat

clear; clc; close all;

%% ============================================================
%  LOCALISATION DES FICHIERS
%% ============================================================

script_dir = fileparts(mfilename('fullpath'));
th_file = fullfile(script_dir, 'pmerge_th_results.mat');
emp_file = fullfile(script_dir, 'pmerge_emp_results.mat');

if isempty(th_file)
    error('Fichier pmerge_th_results.mat introuvable.');
end

if isempty(emp_file)
    error('Fichier pmerge_emp_results.mat introuvable.');
end

fprintf('Fichier theorique : %s\n', th_file);
fprintf('Fichier empirique : %s\n', emp_file);

Sth = load(th_file);
Semp = load(emp_file);

%% ============================================================
%  VERIFICATION ET EXTRACTION DES VARIABLES
%% ============================================================

required_th = {'t_transition', 'p_merge_t'};
required_emp = {'time_transition', 'p_merge'};

for k = 1:numel(required_th)
    if ~isfield(Sth, required_th{k})
        error('Variable %s absente du fichier theorique.', ...
            required_th{k});
    end
end

for k = 1:numel(required_emp)
    if ~isfield(Semp, required_emp{k})
        error('Variable %s absente du fichier empirique.', ...
            required_emp{k});
    end
end

t_th = double(Sth.t_transition(:));
p_merge_th = double(Sth.p_merge_t(:));

t_emp = double(Semp.time_transition(:));
p_merge_emp = double(Semp.p_merge(:));

if numel(t_th) ~= numel(p_merge_th)
    error('Tailles incompatibles entre t_transition et p_merge_t.');
end

if numel(t_emp) ~= numel(p_merge_emp)
    error('Tailles incompatibles entre time_transition et p_merge.');
end

%% Courbe empirique lissee
if isfield(Semp, 'p_merge_moving')
    p_merge_emp_smooth = double(Semp.p_merge_moving(:));
else
    moving_window = 15;
    p_merge_emp_smooth = movmean( ...
        p_merge_emp, moving_window, 'Endpoints', 'shrink');
end

%% ============================================================
%  ALIGNEMENT DE LA THEORIE SUR LA GRILLE EMPIRIQUE
%% ============================================================

p_merge_th_aligned = interp1( ...
    t_th, p_merge_th, t_emp, 'linear', NaN);

valid = isfinite(t_emp) ...
    & isfinite(p_merge_emp) ...
    & isfinite(p_merge_emp_smooth) ...
    & isfinite(p_merge_th_aligned);

if ~any(valid)
    error('Aucun point temporel commun valide entre theorie et empirique.');
end

t_common = t_emp(valid);
p_th = p_merge_th_aligned(valid);
p_emp = p_merge_emp(valid);
p_emp_smooth = p_merge_emp_smooth(valid);

%% ============================================================
%  VALEURS MOYENNES
%% ============================================================

p_merge_th_mean = mean(p_th, 'omitnan');

% On privilegie la moyenne globale ponderee du calcul empirique.
if isfield(Semp, 'p_merge_mean')
    p_merge_emp_mean = double(Semp.p_merge_mean);
    emp_mean_source = 'moyenne globale ponderee du fichier empirique';
else
    p_merge_emp_mean = mean(p_emp, 'omitnan');
    emp_mean_source = 'moyenne temporelle simple';
end

p_merge_emp_time_mean = mean(p_emp, 'omitnan');
p_merge_th_time_mean = mean(p_th, 'omitnan');

%% Indicateurs d'ecart
rmse = sqrt(mean((p_emp_smooth - p_th).^2, 'omitnan'));
mae = mean(abs(p_emp_smooth - p_th), 'omitnan');

%% ============================================================
%  GRAPHE PRINCIPAL
%% ============================================================

figure;
hold on;
grid on;
box on;

plot(t_common, p_emp, '-', ...
    'LineWidth', 0.7, ...
    'DisplayName', 'p_{merge}^{emp}(t) brut');

plot(t_common, p_emp_smooth, ...
    'LineWidth', 2.0, ...
    'DisplayName', 'p_{merge}^{emp}(t) lisse');

plot(t_common, p_th, '--', ...
    'LineWidth', 2.2, ...
    'DisplayName', 'p_{merge}^{th}(t)');

yline(p_merge_emp_mean, ':', ...
    'LineWidth', 1.6, ...
    'DisplayName', sprintf( ...
        'Moyenne empirique = %.4f', p_merge_emp_mean));

yline(p_merge_th_mean, '-.', ...
    'LineWidth', 1.6, ...
    'DisplayName', sprintf( ...
        'Moyenne theorique = %.4f', p_merge_th_mean));

xlabel('Temps (s)');
ylabel('Probabilite de fusion');

if isfield(Semp, 'inc_deg')
    inc_deg = double(Semp.inc_deg);
elseif isfield(Sth, 'inc_deg')
    inc_deg = double(Sth.inc_deg);
else
    inc_deg = NaN;
end

if isfinite(inc_deg)
    title(sprintf( ...
        'Walker-Delta spatial : p_{merge}^{th} et p_{merge}^{emp}, i = %.1f deg', ...
        inc_deg));
else
    title('Walker-Delta spatial : comparaison de p_{merge}^{th} et p_{merge}^{emp}');
end

legend('Location', 'best');
hold off;

%% ============================================================
%  AFFICHAGE CONSOLE
%% ============================================================

fprintf('\n=== Comparaison p_merge theorique / empirique ===\n');
fprintf('Nombre de points compares          : %d\n', numel(t_common));
fprintf('Moyenne theorique                  : %.8f\n', ...
    p_merge_th_mean);
fprintf('Moyenne empirique                  : %.8f\n', ...
    p_merge_emp_mean);
fprintf('Source moyenne empirique           : %s\n', ...
    emp_mean_source);
fprintf('Moyenne temporelle empirique       : %.8f\n', ...
    p_merge_emp_time_mean);
fprintf('RMSE theorie / empirique lisse     : %.8f\n', rmse);
fprintf('MAE theorie / empirique lisse      : %.8f\n', mae);

%% ============================================================
%  SAUVEGARDE
%% ============================================================

save('pmerge_temp_results.mat', ...
    't_common', ...
    'p_th', 'p_emp', 'p_emp_smooth', ...
    'p_merge_th_mean', 'p_merge_emp_mean', ...
    'p_merge_th_time_mean', 'p_merge_emp_time_mean', ...
    'rmse', 'mae', 'emp_mean_source', ...
    'th_file', 'emp_file');

fprintf('\nResultats sauvegardes dans pmerge_temp_results.mat\n');

%% ============================================================
%  FONCTION LOCALE
%% ============================================================

function file = first_existing_file(candidates)
    file = '';

    for k = 1:numel(candidates)
        if isfile(candidates{k})
            file = candidates{k};
            return;
        end
    end
end
