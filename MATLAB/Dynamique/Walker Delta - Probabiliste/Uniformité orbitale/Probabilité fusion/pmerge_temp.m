clear; clc; close all;

%% ============================================================
%  COMPARAISON TEMPORELLE DE p_merge
%
%  Courbe empirique :
%      chargee depuis pmerge_emp_results.mat
%
%  Moyenne theorique :
%      chargee depuis pmerge_th_results.mat
%
%  Sortie :
%      figure avec p_merge empirique, moyenne glissante,
%      moyenne empirique globale et valeur theorique constante.
%% ============================================================

emp_file = 'pmerge_emp_results.mat';
th_file  = 'pmerge_th_results.mat';

if ~isfile(emp_file)
    error('Fichier empirique introuvable : %s', emp_file);
end

if ~isfile(th_file)
    error('Fichier theorique introuvable : %s', th_file);
end

%% Chargement des resultats empiriques
Semp = load(emp_file, ...
    'time_transition', ...
    'p_merge', ...
    'p_merge_moving', ...
    'p_merge_mean', ...
    'p_merge_time_mean', ...
    'moving_window');

time_transition = Semp.time_transition(:);
p_merge_emp = Semp.p_merge(:);
p_merge_emp_moving = Semp.p_merge_moving(:);
p_merge_emp_mean = Semp.p_merge_mean;
p_merge_emp_time_mean = Semp.p_merge_time_mean;
moving_window = Semp.moving_window;

%% Chargement du resultat theorique
Sth = load(th_file, ...
    'p_merge_mean');

p_merge_th = Sth.p_merge_mean;

%% Verifications
if numel(time_transition) ~= numel(p_merge_emp)
    error(['Les dimensions de time_transition et p_merge ', ...
           'ne correspondent pas.']);
end

if numel(p_merge_emp_moving) ~= numel(p_merge_emp)
    error(['Les dimensions de p_merge_moving et p_merge ', ...
           'ne correspondent pas.']);
end

%% ============================================================
%  Trace principal
%% ============================================================

figure;
hold on;
grid on;
box on;

plot(time_transition, p_merge_emp, ...
    '-', ...
    'LineWidth', 0.8, ...
    'DisplayName', 'p_{merge}^{emp}(t)');

plot(time_transition, p_merge_emp_moving, ...
    'LineWidth', 2.2, ...
    'DisplayName', sprintf( ...
        'Moyenne glissante empirique (%d points)', moving_window));

yline(p_merge_emp_mean, ...
    ':', ...
    'LineWidth', 1.8, ...
    'DisplayName', sprintf( ...
        'Moyenne empirique globale = %.4f', p_merge_emp_mean));

yline(p_merge_th, ...
    '--', ...
    'LineWidth', 2.2, ...
    'DisplayName', sprintf( ...
        'p_{merge}^{th} = %.4f', p_merge_th));

xlabel('Temps au milieu de la transition (s)');
ylabel('p_{merge}');
title('Comparaison temporelle de p_{merge} empirique et theorique');

legend('Location', 'best');

valid_values = [ ...
    p_merge_emp(isfinite(p_merge_emp)); ...
    p_merge_emp_moving(isfinite(p_merge_emp_moving)); ...
    p_merge_emp_mean; ...
    p_merge_th];

if ~isempty(valid_values)
    ymax = max(valid_values);
    ylim([0, min(1, 1.15*max(ymax, eps))]);
end

hold off;

%% ============================================================
%  Affichage des valeurs
%% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' COMPARAISON p_merge EMPIRIQUE / THEORIQUE\n');
fprintf('============================================================\n');
fprintf('p_merge theorique probabiliste       : %.8f\n', ...
    p_merge_th);
fprintf('p_merge empirique global pondere     : %.8f\n', ...
    p_merge_emp_mean);
fprintf('moyenne temporelle empirique simple  : %.8f\n', ...
    p_merge_emp_time_mean);
fprintf('ecart absolu theorie / empirique     : %.8f\n', ...
    abs(p_merge_th - p_merge_emp_mean));

if p_merge_emp_mean > 0
    fprintf('rapport theorie / empirique          : %.8f\n', ...
        p_merge_th / p_merge_emp_mean);
end

fprintf('============================================================\n');
