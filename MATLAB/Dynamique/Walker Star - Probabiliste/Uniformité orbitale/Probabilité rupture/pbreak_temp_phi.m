clear; clc; close all;

%% ============================================================
%  COMPARAISON TEMPORELLE DE p_break
%
%  Courbe empirique :
%      chargee depuis pbreak_emp_results.mat
%
%  Moyenne theorique :
%      chargee depuis pbreak_phi_th_results.mat
%
%  Sortie :
%      figure avec p_break empirique, moyenne glissante,
%      moyenne empirique globale et valeur theorique constante.
%% ============================================================
emp_file = 'pbreak_emp_results.mat';
th_file  = 'pbreak_phi_th_results.mat';

if ~isfile(emp_file)
    error('Fichier empirique introuvable : %s', emp_file);
end

if ~isfile(th_file)
    error('Fichier theorique introuvable : %s', th_file);
end

%% Chargement des resultats empiriques
Semp = load(emp_file, ...
    'time_transition', ...
    'p_break', ...
    'p_break_moving', ...
    'p_break_mean', ...
    'p_break_time_mean', ...
    'moving_window');

time_transition = Semp.time_transition(:);
p_break_emp = Semp.p_break(:);
p_break_emp_moving = Semp.p_break_moving(:);
p_break_emp_mean = Semp.p_break_mean;
p_break_emp_time_mean = Semp.p_break_time_mean;
moving_window = Semp.moving_window;

%% Chargement des resultats theoriques
Sth = load(th_file);

% Valeur theorique globale apres deconditionnement local puis
% moyenne selon la loi de latitude des composantes.
if isfield(Sth,'p_break_th')
    p_break_th_deconditioned = Sth.p_break_th;
else
    error(['Le fichier %s doit contenir p_break_th. ', ...
           'Relancez d''abord pbreak_phi_th.m.'],th_file);
end

% Version utilisant beta0(phi) et p_bridge,bord(phi) empiriques.
if isfield(Sth,'p_break_th_corrected')
    p_break_th_corrected = Sth.p_break_th_corrected;
else
    p_break_th_corrected = NaN;
    warning(['Le fichier %s ne contient pas ', ...
             'p_break_th_corrected.'],th_file);
end

% Versions lineaires, uniquement pour l'affichage console.
if isfield(Sth,'p_break_th_linear')
    p_break_th_linear = Sth.p_break_th_linear;
else
    p_break_th_linear = NaN;
end

if isfield(Sth,'p_break_th_linear_corrected')
    p_break_th_linear_corrected = ...
        Sth.p_break_th_linear_corrected;
else
    p_break_th_linear_corrected = NaN;
end

%% Verifications
if numel(time_transition) ~= numel(p_break_emp)
    error(['Les dimensions de time_transition et p_break ', ...
           'ne correspondent pas.']);
end

if numel(p_break_emp_moving) ~= numel(p_break_emp)
    error(['Les dimensions de p_break_moving et p_break ', ...
           'ne correspondent pas.']);
end

%% ============================================================
%  Trace principal
%% ============================================================

figure;
hold on;
grid on;
box on;

plot(time_transition, p_break_emp, ...
    '-', ...
    'LineWidth', 0.8, ...
    'DisplayName', 'p_{break}^{emp}(t)');

plot(time_transition, p_break_emp_moving, ...
    'LineWidth', 2.2, ...
    'DisplayName', sprintf( ...
        'Moyenne glissante empirique (%d points)', moving_window));

yline(p_break_emp_mean, ...
    ':', ...
    'LineWidth', 1.8, ...
    'DisplayName', sprintf( ...
        'Moyenne empirique globale = %.4f', p_break_emp_mean));

yline(p_break_th_deconditioned, ...
    '--', ...
    'LineWidth', 2.2, ...
    'DisplayName', sprintf( ...
        'p_{break Delta,th} deconditionne = %.4f', ...
        p_break_th_deconditioned));

if isfinite(p_break_th_corrected)
    yline(p_break_th_corrected, ...
        '-.', ...
        'LineWidth', 2.2, ...
        'DisplayName', sprintf( ...
            'p_{break Delta,th,corr} = %.4f', ...
            p_break_th_corrected));
end

xlabel('Temps au milieu de la transition (s)');
ylabel('p_{break}');
title('Comparaison temporelle de p_{break} -- modele Delta orbital');

legend('Location', 'best');

valid_values = [ ...
    p_break_emp(isfinite(p_break_emp)); ...
    p_break_emp_moving(isfinite(p_break_emp_moving)); ...
    p_break_emp_mean; ...
    p_break_th_deconditioned; ...
    p_break_th_corrected];

valid_values = valid_values(isfinite(valid_values));

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
fprintf(' COMPARAISON p_break EMPIRIQUE / THEORIQUE\n');
fprintf('============================================================\n');
fprintf('p_break theorique deconditionne       : %.8f\n', ...
    p_break_th_deconditioned);

if isfinite(p_break_th_corrected)
    fprintf('p_break theorique corrige empiriquement : %.8f\n', ...
        p_break_th_corrected);
end

if isfinite(p_break_th_linear)
    fprintf('p_break lineaire deconditionne        : %.8f\n', ...
        p_break_th_linear);
end

if isfinite(p_break_th_linear_corrected)
    fprintf('p_break lineaire corrige              : %.8f\n', ...
        p_break_th_linear_corrected);
end

fprintf('p_break empirique global pondere      : %.8f\n', ...
    p_break_emp_mean);
fprintf('moyenne temporelle empirique simple   : %.8f\n', ...
    p_break_emp_time_mean);

fprintf('ecart absolu deconditionne / empirique: %.8f\n', ...
    abs(p_break_th_deconditioned-p_break_emp_mean));

if isfinite(p_break_th_corrected)
    fprintf('ecart absolu corrige / empirique      : %.8f\n', ...
        abs(p_break_th_corrected-p_break_emp_mean));
end

if p_break_emp_mean > 0
    fprintf('rapport deconditionne / empirique     : %.8f\n', ...
        p_break_th_deconditioned/p_break_emp_mean);

    if isfinite(p_break_th_corrected)
        fprintf('rapport corrige / empirique           : %.8f\n', ...
            p_break_th_corrected/p_break_emp_mean);
    end
end

fprintf('============================================================\n');
