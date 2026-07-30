%% pmerge_th.m
% Probabilite theorique temporelle de fusion pour un Walker-Delta
% a uniformite spatiale initiale.
%
% Entrees :
%   - ../lambda_eff_th_results.mat
%   - ../Vitesse relative/vitesse_rel_temp.mat
%
% Sortie :
%   - pmerge_th_results.mat

clear; clc; close all;

%% Localisation des fichiers
script_dir = fileparts(mfilename('fullpath'));
lambda_file = fullfile(script_dir, '..', 'Paramètres', 'lambda_eff_th_results.mat');
vrel_file = fullfile(script_dir, '..', 'Paramètres', 'Vitesse relative', 'vitesse_rel_temp.mat');

if ~isfile(lambda_file)
    error('Fichier lambda introuvable : %s', lambda_file);
end
if ~isfile(vrel_file)
    error('Fichier de vitesse relative introuvable : %s', vrel_file);
end

fprintf('Fichier lambda : %s\n', lambda_file);
fprintf('Fichier v_rel  : %s\n', vrel_file);

%% Chargement de lambda_eff(t)
Slam = load(lambda_file);
required_lambda = {'time_values','lambda_eff_t','lambda_global_band','dmax'};
for q = 1:numel(required_lambda)
    if ~isfield(Slam, required_lambda{q})
        error('Variable manquante dans %s : %s', lambda_file, required_lambda{q});
    end
end

time_lambda = double(Slam.time_values(:));
lambda_eff_t = double(Slam.lambda_eff_t(:));
lambda_global = double(Slam.lambda_global_band);
dmax = double(Slam.dmax);

if isfield(Slam,'inc_deg')
    inc_deg = double(Slam.inc_deg);
else
    inc_deg = NaN;
end

%% Chargement de v_rel(t)
Sv = load(vrel_file);
required_vrel = {'time_values','vrel_theory','vrel_emp','v_orb'};
for q = 1:numel(required_vrel)
    if ~isfield(Sv, required_vrel{q})
        error('Variable manquante dans %s : %s', vrel_file, required_vrel{q});
    end
end

time_vrel = double(Sv.time_values(:));
vrel_theory_t = double(Sv.vrel_theory(:));
vrel_emp_t = double(Sv.vrel_emp(:));
v_orb = double(Sv.v_orb);

if isfield(Sv,'dt')
    dt = double(Sv.dt);
elseif numel(time_vrel) >= 2
    dt = median(diff(time_vrel));
else
    error('Impossible de determiner dt.');
end

%% Alignement temporel
time_values = time_lambda;
if numel(time_vrel) ~= numel(time_values) || any(abs(time_vrel-time_values) > 1e-9)
    vrel_theory_t = interp1(time_vrel, vrel_theory_t, time_values, 'linear', 'extrap');
    vrel_emp_t = interp1(time_vrel, vrel_emp_t, time_values, 'linear', 'extrap');
end

if numel(time_values) < 2
    error('La grille temporelle doit contenir au moins deux instants.');
end

t_transition = time_values(1:end-1);
lambda_transition = lambda_eff_t(1:end-1);
vrel_theory_transition = vrel_theory_t(1:end-1);
vrel_emp_transition = vrel_emp_t(1:end-1);

%% Calcul theorique
[p_merge_t, chi_merge_t, E_t_merge, beta0_geom_merge] = ...
    calc_p_merge_th(lambda_transition, vrel_theory_transition, dmax, dt);

%% Diagnostic avec vitesse empirique uniquement
[p_merge_vrel_emp_t, chi_merge_vrel_emp_t] = ...
    calc_p_merge_th(lambda_transition, vrel_emp_transition, dmax, dt);

%% Reference constante
vrel_mean = mean(vrel_theory_transition, 'omitnan');
[p_merge_const_t, chi_merge_const_t] = ...
    calc_p_merge_th(lambda_global * ones(size(vrel_theory_transition)), ...
                      vrel_mean * ones(size(vrel_theory_transition)), ...
                      dmax, dt);

p_merge_const = p_merge_const_t(1);
chi_merge_const = chi_merge_const_t(1);

%% Affichage
fprintf('\n=== p_merge theorique Walker-Delta spatial ===\n');
fprintf('Inclinaison                    : %.2f deg\n', inc_deg);
fprintf('dmax                           : %.2f km\n', dmax);
fprintf('dt                             : %.2f s\n', dt);
fprintf('v_rel theorique moyen          : %.6f km/s\n', mean(vrel_theory_transition,'omitnan'));
fprintf('lambda_eff moyen               : %.6e sat/km^2\n', mean(lambda_transition,'omitnan'));
fprintf('p_merge moyen / min / max      : %.6f / %.6f / %.6f\n', ...
    mean(p_merge_t,'omitnan'), min(p_merge_t,[],'omitnan'), max(p_merge_t,[],'omitnan'));
fprintf('Reference constante p_merge    : %.6f\n', p_merge_const);

figure;
hold on; grid on;
plot(t_transition, p_merge_t, 'LineWidth', 1.6, ...
    'DisplayName', 'p_{merge}^{th,\Delta}(t)');
plot(t_transition, p_merge_vrel_emp_t, '--', 'LineWidth', 1.2, ...
    'DisplayName', 'Diagnostic avec v_{rel} empirique');
yline(mean(p_merge_t,'omitnan'), ':', ...
    sprintf('Moyenne = %.4f', mean(p_merge_t,'omitnan')));
yline(p_merge_const, '--', ...
    sprintf('Reference constante = %.4f', p_merge_const));
xlabel('Temps (s)');
ylabel('p_{merge}(t)');
title('Walker-Delta spatial : p_{merge}^{th}(t)');
legend('Location','best');

%% Sauvegarde
save('pmerge_th_results.mat', ...
    'time_values','t_transition', ...
    'lambda_eff_t','lambda_transition','lambda_global', ...
    'vrel_theory_t','vrel_emp_t', ...
    'vrel_theory_transition','vrel_emp_transition', ...
    'dmax','dt','v_orb','inc_deg', ...
    'p_merge_t','p_merge_vrel_emp_t', ...
    'chi_merge_t','chi_merge_vrel_emp_t', ...
    'E_t_merge','beta0_geom_merge', ...
    'p_merge_const','chi_merge_const', ...
    'lambda_file','vrel_file');

fprintf('\nResultats sauvegardes dans pmerge_th_results.mat\n');
