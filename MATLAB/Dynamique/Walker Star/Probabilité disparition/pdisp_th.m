%% pdisp_depuis_lambda_vrel.m
% Calcule p_merge(t), p_break(t), puis p_disp(t) a partir de :
%   - lambda_effective_temp_liens_inter_results.mat
%   - angle_vitesses_lien_temp_corrige_v2_results.mat
%
% Les calculs sont delegues a trois fonctions separees :
%   calc_p_merge_temp.m
%   calc_p_break_temp.m
%   calc_p_disp_temp.m

clear; clc; close all;

%% Chargement de lambda_eff(t)
script_dir = fileparts(mfilename('fullpath'));
lambda_file = fullfile(script_dir, '..', 'lambda_effective_temp_liens_inter_results.mat');

if ~isfile(lambda_file)
    error(['Fichier %s introuvable. Lance d''abord ', ...
           'lambda_temp.m.'], lambda_file);
end
Slam = load(lambda_file);

required_lambda = {'time_values','lambda_eff_t','lambda_global','dmax'};
for q = 1:numel(required_lambda)
    if ~isfield(Slam, required_lambda{q})
        error('Variable manquante dans %s : %s', lambda_file, required_lambda{q});
    end
end

time_lambda = Slam.time_values(:);
lambda_eff_t = Slam.lambda_eff_t(:);
lambda_global = Slam.lambda_global;
dmax = Slam.dmax;

%% Chargement de v_rel(t)
script_dir = fileparts(mfilename('fullpath'));
vrel_file = fullfile(script_dir, '..', 'Vitesse relative', 'angle_vitesses_lien_temp_corrige_v2_results.mat');

if ~isfile(vrel_file)
    error(['Fichier %s introuvable. Lance d''abord le script ', ...
           'vitesse_rel_temp.m.'], vrel_file);
end
Sv = load(vrel_file);

required_vrel = {'time_values','dt','vrel_link_theory_strates','vrel_link_direct','v_orb'};
for q = 1:numel(required_vrel)
    if ~isfield(Sv, required_vrel{q})
        error('Variable manquante dans %s : %s', vrel_file, required_vrel{q});
    end
end

time_vrel = Sv.time_values(:);
dt = Sv.dt;
v_orb = Sv.v_orb;

% On utilise la vitesse relative theorique par strates pour le modele.
vrel_model_t = Sv.vrel_link_theory_strates(:);
vrel_emp_t = Sv.vrel_link_direct(:);

%% Alignement temporel
% Grille de reference : celle de lambda_eff(t).
time_values = time_lambda;

if numel(time_vrel) ~= numel(time_values) || any(abs(time_vrel - time_values) > 1e-9)
    vrel_model_t = interp1(time_vrel, vrel_model_t, time_values, 'linear', 'extrap');
    vrel_emp_t   = interp1(time_vrel, vrel_emp_t,   time_values, 'linear', 'extrap');
end

%% Calculs par fonctions separees
p_merge_t = calc_p_merge_temp(lambda_eff_t, vrel_model_t, dmax, dt);
p_break_t = calc_p_break_temp(vrel_model_t, dmax, dt);
p_disp_t  = calc_p_disp_temp(p_merge_t, p_break_t);

%% References constantes pour comparaison
p_merge_global_const = calc_p_merge_temp(lambda_global * ones(size(vrel_model_t)), ...
                                        mean(vrel_model_t,'omitnan') * ones(size(vrel_model_t)), ...
                                        dmax, dt);
p_break_const = calc_p_break_temp(mean(vrel_model_t,'omitnan') * ones(size(vrel_model_t)), dmax, dt);
p_disp_const  = calc_p_disp_temp(p_merge_global_const, p_break_const);

p_merge_global_const = p_merge_global_const(1);
p_break_const = p_break_const(1);
p_disp_const = p_disp_const(1);

%% Affichage console
fprintf('\n=== Modele temporel p_merge, p_break, p_disp ===\n');
fprintf('dmax                                      : %.2f km\n', dmax);
fprintf('dt                                        : %.2f s\n', dt);
fprintf('v_rel modele moyen                        : %.6f km/s\n', mean(vrel_model_t,'omitnan'));
fprintf('lambda_eff moyen                          : %.6e sat/km^2\n', mean(lambda_eff_t,'omitnan'));
fprintf('lambda_global                             : %.6e sat/km^2\n', lambda_global);
fprintf('\n');
fprintf('p_merge(t) moyen / min / max              : %.6f / %.6f / %.6f\n', ...
    mean(p_merge_t,'omitnan'), min(p_merge_t,[],'omitnan'), max(p_merge_t,[],'omitnan'));
fprintf('p_break(t) moyen / min / max              : %.6f / %.6f / %.6f\n', ...
    mean(p_break_t,'omitnan'), min(p_break_t,[],'omitnan'), max(p_break_t,[],'omitnan'));
fprintf('p_disp(t) moyen / min / max               : %.6f / %.6f / %.6f\n', ...
    mean(p_disp_t,'omitnan'), min(p_disp_t,[],'omitnan'), max(p_disp_t,[],'omitnan'));
fprintf('\n');
fprintf('Reference constante p_merge global         : %.6f\n', p_merge_global_const);
fprintf('Reference constante p_break mean vrel      : %.6f\n', p_break_const);
fprintf('Reference constante p_disp                 : %.6f\n', p_disp_const);

%% Figures
figure;
plot(time_values, p_merge_t, 'LineWidth', 1.5); hold on;
yline(mean(p_merge_t,'omitnan'), ':', sprintf('moyenne = %.4f', mean(p_merge_t,'omitnan')), ...
    'LabelHorizontalAlignment','left');
yline(p_merge_global_const, '--', sprintf('ref constante = %.4f', p_merge_global_const), ...
    'LabelHorizontalAlignment','left');
grid on;
xlabel('Temps (s)');
ylabel('p_{merge}(t)');
title('Evolution temporelle de p_{merge}(t)');
legend('p_{merge}(t)', 'Moyenne temporelle', 'Reference constante', 'Location','best');

figure;
plot(time_values, p_break_t, 'LineWidth', 1.5); hold on;
yline(mean(p_break_t,'omitnan'), ':', sprintf('moyenne = %.4f', mean(p_break_t,'omitnan')), ...
    'LabelHorizontalAlignment','left');
yline(p_break_const, '--', sprintf('ref constante = %.4f', p_break_const), ...
    'LabelHorizontalAlignment','left');
grid on;
xlabel('Temps (s)');
ylabel('p_{break}(t)');
title('Evolution temporelle de p_{break}(t)');
legend('p_{break}(t)', 'Moyenne temporelle', 'Reference constante', 'Location','best');

figure;
plot(time_values, p_disp_t, 'LineWidth', 1.5); hold on;
yline(mean(p_disp_t,'omitnan'), ':', sprintf('moyenne = %.4f', mean(p_disp_t,'omitnan')), ...
    'LabelHorizontalAlignment','left');
yline(p_disp_const, '--', sprintf('ref constante = %.4f', p_disp_const), ...
    'LabelHorizontalAlignment','left');
grid on;
xlabel('Temps (s)');
ylabel('p_{disp}(t)');
title('Evolution temporelle de p_{disp}(t)');
legend('p_{disp}(t)', 'Moyenne temporelle', 'Reference constante', 'Location','best');

figure;
plot(time_values, p_merge_t, 'LineWidth', 1.3); hold on;
plot(time_values, p_break_t, 'LineWidth', 1.3);
plot(time_values, p_disp_t, 'LineWidth', 1.6);
grid on;
xlabel('Temps (s)');
ylabel('Probabilite par pas');
title('Comparaison p_{merge}(t), p_{break}(t), p_{disp}(t)');
legend('p_{merge}', 'p_{break}', 'p_{disp}', 'Location','best');

%% Sauvegarde
save('pdisp_modele_lambda_vrel_results.mat', ...
    'time_values', 'lambda_eff_t', 'lambda_global', 'vrel_model_t', 'vrel_emp_t', ...
    'dmax', 'dt', 'v_orb', ...
    'p_merge_t', 'p_break_t', 'p_disp_t', ...
    'p_merge_global_const', 'p_break_const', 'p_disp_const');

fprintf('\nResultats sauvegardes dans pdisp_modele_lambda_vrel_results.mat\n');
