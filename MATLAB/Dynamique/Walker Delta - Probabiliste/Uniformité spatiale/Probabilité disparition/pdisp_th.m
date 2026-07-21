%% pdisp_modele_delta_spatial.m
% Modèle temporel de p_merge, p_break et p_disp pour un Walker-Delta
% à uniformité spatiale initiale.
%
% Entrées :
%   - lambda_effective_temp_delta_spatial_results.mat
%   - vitesse_rel_temp_delta_spatial_theorique_strates.mat
%
% Le modèle utilise :
%
%   lambda_eff^Delta(t)
%   v_rel^Delta,th(t)
%
% puis appelle :
%   calc_p_merge_temp.m
%   calc_p_break_temp.m
%   calc_p_disp_temp.m

clear; clc; close all;

%% ============================================================
%  LOCALISATION DES FICHIERS
%% ============================================================

script_dir = fileparts(mfilename('fullpath'));
lambda_file = fullfile(script_dir, '..', 'lambda_effective_temp_delta_spatial_pond_satellites_results.mat');
vrel_file = fullfile(script_dir, '..', 'Vitesse relative', 'vitesse_rel_temp.mat');

if isempty(lambda_file)
    error(['Fichier lambda Delta introuvable. Lance d''abord ', ...
        'lambda_effective_temp_delta_spatial.m.']);
end

if isempty(vrel_file)
    error(['Fichier v_rel Delta introuvable. Lance d''abord ', ...
        'vitesse_rel_temp_delta_spatial_theorique_strates.m.']);
end

fprintf('Fichier lambda : %s\n',lambda_file);
fprintf('Fichier v_rel  : %s\n',vrel_file);

%% ============================================================
%  CHARGEMENT DE lambda_eff(t)
%% ============================================================

Slam = load(lambda_file);

required_lambda = {'time_values','lambda_eff_t','lambda_global_band','dmax'};

for q = 1:numel(required_lambda)
    if ~isfield(Slam,required_lambda{q})
        error('Variable manquante dans %s : %s', ...
            lambda_file,required_lambda{q});
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

%% ============================================================
%  CHARGEMENT DE v_rel(t)
%% ============================================================

Sv = load(vrel_file);

required_vrel = {'time_values','vrel_theory','vrel_emp','v_orb'};

for q = 1:numel(required_vrel)
    if ~isfield(Sv,required_vrel{q})
        error('Variable manquante dans %s : %s', ...
            vrel_file,required_vrel{q});
    end
end

time_vrel = double(Sv.time_values(:));
vrel_model_t = double(Sv.vrel_theory(:));
vrel_emp_t = double(Sv.vrel_emp(:));
v_orb = double(Sv.v_orb);

if isfield(Sv,'dt')
    dt = double(Sv.dt);
elseif numel(time_vrel)>=2
    dt = median(diff(time_vrel));
else
    error('Impossible de déterminer dt.');
end

%% ============================================================
%  ALIGNEMENT TEMPOREL
%% ============================================================

time_values = time_lambda;

if numel(time_vrel) ~= numel(time_values) || ...
        any(abs(time_vrel-time_values)>1e-9)

    vrel_model_t = interp1(time_vrel,vrel_model_t, ...
        time_values,'linear','extrap');

    vrel_emp_t = interp1(time_vrel,vrel_emp_t, ...
        time_values,'linear','extrap');
end

% Les probabilités sont définies sur les transitions [t_k,t_{k+1}].
% On associe à chaque transition les valeurs au début de l'intervalle.
if numel(time_values)>=2
    t_transition = time_values(1:end-1);
    lambda_transition = lambda_eff_t(1:end-1);
    vrel_model_transition = vrel_model_t(1:end-1);
    vrel_emp_transition = vrel_emp_t(1:end-1);
else
    error('La grille temporelle doit contenir au moins deux instants.');
end

%% ============================================================
%  CALCUL THEORIQUE
%% ============================================================

[p_merge_t,chi_merge_t,E_t_merge,beta0_geom_merge] = ...
    calc_p_merge_temp( ...
        lambda_transition, ...
        vrel_model_transition, ...
        dmax,dt);

p_break_t = calc_p_break_temp( ...
    vrel_model_transition,dmax,dt);

p_disp_t = calc_p_disp_temp(p_merge_t,p_break_t);

%% Référence utilisant v_rel empirique, pour diagnostic uniquement
[p_merge_vrel_emp_t,~,~,~] = ...
    calc_p_merge_temp( ...
        lambda_transition, ...
        vrel_emp_transition, ...
        dmax,dt);

p_break_vrel_emp_t = calc_p_break_temp( ...
    vrel_emp_transition,dmax,dt);

p_disp_vrel_emp_t = calc_p_disp_temp( ...
    p_merge_vrel_emp_t,p_break_vrel_emp_t);

%% ============================================================
%  REFERENCES CONSTANTES
%% ============================================================

vrel_mean = mean(vrel_model_transition,'omitnan');

[p_merge_const_t,chi_merge_const_t] = ...
    calc_p_merge_temp( ...
        lambda_global*ones(size(vrel_model_transition)), ...
        vrel_mean*ones(size(vrel_model_transition)), ...
        dmax,dt);

p_break_const_t = calc_p_break_temp( ...
    vrel_mean*ones(size(vrel_model_transition)), ...
    dmax,dt);

p_disp_const_t = calc_p_disp_temp( ...
    p_merge_const_t,p_break_const_t);

p_merge_const = p_merge_const_t(1);
chi_merge_const = chi_merge_const_t(1);
p_break_const = p_break_const_t(1);
p_disp_const = p_disp_const_t(1);

%% ============================================================
%  AFFICHAGE CONSOLE
%% ============================================================

fprintf('\n=== Modèle temporel Walker-Delta spatial ===\n');
fprintf('Inclinaison                                 : %.2f deg\n',inc_deg);
fprintf('dmax                                       : %.2f km\n',dmax);
fprintf('dt                                         : %.2f s\n',dt);
fprintf('v_rel théorique moyen                      : %.6f km/s\n', ...
    mean(vrel_model_transition,'omitnan'));
fprintf('lambda_eff moyen                           : %.6e sat/km^2\n', ...
    mean(lambda_transition,'omitnan'));
fprintf('lambda globale dans la bande               : %.6e sat/km^2\n', ...
    lambda_global);

fprintf('\np_merge moyen / min / max : %.6f / %.6f / %.6f\n', ...
    mean(p_merge_t,'omitnan'), ...
    min(p_merge_t,[],'omitnan'), ...
    max(p_merge_t,[],'omitnan'));

fprintf('p_break moyen / min / max : %.6f / %.6f / %.6f\n', ...
    mean(p_break_t,'omitnan'), ...
    min(p_break_t,[],'omitnan'), ...
    max(p_break_t,[],'omitnan'));

fprintf('p_disp moyen / min / max  : %.6f / %.6f / %.6f\n', ...
    mean(p_disp_t,'omitnan'), ...
    min(p_disp_t,[],'omitnan'), ...
    max(p_disp_t,[],'omitnan'));

fprintf('\nRéférence constante p_disp : %.6f\n',p_disp_const);

%% ============================================================
%  FIGURES
%% ============================================================

figure;
hold on; grid on;
plot(t_transition,p_merge_t,'LineWidth',1.5, ...
    'DisplayName','p_{merge}^{th,\Delta}(t)');
yline(mean(p_merge_t,'omitnan'),':', ...
    sprintf('Moyenne = %.4f',mean(p_merge_t,'omitnan')));
yline(p_merge_const,'--', ...
    sprintf('Référence constante = %.4f',p_merge_const));
xlabel('Temps (s)');
ylabel('p_{merge}(t)');
title('Walker-Delta spatial : p_{merge}^{th}(t)');
legend('Location','best');

figure;
hold on; grid on;
plot(t_transition,p_break_t,'LineWidth',1.5, ...
    'DisplayName','p_{break}^{th,\Delta}(t)');
yline(mean(p_break_t,'omitnan'),':', ...
    sprintf('Moyenne = %.4f',mean(p_break_t,'omitnan')));
yline(p_break_const,'--', ...
    sprintf('Référence constante = %.4f',p_break_const));
xlabel('Temps (s)');
ylabel('p_{break}(t)');
title('Walker-Delta spatial : p_{break}^{th}(t)');
legend('Location','best');

figure;
hold on; grid on;
plot(t_transition,p_disp_t,'LineWidth',1.7, ...
    'DisplayName','Modèle entièrement théorique');
plot(t_transition,p_disp_vrel_emp_t,'--','LineWidth',1.2, ...
    'DisplayName','Diagnostic avec v_{rel} empirique');
yline(mean(p_disp_t,'omitnan'),':', ...
    sprintf('Moyenne = %.4f',mean(p_disp_t,'omitnan')));
yline(p_disp_const,'--', ...
    sprintf('Référence constante = %.4f',p_disp_const));
xlabel('Temps (s)');
ylabel('p_{disp}(t)');
title('Walker-Delta spatial : p_{disp}^{th}(t)');
legend('Location','best');

figure;
hold on; grid on;
plot(t_transition,p_merge_t,'LineWidth',1.3);
plot(t_transition,p_break_t,'LineWidth',1.3);
plot(t_transition,p_disp_t,'LineWidth',1.7);
xlabel('Temps (s)');
ylabel('Probabilité par pas');
title('Décomposition du modèle Delta spatial');
legend('p_{merge}^{th}','p_{break}^{th}','p_{disp}^{th}', ...
    'Location','best');

%% ============================================================
%  SAUVEGARDE
%% ============================================================

save('pdisp_modele_delta_spatial_results.mat', ...
    'time_values','t_transition', ...
    'lambda_eff_t','lambda_transition','lambda_global', ...
    'vrel_model_t','vrel_emp_t', ...
    'vrel_model_transition','vrel_emp_transition', ...
    'dmax','dt','v_orb','inc_deg', ...
    'p_merge_t','p_break_t','p_disp_t', ...
    'p_merge_vrel_emp_t','p_break_vrel_emp_t','p_disp_vrel_emp_t', ...
    'chi_merge_t','E_t_merge','beta0_geom_merge', ...
    'p_merge_const','chi_merge_const','p_break_const','p_disp_const', ...
    'lambda_file','vrel_file');

fprintf('\nRésultats sauvegardés dans ');
fprintf('pdisp_modele_delta_spatial_results.mat\n');

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
