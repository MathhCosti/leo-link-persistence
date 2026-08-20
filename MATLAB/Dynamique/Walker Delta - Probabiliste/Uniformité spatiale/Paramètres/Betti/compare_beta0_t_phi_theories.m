%% compare_beta0_t_phi_theories.m
% Compare les deux modeles theoriques de beta_0(t,phi) avec le meme
% empirique complet.
%
% Modele 1 : ancienne theorie
%   beta0_old(t,phi) = N1 + N2 + N3 + contribution macro
%   -> betti_t_phi_results.mat
%
% Modele 2 : nouvelle theorie PPP par expansion en tailles de composantes
%   beta0_new(t,phi) = sum_s C_s(t,phi)
%   -> beta0_t_phi_results.mat
%
% Empirique :
%   beta0_emp(t,phi) compte toutes les composantes connexes.
%
% Le script produit :
%   - heatmaps des deux theories et de l'empirique ;
%   - heatmaps des deux erreurs ;
%   - coupes en latitude a plusieurs instants ;
%   - comparaison de beta0 total(t) ;
%   - moyenne temporelle en fonction de phi ;
%   - RMSE/MAE/biais selon la latitude ;
%   - comparaison globale des erreurs.

clear; clc; close all;

%% ============================================================
%  1. CHARGEMENT
%% ============================================================

script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end

old_candidates = {
    fullfile(script_dir,'betti_t_phi_results.mat')
    fullfile(script_dir,'betti_t_phi_results(1).mat')
    fullfile(script_dir,'..','betti_t_phi_results.mat')
};

new_candidates = {
    fullfile(script_dir,'beta0_t_phi_results.mat')
    fullfile(script_dir,'..','beta0_t_phi_results.mat')
};

old_file = find_first_file(old_candidates,'betti_t_phi_results.mat');
new_file = find_first_file(new_candidates,'beta0_t_phi_results.mat');

Sold = load(old_file);
Snew = load(new_file);

required_old = { ...
    'N','R','dmax','inc', ...
    'time_values','phi_vals','phi_edges', ...
    'beta0_th','beta0_emp_true','beta0_emp_true_sem', ...
    'beta0_total_th','beta0_total_emp_true', ...
    'beta0_total_emp_true_sem'};

required_new = { ...
    'N','R','dmax','inc', ...
    'time_values','phi_vals','phi_edges', ...
    'beta0_bin_th','beta0_bin_emp','beta0_bin_emp_sem', ...
    'beta0_total_th_t','beta0_total_emp_t', ...
    'beta0_total_emp_sem'};

require_fields(Sold,required_old,old_file);
require_fields(Snew,required_new,new_file);

%% ============================================================
%  2. VERIFICATION DES PARAMETRES ET DES GRILLES
%% ============================================================

N_old = double(Sold.N);
N_new = double(Snew.N);

R_old = double(Sold.R);
R_new = double(Snew.R);

dmax_old = double(Sold.dmax);
dmax_new = double(Snew.dmax);

inc_old = double(Sold.inc);
inc_new = double(Snew.inc);

assert_close(N_old,N_new,1e-10,'N');
assert_close(R_old,R_new,1e-10,'R');
assert_close(dmax_old,dmax_new,1e-10,'dmax');
assert_close(inc_old,inc_new,1e-10,'inc');

time_old = col_vector(Sold.time_values);
time_new = col_vector(Snew.time_values);

phi_old = row_vector(Sold.phi_vals);
phi_new = row_vector(Snew.phi_vals);

if numel(time_old) ~= numel(time_new) || ...
        max(abs(time_old-time_new)) > 1e-9
    error('Les grilles temporelles des deux modeles sont incompatibles.');
end

if numel(phi_old) ~= numel(phi_new) || ...
        max(abs(phi_old-phi_new)) > 1e-12
    error('Les grilles de latitude des deux modeles sont incompatibles.');
end

time_values = time_new;
phi_vals = phi_new;

Nt = numel(time_values);
Nb = numel(phi_vals);

%% ============================================================
%  3. RECUPERATION DES TROIS COURBES
%% ============================================================

beta0_old = double(Sold.beta0_th);
beta0_new = double(Snew.beta0_bin_th);

% Les deux fichiers ont ete produits avec la meme procedure empirique.
% On verifie qu'ils contiennent bien le meme resultat.
beta0_emp_old_file = double(Sold.beta0_emp_true);
beta0_emp_new_file = double(Snew.beta0_bin_emp);

check_size(beta0_old,Nt,Nb,'beta0_old');
check_size(beta0_new,Nt,Nb,'beta0_new');
check_size(beta0_emp_old_file,Nt,Nb,'beta0_emp_old_file');
check_size(beta0_emp_new_file,Nt,Nb,'beta0_emp_new_file');

max_empirical_difference = max(abs( ...
    beta0_emp_old_file(:)-beta0_emp_new_file(:)));

if max_empirical_difference > 1e-10
    warning(['Les deux fichiers ne contiennent pas exactement la meme ', ...
             'mesure empirique. La comparaison utilisera celle de ', ...
             'beta0_t_phi_results.mat. Ecart max = %.3e.'], ...
             max_empirical_difference);
end

beta0_emp = beta0_emp_new_file;
beta0_emp_sem = double(Snew.beta0_bin_emp_sem);

%% ============================================================
%  4. ERREURS LOCALES
%% ============================================================

error_old = beta0_old-beta0_emp;
error_new = beta0_new-beta0_emp;

valid_old = isfinite(error_old);
valid_new = isfinite(error_new);

rmse_old = sqrt(mean(error_old(valid_old).^2));
rmse_new = sqrt(mean(error_new(valid_new).^2));

mae_old = mean(abs(error_old(valid_old)));
mae_new = mean(abs(error_new(valid_new)));

bias_old = mean(error_old(valid_old));
bias_new = mean(error_new(valid_new));

% Erreurs par latitude.
rmse_old_by_phi = nan(1,Nb);
rmse_new_by_phi = nan(1,Nb);

mae_old_by_phi = nan(1,Nb);
mae_new_by_phi = nan(1,Nb);

bias_old_by_phi = nan(1,Nb);
bias_new_by_phi = nan(1,Nb);

for b = 1:Nb

    eo = error_old(:,b);
    en = error_new(:,b);

    vo = isfinite(eo);
    vn = isfinite(en);

    if any(vo)
        rmse_old_by_phi(b) = sqrt(mean(eo(vo).^2));
        mae_old_by_phi(b) = mean(abs(eo(vo)));
        bias_old_by_phi(b) = mean(eo(vo));
    end

    if any(vn)
        rmse_new_by_phi(b) = sqrt(mean(en(vn).^2));
        mae_new_by_phi(b) = mean(abs(en(vn)));
        bias_new_by_phi(b) = mean(en(vn));
    end
end

%% ============================================================
%  5. beta0 GLOBAL(t)
%% ============================================================

beta0_total_old = col_vector(Sold.beta0_total_th);
beta0_total_new = col_vector(Snew.beta0_total_th_t);
beta0_total_emp = col_vector(Snew.beta0_total_emp_t);
beta0_total_emp_sem = col_vector(Snew.beta0_total_emp_sem);

if numel(beta0_total_old) ~= Nt || ...
        numel(beta0_total_new) ~= Nt || ...
        numel(beta0_total_emp) ~= Nt
    error('Dimensions incompatibles pour beta0 total(t).');
end

error_total_old = beta0_total_old-beta0_total_emp;
error_total_new = beta0_total_new-beta0_total_emp;

rmse_total_old = sqrt(mean(error_total_old.^2,'omitnan'));
rmse_total_new = sqrt(mean(error_total_new.^2,'omitnan'));

mae_total_old = mean(abs(error_total_old),'omitnan');
mae_total_new = mean(abs(error_total_new),'omitnan');

bias_total_old = mean(error_total_old,'omitnan');
bias_total_new = mean(error_total_new,'omitnan');

mean_total_old = mean(beta0_total_old,'omitnan');
mean_total_new = mean(beta0_total_new,'omitnan');
mean_total_emp = mean(beta0_total_emp,'omitnan');

relative_error_mean_old = ...
    abs(mean_total_old-mean_total_emp)/max(abs(mean_total_emp),eps);

relative_error_mean_new = ...
    abs(mean_total_new-mean_total_emp)/max(abs(mean_total_emp),eps);

%% ============================================================
%  6. AFFICHAGE CONSOLE
%% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' COMPARAISON beta0(t,phi) : DEUX THEORIES / EMPIRIQUE\n');
fprintf('============================================================\n');
fprintf('N                                   : %.0f\n',N_new);
fprintf('R                                   : %.6f km\n',R_new);
fprintf('Inclinaison                         : %.6f deg\n',rad2deg(inc_new));
fprintf('dmax                                : %.6f km\n',dmax_new);
fprintf('Nt / Nb                             : %d / %d\n',Nt,Nb);
fprintf('Ecart max entre empiriques fichiers : %.3e\n', ...
    max_empirical_difference);
fprintf('------------------------------------------------------------\n');
fprintf('RMSE locale ancienne theorie        : %.8f\n',rmse_old);
fprintf('RMSE locale nouvelle theorie PPP    : %.8f\n',rmse_new);
fprintf('MAE locale ancienne theorie         : %.8f\n',mae_old);
fprintf('MAE locale nouvelle theorie PPP     : %.8f\n',mae_new);
fprintf('Biais local ancienne theorie        : %+.8f\n',bias_old);
fprintf('Biais local nouvelle theorie PPP    : %+.8f\n',bias_new);
fprintf('------------------------------------------------------------\n');
fprintf('beta0 moyen empirique               : %.8f\n',mean_total_emp);
fprintf('beta0 moyen ancienne theorie        : %.8f\n',mean_total_old);
fprintf('beta0 moyen nouvelle theorie PPP    : %.8f\n',mean_total_new);
fprintf('Erreur relative moyenne ancienne    : %.3f %%\n', ...
    100*relative_error_mean_old);
fprintf('Erreur relative moyenne nouvelle    : %.3f %%\n', ...
    100*relative_error_mean_new);
fprintf('------------------------------------------------------------\n');
fprintf('RMSE beta0 total ancienne           : %.8f\n',rmse_total_old);
fprintf('RMSE beta0 total nouvelle           : %.8f\n',rmse_total_new);
fprintf('MAE beta0 total ancienne            : %.8f\n',mae_total_old);
fprintf('MAE beta0 total nouvelle            : %.8f\n',mae_total_new);
fprintf('============================================================\n');

%% ============================================================
%  7. HEATMAPS
%% ============================================================

% Empirique
figure;
imagesc(time_values,rad2deg(phi_vals),beta0_emp.');
axis xy;
colorbar;
xlabel('Temps (s)');
ylabel('Latitude \phi (deg)');
title('\beta_0^{emp}(t,\phi) : toutes les composantes');

% Ancienne theorie
figure;
imagesc(time_values,rad2deg(phi_vals),beta0_old.');
axis xy;
colorbar;
xlabel('Temps (s)');
ylabel('Latitude \phi (deg)');
title('\beta_0^{th}(t,\phi) : ancienne theorie N_1+N_2+N_3+macro');

% Nouvelle theorie
figure;
imagesc(time_values,rad2deg(phi_vals),beta0_new.');
axis xy;
colorbar;
xlabel('Temps (s)');
ylabel('Latitude \phi (deg)');
title('\beta_0^{th}(t,\phi) : nouvelle theorie PPP \Sigma_s C_s');

% Erreur ancienne
figure;
imagesc(time_values,rad2deg(phi_vals),error_old.');
axis xy;
colorbar;
xlabel('Temps (s)');
ylabel('Latitude \phi (deg)');
title('Erreur ancienne theorie : \beta_0^{th}-\beta_0^{emp}');

% Erreur nouvelle
figure;
imagesc(time_values,rad2deg(phi_vals),error_new.');
axis xy;
colorbar;
xlabel('Temps (s)');
ylabel('Latitude \phi (deg)');
title('Erreur nouvelle theorie PPP : \beta_0^{th}-\beta_0^{emp}');

%% ============================================================
%  8. COUPES EN LATITUDE A PLUSIEURS INSTANTS
%% ============================================================

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

    plot(rad2deg(phi_vals),beta0_old(it,:), ...
        'LineWidth',1.8, ...
        'DisplayName','Ancienne theorie');

    plot(rad2deg(phi_vals),beta0_new(it,:), ...
        'LineWidth',2.0, ...
        'DisplayName','Nouvelle theorie PPP');

    errorbar(rad2deg(phi_vals), ...
        beta0_emp(it,:), ...
        beta0_emp_sem(it,:), ...
        'o--','LineWidth',1.1, ...
        'DisplayName','Empirique \pm SEM');

    grid on;
    ylabel('\beta_0 par tranche');
    title(sprintf('t = %.1f s',time_values(it)));

    if k == 1
        legend('Location','best');
    end

    if k == numel(selected_indices)
        xlabel('Latitude \phi (deg)');
    end

    hold off;
end

%% ============================================================
%  9. COMPARAISON GLOBALE beta0(t)
%% ============================================================

figure;
hold on;

plot(time_values,beta0_total_old, ...
    'LineWidth',1.8, ...
    'DisplayName','Ancienne theorie');

plot(time_values,beta0_total_new, ...
    'LineWidth',2.0, ...
    'DisplayName','Nouvelle theorie PPP');

plot(time_values,beta0_total_emp, ...
    '--','LineWidth',1.5, ...
    'DisplayName','Empirique complet');

plot(time_values, ...
    beta0_total_emp+beta0_total_emp_sem, ...
    ':','LineWidth',0.9, ...
    'DisplayName','Empirique + SEM');

plot(time_values, ...
    beta0_total_emp-beta0_total_emp_sem, ...
    ':','LineWidth',0.9, ...
    'DisplayName','Empirique - SEM');

grid on;
xlabel('Temps (s)');
ylabel('\beta_0(t)');
title('\beta_0(t) : comparaison des deux theories avec l''empirique');
legend('Location','best');
hold off;

%% ============================================================
%  10. MOYENNE TEMPORELLE SELON LA LATITUDE
%% ============================================================

mean_old_phi = mean(beta0_old,1,'omitnan');
mean_new_phi = mean(beta0_new,1,'omitnan');
mean_emp_phi = mean(beta0_emp,1,'omitnan');

sem_emp_phi_time = ...
    std(beta0_emp,0,1,'omitnan')/sqrt(Nt);

figure;
hold on;

plot(rad2deg(phi_vals),mean_old_phi, ...
    'LineWidth',1.8, ...
    'DisplayName','Ancienne theorie');

plot(rad2deg(phi_vals),mean_new_phi, ...
    'LineWidth',2.0, ...
    'DisplayName','Nouvelle theorie PPP');

errorbar(rad2deg(phi_vals), ...
    mean_emp_phi,sem_emp_phi_time, ...
    'o--','LineWidth',1.1, ...
    'DisplayName','Empirique moyen temporel');

grid on;
xlabel('Latitude \phi (deg)');
ylabel('\beta_0 moyen par tranche');
title('Moyenne temporelle de \beta_0(t,\phi)');
legend('Location','best');
hold off;

%% ============================================================
%  11. ERREUR SELON LA LATITUDE
%% ============================================================

figure;
hold on;

plot(rad2deg(phi_vals),rmse_old_by_phi, ...
    'LineWidth',1.8, ...
    'DisplayName','RMSE ancienne theorie');

plot(rad2deg(phi_vals),rmse_new_by_phi, ...
    'LineWidth',2.0, ...
    'DisplayName','RMSE nouvelle theorie PPP');

grid on;
xlabel('Latitude \phi (deg)');
ylabel('RMSE sur \beta_0 par tranche');
title('Erreur locale selon la latitude');
legend('Location','best');
hold off;

figure;
hold on;

plot(rad2deg(phi_vals),bias_old_by_phi, ...
    'LineWidth',1.8, ...
    'DisplayName','Biais ancienne theorie');

plot(rad2deg(phi_vals),bias_new_by_phi, ...
    'LineWidth',2.0, ...
    'DisplayName','Biais nouvelle theorie PPP');

yline(0,'--','Accord parfait');

grid on;
xlabel('Latitude \phi (deg)');
ylabel('Biais moyen th-emp');
title('Biais des deux theories selon la latitude');
legend('Location','best');
hold off;

%% ============================================================
%  12. ERREUR EN FONCTION DU TEMPS
%% ============================================================

rmse_old_by_t = sqrt(mean(error_old.^2,2,'omitnan'));
rmse_new_by_t = sqrt(mean(error_new.^2,2,'omitnan'));

figure;
hold on;

plot(time_values,rmse_old_by_t, ...
    'LineWidth',1.8, ...
    'DisplayName','Ancienne theorie');

plot(time_values,rmse_new_by_t, ...
    'LineWidth',2.0, ...
    'DisplayName','Nouvelle theorie PPP');

grid on;
xlabel('Temps (s)');
ylabel('RMSE sur les tranches de latitude');
title('Erreur locale des deux theories en fonction du temps');
legend('Location','best');
hold off;

%% ============================================================
%  13. SAUVEGARDE
%% ============================================================

output_file = fullfile(script_dir, ...
    'compare_beta0_t_phi_theories_results.mat');

save(output_file, ...
    'old_file','new_file', ...
    'N_new','R_new','dmax_new','inc_new', ...
    'time_values','phi_vals','Nt','Nb', ...
    'beta0_old','beta0_new','beta0_emp','beta0_emp_sem', ...
    'error_old','error_new', ...
    'rmse_old','rmse_new','mae_old','mae_new', ...
    'bias_old','bias_new', ...
    'rmse_old_by_phi','rmse_new_by_phi', ...
    'mae_old_by_phi','mae_new_by_phi', ...
    'bias_old_by_phi','bias_new_by_phi', ...
    'rmse_old_by_t','rmse_new_by_t', ...
    'beta0_total_old','beta0_total_new','beta0_total_emp', ...
    'beta0_total_emp_sem', ...
    'error_total_old','error_total_new', ...
    'rmse_total_old','rmse_total_new', ...
    'mae_total_old','mae_total_new', ...
    'bias_total_old','bias_total_new', ...
    'mean_total_old','mean_total_new','mean_total_emp', ...
    'relative_error_mean_old','relative_error_mean_new', ...
    'mean_old_phi','mean_new_phi','mean_emp_phi', ...
    'max_empirical_difference', ...
    '-v7.3');

fprintf('Resultats sauvegardes dans %s\n',output_file);

%% ============================================================
%  FONCTIONS LOCALES
%% ============================================================

function path_out = find_first_file(candidates,label)

    path_out = '';

    for k = 1:numel(candidates)
        if isfile(candidates{k})
            path_out = candidates{k};
            return;
        end
    end

    error('Fichier %s introuvable.',label);
end

function require_fields(S,names,file_name)

    for k = 1:numel(names)
        if ~isfield(S,names{k})
            error('Le fichier %s doit contenir %s.', ...
                file_name,names{k});
        end
    end
end

function assert_close(a,b,tol,name)

    scale = max([abs(a),abs(b),1]);

    if abs(a-b) > tol*scale
        error('%s differe : %.12g contre %.12g.', ...
            name,a,b);
    end
end

function x = row_vector(x)
    x = double(x(:).');
end

function x = col_vector(x)
    x = double(x(:));
end

function check_size(X,Nt,Nb,name)

    if ~isequal(size(X),[Nt,Nb])
        error('%s doit etre de taille %d x %d.', ...
            name,Nt,Nb);
    end
end
