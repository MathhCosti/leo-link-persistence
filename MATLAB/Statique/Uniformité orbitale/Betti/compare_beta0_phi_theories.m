%% compare_beta0_phi_theories.m
% Compare sur un meme graphe :
%   1) ancienne theorie beta0(phi) jusqu'aux trimeres + composantes macro ;
%   2) nouvelle theorie PPP par expansion en tailles de composantes ;
%   3) beta0(phi) empirique direct.
%
% Entrees :
%   betti_phi_results.mat
%   beta0_phi_results.mat
%
% Le code utilise comme grille de comparaison celle de beta0_phi_results.mat
% et interpole l'ancienne theorie si les grilles different.

clear; clc; close all;

%% ============================================================
%  1. CHARGEMENT
%% ============================================================

script_dir = fileparts(mfilename('fullpath'));

old_candidates = {
    fullfile(script_dir,'betti_phi_results.mat')
    fullfile(script_dir,'betti_phi_results(1).mat')
    fullfile(script_dir,'..','betti_phi_results.mat')
};

new_candidates = {
    fullfile(script_dir,'beta0_phi_results.mat')
    fullfile(script_dir,'beta0_phi_th_results.mat')
    fullfile(script_dir,'..','beta0_phi_results.mat')
    fullfile(script_dir,'..','beta0_phi_th_results.mat')
};

old_file = find_first_file(old_candidates,'betti_phi_results.mat');
new_file = find_first_file(new_candidates,'beta0_phi_results.mat');

Sold = load(old_file);
Snew = load(new_file);

%% ============================================================
%  2. VERIFICATION DES PARAMETRES
%% ============================================================

R_old = double(Sold.R);
R_new = double(Snew.R);

inc_old = double(Sold.inc);
inc_new = double(Snew.inc);

dmax_old = double(Sold.dmax);
dmax_new = double(Snew.dmax);

N_old = double(Sold.N);
N_new = double(Snew.N);

assert_close(R_old,R_new,1e-10,'R');
assert_close(inc_old,inc_new,1e-10,'inc');
assert_close(dmax_old,dmax_new,1e-10,'dmax');

if abs(N_old-N_new) > 0.02*max([N_old,N_new,1])
    warning('N differe entre les deux fichiers : %.3f contre %.3f.', ...
        N_old,N_new);
end

%% ============================================================
%  3. GRILLE DE REFERENCE : NOUVELLE THEORIE
%% ============================================================

phi = rowvec(Snew.phi);

if isfield(Snew,'dphi_bins')
    dphi = rowvec(Snew.dphi_bins);
else
    error('beta0_phi_results.mat doit contenir dphi_bins.');
end

phi_deg = rad2deg(phi);

%% ============================================================
%  4. NOUVELLE THEORIE PPP
%% ============================================================

if isfield(Snew,'beta0_density_per_rad_th')
    beta_new = rowvec(Snew.beta0_density_per_rad_th);
else
    error(['beta0_phi_results.mat doit contenir ', ...
           'beta0_density_per_rad_th.']);
end

%% ============================================================
%  5. EMPIRIQUE
%
% On prend de preference l'empirique du nouveau fichier car il utilise
% exactement la meme grille que la nouvelle theorie.
%% ============================================================

if isfield(Snew,'beta0_density_per_rad_emp')
    beta_emp = rowvec(Snew.beta0_density_per_rad_emp);
elseif isfield(Sold,'betti0_density_phi_emp')
    phi_emp_old = rowvec(Sold.phi_vals);
    beta_emp_old = rowvec(Sold.betti0_density_phi_emp);
    beta_emp = interp1(phi_emp_old,beta_emp_old,phi,'linear','extrap');
else
    error('Aucune courbe empirique beta0(phi) reconnue.');
end

%% ============================================================
%  6. ANCIENNE THEORIE : TRIMERES + MACRO
%% ============================================================

if ~isfield(Sold,'phi_vals') || ~isfield(Sold,'betti0_density_phi_th')
    error(['betti_phi_results.mat doit contenir phi_vals et ', ...
           'betti0_density_phi_th.']);
end

phi_old = rowvec(Sold.phi_vals);
beta_old_raw = rowvec(Sold.betti0_density_phi_th);

% Re-echantillonnage sur la grille de la nouvelle theorie.
beta_old = interp1(phi_old,beta_old_raw,phi,'linear','extrap');

%% ============================================================
%  7. DIAGNOSTICS
%% ============================================================

valid = isfinite(beta_emp) & isfinite(beta_old) & isfinite(beta_new);

rmse_old = sqrt(mean((beta_old(valid)-beta_emp(valid)).^2));
rmse_new = sqrt(mean((beta_new(valid)-beta_emp(valid)).^2));

mae_old = mean(abs(beta_old(valid)-beta_emp(valid)));
mae_new = mean(abs(beta_new(valid)-beta_emp(valid)));

% Totaux
if isfield(Sold,'betti0_total_th')
    beta0_total_old = double(Sold.betti0_total_th);
else
    beta0_total_old = trapz(phi,beta_old);
end

if isfield(Snew,'beta0_integrated_phi')
    beta0_total_new = double(Snew.beta0_integrated_phi);
else
    beta0_total_new = trapz(phi,beta_new);
end

if isfield(Snew,'beta0_total_emp')
    beta0_total_emp = double(Snew.beta0_total_emp);
elseif isfield(Sold,'betti0_total_emp')
    beta0_total_emp = double(Sold.betti0_total_emp);
else
    beta0_total_emp = trapz(phi,beta_emp);
end

err_total_old = ...
    abs(beta0_total_old-beta0_total_emp)/max(abs(beta0_total_emp),eps);

err_total_new = ...
    abs(beta0_total_new-beta0_total_emp)/max(abs(beta0_total_emp),eps);

fprintf('\n============================================================\n');
fprintf(' COMPARAISON beta0(phi) : DEUX THEORIES / EMPIRIQUE\n');
fprintf('============================================================\n');
fprintf('N ancien / nouveau                  : %.3f / %.3f\n', ...
    N_old,N_new);
fprintf('Inclinaison                         : %.3f deg\n',rad2deg(inc_new));
fprintf('dmax                                : %.3f km\n',dmax_new);
fprintf('------------------------------------------------------------\n');
fprintf('beta0 total empirique               : %.8f\n',beta0_total_emp);
fprintf('beta0 total ancienne theorie        : %.8f\n',beta0_total_old);
fprintf('beta0 total nouvelle theorie PPP    : %.8f\n',beta0_total_new);
fprintf('Erreur totale ancienne theorie      : %.2f %%\n',100*err_total_old);
fprintf('Erreur totale nouvelle theorie      : %.2f %%\n',100*err_total_new);
fprintf('------------------------------------------------------------\n');
fprintf('RMSE locale ancienne theorie        : %.8f\n',rmse_old);
fprintf('RMSE locale nouvelle theorie        : %.8f\n',rmse_new);
fprintf('MAE locale ancienne theorie         : %.8f\n',mae_old);
fprintf('MAE locale nouvelle theorie         : %.8f\n',mae_new);
fprintf('============================================================\n');

%% ============================================================
%  8. GRAPHE PRINCIPAL
%% ============================================================

figure;
plot(phi_deg,beta_emp,'k--','LineWidth',1.9, ...
    'DisplayName','Empirique');
hold on;

plot(phi_deg,beta_old,'LineWidth',1.8, ...
    'DisplayName','Ancienne theorie : macro + N_1+N_2+N_3');

plot(phi_deg,beta_new,'LineWidth',2.0, ...
    'DisplayName','Nouvelle theorie PPP : \Sigma_s C_s');

grid on;
xlabel('Latitude \phi (deg)');
ylabel('Densite moyenne de composantes par radian');
title('\beta_0(\phi) : comparaison des deux theories avec l''empirique');
legend('Location','best');
hold off;

%% ============================================================
%  9. GRAPHE DES ERREURS LOCALES
%% ============================================================

figure;
plot(phi_deg,beta_old-beta_emp,'LineWidth',1.6, ...
    'DisplayName','Ancienne theorie - empirique');
hold on;
plot(phi_deg,beta_new-beta_emp,'LineWidth',1.8, ...
    'DisplayName','Nouvelle theorie - empirique');
yline(0,'k--','Accord parfait');
grid on;
xlabel('Latitude \phi (deg)');
ylabel('Erreur locale (composantes/rad)');
title('Erreur locale des deux approximations de \beta_0');
legend('Location','best');
hold off;

%% ============================================================
%  10. SAUVEGARDE
%% ============================================================

save('compare_beta0_phi_theories_results.mat', ...
    'phi','phi_deg','beta_emp','beta_old','beta_new', ...
    'beta0_total_emp','beta0_total_old','beta0_total_new', ...
    'err_total_old','err_total_new', ...
    'rmse_old','rmse_new','mae_old','mae_new', ...
    'old_file','new_file');

fprintf('Resultats sauvegardes dans compare_beta0_phi_theories_results.mat\n');

%% ============================================================
%  FONCTIONS LOCALES
%% ============================================================

function f = find_first_file(candidates,label)
    f = '';
    for k = 1:numel(candidates)
        if isfile(candidates{k})
            f = candidates{k};
            return;
        end
    end
    error('Fichier %s introuvable.',label);
end

function x = rowvec(x)
    x = double(x(:).');
end

function assert_close(a,b,tol,name)
    scale = max([abs(a),abs(b),1]);
    if abs(a-b) > tol*scale
        error('%s differe entre les deux fichiers : %.12g contre %.12g.', ...
            name,a,b);
    end
end
