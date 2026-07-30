%% liens_bridge_lambda.m
% Comparaison temporelle du facteur :
%
%       chi_merge(t) = (N - beta0(t)) / E(t)
%
% pour un Walker-Delta a uniformite spatiale initiale.
%
% Courbes affichees :
%   1. Empirique :
%        N, beta0(t) et E(t) issus de analysis_temp_results.mat.
%
%   2. Theorie avec beta0 constant "isoles" :
%        beta0 = 1 + E[N1(0)].
%
%   3. Theorie avec beta0 constant "isoles + dimeres" :
%        beta0 = 1 + E[N1(0)] + E[N2(0)].
%
%   4. Theorie avec beta0 constant "isoles + dimeres + trimeres" :
%        beta0 = 1 + E[N1(0)] + E[N2(0)] + E[N3(0)].
%
% Pour les trois courbes theoriques, E(t) est pris dans
% liens_quadrature_results.mat.
%
% Sortie :
%   - figure unique de comparaison
%   - liens_bridge_temp_results.mat

clear; clc; close all;

%% ============================================================
%  LOCALISATION DES FICHIERS
%% ============================================================

script_dir = fileparts(mfilename('fullpath'));
links_file = fullfile(script_dir, '..', 'Paramètres', 'Nombre liens', 'liens_quadrature_results.mat');
analysis_file = fullfile(script_dir, '..', 'analysis_temp_results.mat');

if isempty(analysis_file)
    error('Fichier analysis_temp_results.mat introuvable.');
end

if isempty(links_file)
    error('Fichier liens_quadrature_results.mat introuvable.');
end

fprintf('Fichier empirique  : %s\n',analysis_file);
fprintf('Fichier theorique  : %s\n',links_file);

Semp = load(analysis_file);
Sth = load(links_file);

%% ============================================================
%  DONNEES EMPIRIQUES
%% ============================================================

required_emp = {'time_values','N','beta0','num_edges'};

for k = 1:numel(required_emp)
    if ~isfield(Semp,required_emp{k})
        error('Variable %s absente de %s.', ...
            required_emp{k},analysis_file);
    end
end

t_emp = double(Semp.time_values(:));
N_emp = double(Semp.N);
beta0_emp = double(Semp.beta0(:));
E_emp = double(Semp.num_edges(:));

if numel(t_emp) ~= numel(beta0_emp) || numel(t_emp) ~= numel(E_emp)
    error('Tailles incompatibles dans les donnees empiriques.');
end

chi_emp = NaN(size(E_emp));
valid_emp = E_emp > 0;

chi_emp(valid_emp) = ...
    (N_emp-beta0_emp(valid_emp)) ./ E_emp(valid_emp);

chi_emp = min(max(chi_emp,0),1);

%% ============================================================
%  DONNEES THEORIQUES
%% ============================================================

required_th = {'time_values','L_theory_time','p_link_theory_time'};

for k = 1:numel(required_th)
    if ~isfield(Sth,required_th{k})
        error('Variable %s absente de %s.', ...
            required_th{k},links_file);
    end
end

t_th = double(Sth.time_values(:));
E_th = double(Sth.L_theory_time(:));
p_link_th = double(Sth.p_link_theory_time(:));

if isfield(Sth,'N_mean_theory')
    N_th = double(Sth.N_mean_theory);
elseif isfield(Sth,'N_all')
    N_th = mean(double(Sth.N_all(:)),'omitnan');
else
    error('Impossible de determiner N theorique.');
end

if numel(t_th) ~= numel(E_th) || numel(t_th) ~= numel(p_link_th)
    error('Tailles incompatibles dans les donnees theoriques.');
end

%% ============================================================
%  beta0 THEORIQUES CONSTANTS, EVALUES A t = 0
%% ============================================================

p_link_0 = p_link_th(1);

% Coefficients geometriques.
if isfield(Semp,'c2')
    c2_union = double(Semp.c2);
else
    c2_union = 1 + 3*sqrt(3)/(4*pi);
end

if isfield(Semp,'c3_conn')
    c3_conn = double(Semp.c3_conn);
else
    c3_conn = 1 + 3*sqrt(3)/(2*pi);
end

if isfield(Semp,'c3_union')
    c3_union = double(Semp.c3_union);
else
    c3_union = 1.80;
end

% Isoles.
N1_0 = N_th * max(1-p_link_0,0)^(N_th-1);

% Dimeres.
if N_th >= 2
    C2 = exp(gammaln(N_th+1)-gammaln(3)-gammaln(N_th-1));

    N2_0 = C2 * p_link_0 ...
        * max(1-c2_union*p_link_0,0)^(N_th-2);
else
    N2_0 = 0;
end

% Trimeres.
if N_th >= 3
    C3 = exp(gammaln(N_th+1)-gammaln(4)-gammaln(N_th-2));

    p_conn_3_0 = min(max(c3_conn*p_link_0^2,0),1);

    N3_0 = C3 * p_conn_3_0 ...
        * max(1-c3_union*p_link_0,0)^(N_th-3);
else
    N3_0 = 0;
end

% Une composante macroscopique est ajoutee dans chaque approximation.
beta0_iso = min(max(1+N1_0,1),N_th);
beta0_iso_dim = min(max(1+N1_0+N2_0,1),N_th);
beta0_iso_dim_tri = min(max(1+N1_0+N2_0+N3_0,1),N_th);

%% ============================================================
%  COURBES THEORIQUES DE (N-beta0)/E(t)
%% ============================================================

chi_th_iso = NaN(size(E_th));
chi_th_iso_dim = NaN(size(E_th));
chi_th_iso_dim_tri = NaN(size(E_th));

valid_th = E_th > 0;

chi_th_iso(valid_th) = ...
    (N_th-beta0_iso) ./ E_th(valid_th);

chi_th_iso_dim(valid_th) = ...
    (N_th-beta0_iso_dim) ./ E_th(valid_th);

chi_th_iso_dim_tri(valid_th) = ...
    (N_th-beta0_iso_dim_tri) ./ E_th(valid_th);

chi_th_iso = min(max(chi_th_iso,0),1);
chi_th_iso_dim = min(max(chi_th_iso_dim,0),1);
chi_th_iso_dim_tri = min(max(chi_th_iso_dim_tri,0),1);

%% ============================================================
%  ALIGNEMENT DES COURBES SUR UNE MEME GRILLE
%% ============================================================

chi_th_iso_on_emp = interp1( ...
    t_th,chi_th_iso,t_emp,'linear',NaN);

chi_th_iso_dim_on_emp = interp1( ...
    t_th,chi_th_iso_dim,t_emp,'linear',NaN);

chi_th_iso_dim_tri_on_emp = interp1( ...
    t_th,chi_th_iso_dim_tri,t_emp,'linear',NaN);

%% ============================================================
%  GRAPHE UNIQUE
%% ============================================================

figure;
hold on;
grid on;
box on;

plot(t_emp,chi_emp,'LineWidth',1.8, ...
    'DisplayName','Empirique');

plot(t_emp,chi_th_iso_on_emp,'--','LineWidth',2.0, ...
    'DisplayName','Théorie : isolés');

plot(t_emp,chi_th_iso_dim_on_emp,'-.','LineWidth',2.0, ...
    'DisplayName','Théorie : isolés + dimères');

plot(t_emp,chi_th_iso_dim_tri_on_emp,':','LineWidth',2.4, ...
    'DisplayName','Théorie : isolés + dimères + trimères');

xlabel('Temps (s)');
ylabel('(N-\beta_0)/E');
title('Walker-Delta spatial : comparaison du facteur (N-\beta_0)/E');
legend('Location','best');
ylim([0,1.05]);
hold off;

%% ============================================================
%  AFFICHAGE CONSOLE
%% ============================================================

fprintf('\n=== Facteur (N-beta0)/E ===\n');
fprintf('N empirique                           : %.6f\n',N_emp);
fprintf('N theorique                           : %.6f\n',N_th);
fprintf('p_link theorique initial              : %.10f\n',p_link_0);
fprintf('E[N1(0)]                              : %.6f\n',N1_0);
fprintf('E[N2(0)]                              : %.6f\n',N2_0);
fprintf('E[N3(0)]                              : %.6f\n',N3_0);
fprintf('beta0 constant, isoles                : %.6f\n',beta0_iso);
fprintf('beta0 constant, isoles+dimeres        : %.6f\n',beta0_iso_dim);
fprintf('beta0 constant, isoles+dimeres+trim.  : %.6f\n',beta0_iso_dim_tri);
fprintf('Moyenne empirique de (N-beta0)/E      : %.6f\n', ...
    mean(chi_emp,'omitnan'));
fprintf('Moyenne theorique, isoles             : %.6f\n', ...
    mean(chi_th_iso,'omitnan'));
fprintf('Moyenne theorique, isoles+dimeres     : %.6f\n', ...
    mean(chi_th_iso_dim,'omitnan'));
fprintf('Moyenne theorique, jusqu''aux trimeres : %.6f\n', ...
    mean(chi_th_iso_dim_tri,'omitnan'));

%% ============================================================
%  SAUVEGARDE
%% ============================================================

save('liens_bridge_temp_results.mat', ...
    't_emp','t_th', ...
    'N_emp','N_th', ...
    'beta0_emp','E_emp','chi_emp', ...
    'E_th','p_link_th','p_link_0', ...
    'N1_0','N2_0','N3_0', ...
    'beta0_iso','beta0_iso_dim','beta0_iso_dim_tri', ...
    'chi_th_iso','chi_th_iso_dim','chi_th_iso_dim_tri', ...
    'chi_th_iso_on_emp','chi_th_iso_dim_on_emp', ...
    'chi_th_iso_dim_tri_on_emp', ...
    'analysis_file','links_file');

fprintf('\nResultats sauvegardes dans liens_bridge_temp_results.mat\n');

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
