%% N1_phi.m
% Nombre moyen de satellites isoles en fonction de la latitude
% Walker Delta a uniformite orbitale.
%
% La valeur N, la geometrie orbitale, la grille de latitude et le nombre
% de simulations sont recuperes depuis plink_phi_results.mat afin de
% rester coherents avec plink_phi.m, N2_phi.m et N3_phi.m.
%
% Theorie locale :
%   N1_b^th = N * P(Phi dans b) * (1-p_link(phi_b))^(N-1)
%
% Entree : plink_phi_results.mat
% Sortie : N1_phi_results.mat

clear; clc; close all;

%% 1. Chargement des parametres de plink_phi
script_dir = fileparts(mfilename('fullpath'));

candidate_files = {
    fullfile(script_dir,'..','Valeurs locales','plink_phi_results.mat')
    fullfile(script_dir,'plink_phi_results.mat')
    fullfile(script_dir,'..','plink_phi_results.mat')
};

input_file = '';
for k = 1:numel(candidate_files)
    if isfile(candidate_files{k})
        input_file = candidate_files{k};
        break;
    end
end

if isempty(input_file)
    error('Fichier plink_phi_results.mat introuvable.');
end

S = load(input_file);
required_fields = {'R','inc','dmax','N','phi_vals','phi_edges', ...
    'dphi_bins','phi_bin_probability','p_link_phi'};

for k = 1:numel(required_fields)
    if ~isfield(S,required_fields{k})
        error(['La variable %s est absente. Relancez d''abord ', ...
            'plink_phi_mc_corrige.m.'],required_fields{k});
    end
end

R = double(S.R);
inc = double(S.inc);
inc_deg = rad2deg(inc);
dmax = double(S.dmax);
N = double(S.N); % moyenne des tirages issue de plink_phi

phi_vals = S.phi_vals(:).';
phi_edges = S.phi_edges(:).';
dphi_bins = S.dphi_bins(:).';
phi_bin_probability = S.phi_bin_probability(:).';
p_link_phi = S.p_link_phi(:).';

n_phi = numel(phi_vals);
cmax = 1-dmax^2/(2*R^2);

if isfield(S,'n_sim')
    n_realizations = double(S.n_sim);
else
    n_realizations = 500;
end

if isfield(S,'N_mean_target')
    N_mean_target = double(S.N_mean_target);
else
    N_mean_target = N;
end

rng_seed = 2;
rng(rng_seed);

%% 2. Prediction theorique locale
p_isolated_phi_th = (1-p_link_phi).^(N-1);

N1_per_bin_th = ...
    N .* phi_bin_probability .* p_isolated_phi_th;

N1_density_phi_th = N1_per_bin_th ./ dphi_bins;
N1_total_th_local = sum(N1_per_bin_th);

if isfield(S,'p_link_global_from_phi')
    p_link_global = double(S.p_link_global_from_phi);
else
    p_link_global = sum(phi_bin_probability.*p_link_phi);
end

N1_total_th_global = N*(1-p_link_global)^(N-1);

%% 3. Monte-Carlo empirique
N_draws_emp = zeros(n_realizations,1);
N1_total_sim = zeros(n_realizations,1);
N1_per_bin_sim = zeros(n_realizations,n_phi);

nodes_per_bin = zeros(1,n_phi);
isolated_sum_bin = zeros(1,n_phi);

for r = 1:n_realizations
    Ns = max(draw_poisson(N_mean_target),1);
    N_draws_emp(r) = Ns;

    Omega = 2*pi*rand(Ns,1);
    u = 2*pi*rand(Ns,1);
    points = orbital_position(Omega,u,inc);

    gram = points*points.';
    adjacency = gram >= cmax;
    adjacency(1:Ns+1:end) = false;
    adjacency = adjacency | adjacency.';

    degree = sum(adjacency,2);
    is_isolated = degree == 0;
    N1_total_sim(r) = nnz(is_isolated);

    latitudes = asin(max(min(points(:,3),1),-1));
    bin_id = discretize(latitudes,phi_edges);

    for b = 1:n_phi
        mask = bin_id == b;
        nodes_per_bin(b) = nodes_per_bin(b)+nnz(mask);
        isolated_sum_bin(b) = isolated_sum_bin(b)+nnz(mask & is_isolated);
        N1_per_bin_sim(r,b) = nnz(mask & is_isolated);
    end
end

N_emp = mean(N_draws_emp);
N1_per_bin_emp = mean(N1_per_bin_sim,1);
N1_density_phi_emp = N1_per_bin_emp ./ dphi_bins;
N1_total_emp = mean(N1_total_sim);

p_isolated_phi_emp = isolated_sum_bin ./ max(nodes_per_bin,1);
N1_per_bin_emp_std = std(N1_per_bin_sim,0,1);
N1_per_bin_emp_sem = N1_per_bin_emp_std/sqrt(n_realizations);

%% 4. Diagnostic
valid = isfinite(N1_density_phi_th) & ...
    isfinite(N1_density_phi_emp) & N1_density_phi_th > 0;

ratio_emp_th_phi = nan(1,n_phi);
ratio_emp_th_phi(valid) = ...
    N1_density_phi_emp(valid)./N1_density_phi_th(valid);

rmse_density = sqrt(mean((N1_density_phi_emp(valid) ...
    -N1_density_phi_th(valid)).^2));

relative_error_total = abs(N1_total_emp-N1_total_th_local) ...
    / max(abs(N1_total_th_local),eps);

%% 5. Figures
figure;
plot(rad2deg(phi_vals),p_isolated_phi_th,'LineWidth',2);
hold on;
plot(rad2deg(phi_vals),p_isolated_phi_emp,'--','LineWidth',1.8);
grid on;
xlabel('Latitude \phi (deg)');
ylabel('Probabilite d''etre isole');
title('Probabilite locale d''isolement');
legend('Theorie','Monte-Carlo','Location','best');
hold off;

figure;
plot(rad2deg(phi_vals),N1_density_phi_th,'LineWidth',2);
hold on;
plot(rad2deg(phi_vals),N1_density_phi_emp,'--','LineWidth',1.8);
grid on;
xlabel('Latitude \phi (deg)');
ylabel('Densite moyenne d''isoles par radian');
title('Satellites isoles selon la latitude');
legend('Theorie','Monte-Carlo','Location','best');
hold off;

figure;
plot(rad2deg(phi_vals),N1_per_bin_th,'LineWidth',2);
hold on;
plot(rad2deg(phi_vals),N1_per_bin_emp,'--','LineWidth',1.8);
grid on;
xlabel('Latitude \phi (deg)');
ylabel('Nombre moyen d''isoles dans la tranche');
title('Nombre moyen de satellites isoles par tranche');
legend('Theorie','Monte-Carlo','Location','best');
hold off;

figure;
plot(rad2deg(phi_vals),ratio_emp_th_phi,'LineWidth',1.8);
hold on;
yline(1,'--','Accord parfait');
grid on;
xlabel('Latitude \phi (deg)');
ylabel('Empirique / theorie');
title('Qualite du modele local de N_1');
hold off;

figure;
histogram(N1_total_sim,30,'Normalization','pdf');
hold on;
xline(N1_total_th_local,'--','Prediction locale','LineWidth',1.8);
xline(N1_total_th_global,':','Prediction globale','LineWidth',1.8);
grid on;
xlabel('Nombre de satellites isoles');
ylabel('Densite');
title('Distribution du nombre total de satellites isoles');
hold off;

%% 6. Affichage
fprintf('\n');
fprintf('============================================================\n');
fprintf(' N1(phi) - ISOLES THEORIE / MONTE-CARLO\n');
fprintf('============================================================\n');
fprintf('Fichier source                       : %s\n',input_file);
fprintf('N moyen venant de p_link             : %.8f\n',N);
fprintf('N moyen tire dans ce script          : %.8f\n',N_emp);
fprintf('Moyenne cible de N                   : %.8f\n',N_mean_target);
fprintf('Nombre de simulations                : %d\n',n_realizations);
fprintf('------------------------------------------------------------\n');
fprintf('N1 total theorique local             : %.10f\n',N1_total_th_local);
fprintf('N1 total theorique global            : %.10f\n',N1_total_th_global);
fprintf('N1 total empirique                   : %.10f\n',N1_total_emp);
fprintf('Erreur relative theorie locale       : %.3e\n',relative_error_total);
fprintf('RMSE densite locale                  : %.3e\n',rmse_density);
fprintf('============================================================\n');

%% 7. Sauvegarde
output_file = fullfile(script_dir,'N1_phi_results.mat');

save(output_file, ...
    'input_file','R','inc','inc_deg','dmax','N','N_emp', ...
    'N_mean_target','n_realizations','rng_seed', ...
    'phi_vals','phi_edges','dphi_bins','phi_bin_probability', ...
    'p_link_phi','p_link_global', ...
    'p_isolated_phi_th','p_isolated_phi_emp', ...
    'N1_density_phi_th','N1_density_phi_emp', ...
    'N1_per_bin_th','N1_per_bin_emp','N1_per_bin_sim', ...
    'N1_per_bin_emp_std','N1_per_bin_emp_sem', ...
    'N1_total_th_local','N1_total_th_global', ...
    'N1_total_emp','N1_total_sim', ...
    'N_draws_emp','nodes_per_bin','isolated_sum_bin', ...
    'ratio_emp_th_phi','rmse_density','relative_error_total');

fprintf('Resultats sauvegardes dans %s\n',output_file);

%% Fonctions locales
function X = orbital_position(Omega,u,inc)
    Omega = Omega(:);
    u = u(:);

    cO = cos(Omega);
    sO = sin(Omega);
    cu = cos(u);
    su = sin(u);

    X = [ ...
        cO.*cu-sO.*su*cos(inc), ...
        sO.*cu+cO.*su*cos(inc), ...
        su*sin(inc)];
end

function n = draw_poisson(lambda)
    L = exp(-lambda);
    p = 1;
    k = 0;

    while p > L
        k = k+1;
        p = p*rand;
    end

    n = k-1;
end
