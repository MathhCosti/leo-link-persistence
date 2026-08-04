%% N2_phi.m
% Nombre moyen de dimeres isoles en fonction de la latitude
% Walker Delta a uniformite orbitale.
%
% Comparaison :
%   - theorie locale estimee par Monte-Carlo conditionnel ;
%   - mesure empirique sur des graphes complets.
%
% Convention locale :
% chaque dimere contribue pour 1/2 dans la tranche de latitude de
% chacune de ses deux extremites. L'integrale sur la latitude redonne
% ainsi le nombre total de dimeres.
%
% Entree :
%   plink_phi_results.mat produit par plink_phi_mc_corrige.m
%
% Sortie :
%   N2_phi_results.mat

clear; clc; close all;

%% ============================================================
%  1. Chargement et parametres
%% ============================================================

script_dir = fileparts(mfilename('fullpath'));
input_file = fullfile(script_dir,'..', 'Valeurs locales', 'plink_phi_results.mat');

if ~isfile(input_file)
    input_file = fullfile(script_dir,'..','plink_phi_results.mat');
end

if ~isfile(input_file)
    error('Fichier plink_phi_results.mat introuvable.');
end

S = load(input_file);

required_fields = { ...
    'R','inc','dmax','N','phi_vals','phi_edges', ...
    'dphi_bins','phi_bin_probability'};

for k = 1:numel(required_fields)
    if ~isfield(S,required_fields{k})
        error(['La variable %s est absente. Relancez d''abord ', ...
               'plink_phi_mc_corrige.m.'],required_fields{k});
    end
end

R = double(S.R);
inc = double(S.inc);
dmax = double(S.dmax);
N = double(S.N);

phi_vals = S.phi_vals(:).';
phi_edges = S.phi_edges(:).';
dphi_bins = S.dphi_bins(:).';
phi_bin_probability = S.phi_bin_probability(:).';

n_phi = numel(phi_vals);
cmax = 1-dmax^2/(2*R^2);

% Monte-Carlo theorique conditionnel
n_outer_theory = 5000;
n_probe_union = 3000;

% Monte-Carlo empirique sur graphes complets
if isfield(S,'n_sim')
    n_sim_emp = double(S.n_sim);
else
    n_sim_emp = 500;
end

if isfield(S,'N_mean_target')
    N_mean_target = double(S.N_mean_target);
else
    N_mean_target = N;
end

rng_seed = 21;
rng(rng_seed);

%% ============================================================
%  2. Theorie locale des dimeres
%
% Pour X1 fixe a la latitude phi :
%
% h2(phi) = E[
%   1_{X1~X2} (1-q_union,2(X1,X2))^(N-2)
%   | Phi1=phi
% ]
%
% Puis :
% N2_b^th = C(N,2) P(Phi dans b) h2(phi_b).
%% ============================================================

h2_phi = zeros(1,n_phi);
h2_phi_sem = zeros(1,n_phi);

for b = 1:n_phi

    phi1 = phi_vals(b);
    contributions = zeros(n_outer_theory,1);

    for m = 1:n_outer_theory

        X1 = sample_satellite_given_latitude(phi1,inc);

        Omega2 = 2*pi*rand;
        u2 = 2*pi*rand;
        X2 = orbital_position(Omega2,u2,inc);

        linked12 = dot(X1,X2) >= cmax;

        if ~linked12
            continue;
        end

        probes = sample_orbital_points(n_probe_union,inc);

        linked_to_1 = probes*X1.' >= cmax;
        linked_to_2 = probes*X2.' >= cmax;

        q_union_2 = mean(linked_to_1 | linked_to_2);

        contributions(m) = ...
            (1-q_union_2)^max(N-2,0);
    end

    h2_phi(b) = mean(contributions);
    h2_phi_sem(b) = std(contributions)/sqrt(n_outer_theory);
end

comb_N_2 = N*(N-1)/2;

N2_per_bin_th = ...
    comb_N_2 .* phi_bin_probability .* h2_phi;

N2_density_phi_th = ...
    N2_per_bin_th ./ dphi_bins;

N2_total_th = sum(N2_per_bin_th);

%% ============================================================
%  3. Mesure empirique sur des graphes complets
%% ============================================================

N_draws_emp = zeros(n_sim_emp,1);
N2_total_sim = zeros(n_sim_emp,1);
N2_per_bin_sim = zeros(n_sim_emp,n_phi);

for s = 1:n_sim_emp

    Ns = max(draw_poisson(N_mean_target),2);
    N_draws_emp(s) = Ns;

    Omega = 2*pi*rand(Ns,1);
    u = 2*pi*rand(Ns,1);
    points = orbital_position(Omega,u,inc);

    gram = points*points.';
    adjacency = gram >= cmax;
    adjacency(1:Ns+1:end) = false;
    adjacency = adjacency | adjacency.';

    G = graph(adjacency);
    labels = conncomp(G);
    component_sizes = accumarray(labels(:),1);

    latitudes = asin(max(min(points(:,3),1),-1));
    bin_id = discretize(latitudes,phi_edges);

    dimer_labels = find(component_sizes == 2);
    N2_total_sim(s) = numel(dimer_labels);

    for q = 1:numel(dimer_labels)
        members = find(labels == dimer_labels(q));

        % Attribution symetrique : 1/2 par extremite.
        for r = 1:2
            bin = bin_id(members(r));
            if ~isnan(bin)
                N2_per_bin_sim(s,bin) = ...
                    N2_per_bin_sim(s,bin)+1/2;
            end
        end
    end
end

N_emp = mean(N_draws_emp);

N2_per_bin_emp = mean(N2_per_bin_sim,1);
N2_density_phi_emp = N2_per_bin_emp./dphi_bins;
N2_total_emp = mean(N2_total_sim);

N2_per_bin_emp_std = std(N2_per_bin_sim,0,1);
N2_per_bin_emp_sem = N2_per_bin_emp_std/sqrt(n_sim_emp);

%% ============================================================
%  4. Diagnostic
%% ============================================================

valid = isfinite(N2_density_phi_th) ...
    & isfinite(N2_density_phi_emp) ...
    & N2_density_phi_th > 0;

ratio_emp_th_phi = nan(1,n_phi);
ratio_emp_th_phi(valid) = ...
    N2_density_phi_emp(valid)./N2_density_phi_th(valid);

rmse_density = sqrt(mean( ...
    (N2_density_phi_emp(valid)-N2_density_phi_th(valid)).^2));

relative_error_total = ...
    abs(N2_total_emp-N2_total_th)/max(abs(N2_total_th),eps);

%% ============================================================
%  5. Figures
%% ============================================================

figure;
plot(rad2deg(phi_vals),N2_density_phi_th,'LineWidth',2);
hold on;
plot(rad2deg(phi_vals),N2_density_phi_emp,'--','LineWidth',1.8);
grid on;
xlabel('Latitude \phi (deg)');
ylabel('Densite moyenne de dimeres par radian');
title('Dimeres isoles selon la latitude');
legend('Theorie','Monte-Carlo','Location','best');
hold off;

figure;
plot(rad2deg(phi_vals),N2_per_bin_th,'LineWidth',2);
hold on;
plot(rad2deg(phi_vals),N2_per_bin_emp,'--','LineWidth',1.8);
grid on;
xlabel('Latitude \phi (deg)');
ylabel('Nombre moyen de dimeres dans la tranche');
title('Nombre moyen de dimeres par tranche');
legend('Theorie','Monte-Carlo','Location','best');
hold off;

figure;
plot(rad2deg(phi_vals),ratio_emp_th_phi,'LineWidth',1.8);
hold on;
yline(1,'--','Accord parfait');
grid on;
xlabel('Latitude \phi (deg)');
ylabel('Empirique / theorie');
title('Qualite du modele local des dimeres');
hold off;

%% ============================================================
%  6. Affichage
%% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' N2(phi) - DIMERES THEORIE / MONTE-CARLO\n');
fprintf('============================================================\n');
fprintf('N moyen venant de p_link            : %.8f\n',N);
fprintf('N moyen tire dans ce script         : %.8f\n',N_emp);
fprintf('Simulations empiriques              : %d\n',n_sim_emp);
fprintf('Tirages conditionnels par latitude  : %d\n',n_outer_theory);
fprintf('Sondes pour q_union,2               : %d\n',n_probe_union);
fprintf('N2 total theorique                  : %.10f\n',N2_total_th);
fprintf('N2 total empirique                  : %.10f\n',N2_total_emp);
fprintf('Erreur relative totale             : %.3e\n',relative_error_total);
fprintf('RMSE densite locale                : %.3e\n',rmse_density);
fprintf('============================================================\n');

%% ============================================================
%  7. Sauvegarde
%% ============================================================

output_file = fullfile(script_dir,'N2_phi_results.mat');

save(output_file, ...
    'R','inc','dmax','N','N_emp','N_mean_target', ...
    'phi_vals','phi_edges','dphi_bins','phi_bin_probability', ...
    'n_outer_theory','n_probe_union','n_sim_emp','rng_seed', ...
    'h2_phi','h2_phi_sem','comb_N_2', ...
    'N2_density_phi_th','N2_density_phi_emp', ...
    'N2_per_bin_th','N2_per_bin_emp','N2_per_bin_sim', ...
    'N2_per_bin_emp_std','N2_per_bin_emp_sem', ...
    'N2_total_th','N2_total_emp','N2_total_sim', ...
    'ratio_emp_th_phi','rmse_density','relative_error_total');

fprintf('Resultats sauvegardes dans %s\n',output_file);

%% ============================================================
%  Fonctions locales
%% ============================================================

function X = sample_satellite_given_latitude(phi,inc)
    s = sin(phi)/sin(inc);
    s = max(min(s,1),-1);

    u_a = asin(s);
    u_b = pi-u_a;

    if rand < 0.5
        u = u_a;
    else
        u = u_b;
    end

    Omega = 2*pi*rand;
    X = orbital_position(Omega,u,inc);
end

function X = sample_orbital_points(n,inc)
    Omega = 2*pi*rand(n,1);
    u = 2*pi*rand(n,1);
    X = orbital_position(Omega,u,inc);
end

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
