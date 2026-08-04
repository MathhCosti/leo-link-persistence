%% betti_phi.m
% Approximation locale de beta_0 jusqu'aux trimeres
% Walker Delta a uniformite orbitale.
%
% Modele theorique :
%
%   B_0^Delta(phi)
%     ~= B_macro^Delta(phi)
%        + N_1^Delta(phi)
%        + N_2^Delta(phi)
%        + N_3^Delta(phi)
%
% avec deux composantes macroscopiques au total. Leur contribution locale
% est repartie selon la loi de latitude :
%
%   B_macro,b^th = 2 P(Phi dans la tranche b).
%
% Le code compare :
%   - l'approximation theorique jusqu'aux trimeres ;
%   - la valeur empirique de beta_0 mesuree directement sur les graphes
%     complets, toutes tailles de composantes confondues.
%
% Entrees :
%   N1_phi_results.mat
%   N2_phi_results.mat
%   N3_phi_results.mat
%
% Sortie :
%   betti_phi_results.mat
%
% IMPORTANT :
% les trois fichiers doivent correspondre aux memes parametres
% physiques et au meme nombre moyen de satellites. Les grilles de
% latitude peuvent differer : N1 est alors re-echantillonne sur la
% grille commune de N2 et N3 par conservation de la densite.

clear; clc; close all;

%% ============================================================
%  1. Chargement
%% ============================================================

script_dir = fileparts(mfilename('fullpath'));

file_N1 = find_result_file(script_dir,'N1_phi_results.mat');
file_N2 = find_result_file(script_dir,'N2_phi_results.mat');
file_N3 = find_result_file(script_dir,'N3_phi_results.mat');

S1 = load(file_N1);
S2 = load(file_N2);
S3 = load(file_N3);

%% ============================================================
%  2. Verification des parametres communs
%% ============================================================

N1_model = get_scalar_field(S1,{'N'});
N2_model = get_scalar_field(S2,{'N'});
N3_model = get_scalar_field(S3,{'N'});

inc1 = get_scalar_field(S1,{'inc','inc_rad'});
inc2 = get_scalar_field(S2,{'inc','inc_rad'});
inc3 = get_scalar_field(S3,{'inc','inc_rad'});

dmax1 = get_scalar_field(S1,{'dmax'});
dmax2 = get_scalar_field(S2,{'dmax'});
dmax3 = get_scalar_field(S3,{'dmax'});

R1 = get_scalar_field(S1,{'R'});
R2 = get_scalar_field(S2,{'R'});
R3 = get_scalar_field(S3,{'R'});

% Tolerances : N est une moyenne Monte-Carlo et peut legerement varier.
tol_N_rel = 2e-2;
tol_geom_rel = 1e-8;

assert_close_relative(N1_model,N2_model,tol_N_rel, ...
    'N differe entre N1 et N2');
assert_close_relative(N2_model,N3_model,tol_N_rel, ...
    'N differe entre N2 et N3');

assert_close_relative(inc1,inc2,tol_geom_rel, ...
    'L''inclinaison differe entre N1 et N2');
assert_close_relative(inc2,inc3,tol_geom_rel, ...
    'L''inclinaison differe entre N2 et N3');

assert_close_relative(dmax1,dmax2,tol_geom_rel, ...
    'dmax differe entre N1 et N2');
assert_close_relative(dmax2,dmax3,tol_geom_rel, ...
    'dmax differe entre N2 et N3');

assert_close_relative(R1,R2,tol_geom_rel, ...
    'R differe entre N1 et N2');
assert_close_relative(R2,R3,tol_geom_rel, ...
    'R differe entre N2 et N3');

% Valeur commune retenue
N = mean([N1_model,N2_model,N3_model]);
inc = mean([inc1,inc2,inc3]);
dmax = mean([dmax1,dmax2,dmax3]);
R = mean([R1,R2,R3]);

%% ============================================================
%  3. Grille commune : celle de N2 et N3
%% ============================================================

required_N2 = { ...
    'phi_vals','phi_edges','dphi_bins','phi_bin_probability', ...
    'N2_per_bin_th','N2_per_bin_emp'};

required_N3 = { ...
    'phi_vals','phi_edges','dphi_bins','phi_bin_probability', ...
    'N3_per_bin_th','N3_per_bin_emp'};

require_fields(S2,required_N2,'N2_phi_results.mat');
require_fields(S3,required_N3,'N3_phi_results.mat');

phi_vals = row_vector(S2.phi_vals);
phi_edges = row_vector(S2.phi_edges);
dphi_bins = row_vector(S2.dphi_bins);
phi_bin_probability = row_vector(S2.phi_bin_probability);

if max(abs(phi_vals-row_vector(S3.phi_vals))) > 1e-10 || ...
   max(abs(phi_edges-row_vector(S3.phi_edges))) > 1e-10
    error(['N2_phi_results.mat et N3_phi_results.mat doivent utiliser ', ...
           'la meme grille de latitude.']);
end

Nb = numel(phi_vals);

N2_per_bin_th = row_vector(S2.N2_per_bin_th);
N2_per_bin_emp = row_vector(S2.N2_per_bin_emp);

N3_per_bin_th = row_vector(S3.N3_per_bin_th);
N3_per_bin_emp = row_vector(S3.N3_per_bin_emp);

%% ============================================================
%  4. Recuperation de N1 et projection sur la grille commune
%% ============================================================

[N1_edges_source, ...
 N1_per_bin_th_source, ...
 N1_per_bin_emp_source] = extract_N1_per_bin(S1);

% Repartition conservatrice des nombres par tranche.
N1_per_bin_th = rebin_counts( ...
    N1_edges_source,N1_per_bin_th_source,phi_edges);

N1_per_bin_emp = rebin_counts( ...
    N1_edges_source,N1_per_bin_emp_source,phi_edges);

%% ============================================================
%  5. Contribution theorique des composantes macroscopiques
%
% Deux composantes macroscopiques au total, reparties selon la
% probabilite de presence d'un satellite a chaque latitude.
%% ============================================================

n_macro_components = 2;

beta_macro_per_bin_th = ...
    n_macro_components .* phi_bin_probability;

%% ============================================================
%  6. beta_0 theorique local jusqu'aux trimeres
%% ============================================================

betti0_per_bin_th = ...
    beta_macro_per_bin_th ...
    + N1_per_bin_th ...
    + N2_per_bin_th ...
    + N3_per_bin_th;

betti0_density_phi_th = ...
    betti0_per_bin_th ./ dphi_bins;

betti0_total_th = sum(betti0_per_bin_th);

%% ============================================================
%  6.b beta_0 empirique direct sur les graphes complets
%
% Chaque composante connexe de taille s contribue pour 1/s dans
% la tranche de latitude de chacun de ses s satellites. Ainsi,
% la somme sur toutes les tranches vaut exactement 1 pour chaque
% composante, et donc beta_0 pour chaque realisation.
%% ============================================================

% Parametres Monte-Carlo
if isfield(S1,'N_mean_target')
    N_mean_target = double(S1.N_mean_target);
elseif isfield(S2,'N_mean_target')
    N_mean_target = double(S2.N_mean_target);
elseif isfield(S3,'N_mean_target')
    N_mean_target = double(S3.N_mean_target);
else
    N_mean_target = N;
end

n_sim_candidates = [];
if isfield(S1,'n_realizations')
    n_sim_candidates(end+1) = double(S1.n_realizations);
end
if isfield(S2,'n_sim_emp')
    n_sim_candidates(end+1) = double(S2.n_sim_emp);
end
if isfield(S3,'n_sim_emp')
    n_sim_candidates(end+1) = double(S3.n_sim_emp);
end

if isempty(n_sim_candidates)
    n_sim_emp = 500;
else
    n_sim_emp = round(max(n_sim_candidates));
end

rng_seed = 41;
rng(rng_seed);

cmax = 1-dmax^2/(2*R^2);

N_draws_emp = zeros(n_sim_emp,1);
betti0_total_sim = zeros(n_sim_emp,1);
betti0_per_bin_sim = zeros(n_sim_emp,Nb);

for s = 1:n_sim_emp

    Ns = max(draw_poisson(N_mean_target),1);
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

    betti0_total_sim(s) = numel(component_sizes);

    latitudes = asin(max(min(points(:,3),1),-1));
    bin_id = discretize(latitudes,phi_edges);

    for c = 1:numel(component_sizes)
        members = find(labels == c);
        component_size = numel(members);

        % Attribution symetrique : chaque composante contribue au total
        % pour 1, repartie entre les latitudes de ses satellites.
        contribution = 1/component_size;

        for r = 1:component_size
            bin = bin_id(members(r));

            if ~isnan(bin)
                betti0_per_bin_sim(s,bin) = ...
                    betti0_per_bin_sim(s,bin)+contribution;
            end
        end
    end
end

N_emp = mean(N_draws_emp);

betti0_per_bin_emp = mean(betti0_per_bin_sim,1);
betti0_density_phi_emp = ...
    betti0_per_bin_emp ./ dphi_bins;

betti0_total_emp = mean(betti0_total_sim);

betti0_per_bin_emp_std = ...
    std(betti0_per_bin_sim,0,1);

betti0_per_bin_emp_sem = ...
    betti0_per_bin_emp_std/sqrt(n_sim_emp);

%% ============================================================
%  7. Contributions individuelles
%% ============================================================

N1_density_phi_th = N1_per_bin_th./dphi_bins;
N1_density_phi_emp = N1_per_bin_emp./dphi_bins;

N2_density_phi_th = N2_per_bin_th./dphi_bins;
N2_density_phi_emp = N2_per_bin_emp./dphi_bins;

N3_density_phi_th = N3_per_bin_th./dphi_bins;
N3_density_phi_emp = N3_per_bin_emp./dphi_bins;

beta_macro_density_phi = ...
    beta_macro_per_bin_th./dphi_bins;

%% ============================================================
%  8. Diagnostic
%% ============================================================

valid = isfinite(betti0_density_phi_th) ...
    & isfinite(betti0_density_phi_emp) ...
    & betti0_density_phi_th > 0;

ratio_emp_th_phi = nan(1,Nb);
ratio_emp_th_phi(valid) = ...
    betti0_density_phi_emp(valid) ...
    ./ betti0_density_phi_th(valid);

rmse_density = sqrt(mean( ...
    (betti0_density_phi_emp(valid) ...
    -betti0_density_phi_th(valid)).^2));

relative_error_total = ...
    abs(betti0_total_emp-betti0_total_th) ...
    / max(abs(betti0_total_th),eps);

%% ============================================================
%  9. Figures
%% ============================================================

figure;
plot(rad2deg(phi_vals),betti0_density_phi_th,'LineWidth',2);
hold on;
plot(rad2deg(phi_vals),betti0_density_phi_emp,'--','LineWidth',1.8);
grid on;
xlabel('Latitude \phi (deg)');
ylabel('Densite moyenne de composantes par radian');
title('\beta_0 local : theorie jusqu''aux trimeres et empirique total');
legend('Theorie jusqu''aux trimeres','Monte-Carlo total','Location','best');
hold off;

figure;
plot(rad2deg(phi_vals),betti0_per_bin_th,'LineWidth',2);
hold on;
plot(rad2deg(phi_vals),betti0_per_bin_emp,'--','LineWidth',1.8);
grid on;
xlabel('Latitude \phi (deg)');
ylabel('Nombre moyen de composantes dans la tranche');
title('Nombre moyen local de composantes connexes');
legend('Theorie jusqu''aux trimeres','Monte-Carlo total','Location','best');
hold off;

figure;
plot(rad2deg(phi_vals),ratio_emp_th_phi,'LineWidth',1.8);
hold on;
yline(1,'--','Accord parfait');
grid on;
xlabel('Latitude \phi (deg)');
ylabel('Empirique / theorie');
title('Qualite du modele local de \beta_0');
hold off;

figure;
plot(rad2deg(phi_vals),beta_macro_per_bin_th,'LineWidth',1.5);
hold on;
plot(rad2deg(phi_vals),N1_per_bin_th,'LineWidth',1.5);
plot(rad2deg(phi_vals),N2_per_bin_th,'LineWidth',1.5);
plot(rad2deg(phi_vals),N3_per_bin_th,'LineWidth',1.5);
grid on;
xlabel('Latitude \phi (deg)');
ylabel('Contribution theorique par tranche');
title('Decomposition locale du modele de \beta_0');
legend('Composantes macroscopiques','N_1','N_2','N_3', ...
    'Location','best');
hold off;

%% ============================================================
%  10. Affichage
%% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' BETTI_0(phi) : THEORIE TRIMERES / EMPIRIQUE TOTAL\n');
fprintf('============================================================\n');
fprintf('N moyen retenu                     : %.8f\n',N);
fprintf('Inclinaison                        : %.8f deg\n',rad2deg(inc));
fprintf('Nombre de composantes macros       : %d\n', ...
    n_macro_components);
fprintf('------------------------------------------------------------\n');
fprintf('N1 total theorique                 : %.10f\n', ...
    sum(N1_per_bin_th));
fprintf('N2 total theorique                 : %.10f\n', ...
    sum(N2_per_bin_th));
fprintf('N3 total theorique                 : %.10f\n', ...
    sum(N3_per_bin_th));
fprintf('beta0 total theorique              : %.10f\n', ...
    betti0_total_th);
fprintf('------------------------------------------------------------\n');
fprintf('N moyen tire empiriquement         : %.10f\n', ...
    N_emp);
fprintf('Nombre de simulations empiriques  : %d\n', ...
    n_sim_emp);
fprintf('beta0 total Monte-Carlo direct    : %.10f\n', ...
    betti0_total_emp);
fprintf('------------------------------------------------------------\n');
fprintf('Erreur relative totale            : %.3e\n', ...
    relative_error_total);
fprintf('RMSE densite locale               : %.3e\n', ...
    rmse_density);
fprintf('============================================================\n');

%% ============================================================
%  11. Sauvegarde
%% ============================================================

output_file = fullfile(script_dir,'betti_phi_results.mat');

save(output_file, ...
    'R','inc','dmax','N','n_macro_components', ...
    'phi_vals','phi_edges','dphi_bins','phi_bin_probability', ...
    'beta_macro_per_bin_th','beta_macro_density_phi', ...
    'N1_per_bin_th','N1_per_bin_emp', ...
    'N2_per_bin_th','N2_per_bin_emp', ...
    'N3_per_bin_th','N3_per_bin_emp', ...
    'N1_density_phi_th','N1_density_phi_emp', ...
    'N2_density_phi_th','N2_density_phi_emp', ...
    'N3_density_phi_th','N3_density_phi_emp', ...
    'betti0_per_bin_th','betti0_per_bin_emp', ...
    'betti0_per_bin_sim','betti0_per_bin_emp_std', ...
    'betti0_per_bin_emp_sem', ...
    'betti0_density_phi_th','betti0_density_phi_emp', ...
    'betti0_total_th','betti0_total_emp','betti0_total_sim', ...
    'N_mean_target','N_draws_emp','N_emp','n_sim_emp','rng_seed', ...
    'ratio_emp_th_phi','rmse_density','relative_error_total', ...
    'file_N1','file_N2','file_N3');

fprintf('Resultats sauvegardes dans %s\n',output_file);

%% ============================================================
%  Fonctions locales
%% ============================================================

function path_out = find_result_file(script_dir,file_name)
    candidates = {
        fullfile(script_dir,file_name)
        fullfile(script_dir,'..',file_name)
    };

    path_out = '';

    for k = 1:numel(candidates)
        if isfile(candidates{k})
            path_out = candidates{k};
            return;
        end
    end

    error('Fichier %s introuvable.',file_name);
end

function value = get_scalar_field(S,names)
    for k = 1:numel(names)
        if isfield(S,names{k})
            value = double(S.(names{k}));
            if isscalar(value)
                return;
            end
        end
    end

    error('Aucun champ scalaire parmi : %s.',strjoin(names,', '));
end

function require_fields(S,names,file_name)
    for k = 1:numel(names)
        if ~isfield(S,names{k})
            error('Le fichier %s doit contenir %s.',file_name,names{k});
        end
    end
end

function x = row_vector(x)
    x = double(x(:).');
end

function assert_close_relative(a,b,tol,message)
    scale = max([abs(a),abs(b),1]);
    if abs(a-b) > tol*scale
        error('%s : %.10g contre %.10g.',message,a,b);
    end
end

function [edges,N1_th,N1_emp] = extract_N1_per_bin(S)
    % Format recent directement compatible
    if isfield(S,'phi_edges')
        edges = row_vector(S.phi_edges);
    elseif isfield(S,'lat_edges')
        edges = row_vector(S.lat_edges);
    else
        error('N1_phi_results.mat ne contient ni phi_edges ni lat_edges.');
    end

    if isfield(S,'N1_per_bin_th')
        N1_th = row_vector(S.N1_per_bin_th);

    elseif isfield(S,'p_isolated_lat_th') && ...
           isfield(S,'N') && isfield(S,'n_realizations') && ...
           isfield(S,'nodes_per_bin')

        % Probabilite theorique d'isolement appliquee au nombre moyen
        % de satellites observe dans chaque tranche.
        mean_nodes_per_bin = ...
            row_vector(S.nodes_per_bin)/double(S.n_realizations);

        N1_th = mean_nodes_per_bin ...
            .* row_vector(S.p_isolated_lat_th);

    else
        error(['Impossible de reconstruire N1 theorique par tranche ', ...
               'depuis N1_phi_results.mat.']);
    end

    if isfield(S,'N1_per_bin_emp')
        N1_emp = row_vector(S.N1_per_bin_emp);

    elseif isfield(S,'isolated_sum_bin') && isfield(S,'n_realizations')
        N1_emp = row_vector(S.isolated_sum_bin) ...
            / double(S.n_realizations);

    elseif isfield(S,'nodes_per_bin') && ...
           isfield(S,'p_isolated_lat_emp') && ...
           isfield(S,'n_realizations')

        mean_nodes_per_bin = ...
            row_vector(S.nodes_per_bin)/double(S.n_realizations);

        N1_emp = mean_nodes_per_bin ...
            .* row_vector(S.p_isolated_lat_emp);

    else
        error(['Impossible de reconstruire N1 empirique par tranche ', ...
               'depuis N1_phi_results.mat.']);
    end

    if numel(edges) ~= numel(N1_th)+1 || ...
       numel(edges) ~= numel(N1_emp)+1
        error('Dimensions incoherentes dans N1_phi_results.mat.');
    end
end

function target_counts = rebin_counts(source_edges,source_counts,target_edges)
    % Redistribution conservatrice en supposant une densite uniforme
    % a l'interieur de chaque tranche source.

    source_edges = row_vector(source_edges);
    source_counts = row_vector(source_counts);
    target_edges = row_vector(target_edges);

    target_counts = zeros(1,numel(target_edges)-1);

    for s = 1:numel(source_counts)
        a = source_edges(s);
        b = source_edges(s+1);
        width = b-a;

        if width <= 0
            error('Les bords de tranche source ne sont pas croissants.');
        end

        for t = 1:numel(target_counts)
            c = target_edges(t);
            d = target_edges(t+1);

            overlap = max(0,min(b,d)-max(a,c));

            if overlap > 0
                target_counts(t) = target_counts(t) ...
                    + source_counts(s)*overlap/width;
            end
        end
    end
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
