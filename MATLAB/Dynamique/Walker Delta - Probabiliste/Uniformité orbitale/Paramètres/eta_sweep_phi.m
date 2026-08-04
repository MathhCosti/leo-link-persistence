%% eta_sweep_phi_corrige.m
% Comparaison locale du facteur eta_sweep(phi)
% Walker Delta a uniformite orbitale.
%
% Ce script reutilise directement la densite locale deja calculee dans
% densite_phi_results.mat :
%
%   eta_sweep^th(phi_b)
%      = exp[-lambda_theory_bins(phi_b) A_cap(dmax)]
%
%   eta_sweep^{lambda_emp}(phi_b)
%      = exp[-lambda_empirical(phi_b) A_cap(dmax)]
%
% Une mesure directe de eta_sweep est egalement realisee a partir de
% analysis_temp_results.mat en echantillonnant la zone nouvellement
% balayee par les satellites.
%
% Le script mesure aussi la vraie probabilite empirique
%
%   p_bridge,bord^emp(phi)
%     = P(lien pont a t | lien rompu entre t et t+dt, phi),
%
% en identifiant exactement les ponts du graphe G_t, puis en examinant
% les liens presents a t et absents a t+dt.
%
% Entrees :
%   densite_phi_results.mat
%   analysis_temp_results.mat
%
% Sortie :
%   eta_sweep_phi_results.mat

clear; clc; close all;

%% ============================================================
%  1. Parametres numeriques
%% ============================================================

rng_seed = 51;
rng(rng_seed);

% Nombre maximal d'intervalles temporels analyses
n_time_samples = 800;

% Nombre maximal de satellites testes par tranche et par instant
n_sat_per_bin_per_time = 8;

% Nombre de points tires dans la calotte centree en r_i(t+dt)
n_points_cap = 1200;

% Nombre minimal de points dans le croissant pour accepter la mesure
min_crescent_points = 20;

%% ============================================================
%  2. Chargement
%% ============================================================

script_dir = fileparts(mfilename('fullpath'));

density_candidates = {
    fullfile(script_dir,'densite_phi_results.mat')
    fullfile(script_dir,'..','densite_phi_results.mat')
    fullfile(script_dir,'..','Valeurs locales','densite_phi_results.mat')
};

analysis_candidates = {
    fullfile(script_dir,'analysis_temp_results.mat')
    fullfile(script_dir,'..','analysis_temp_results.mat')
};

density_file = find_first_existing(density_candidates);
analysis_file = find_first_existing(analysis_candidates);

if isempty(density_file)
    error('Fichier densite_phi_results.mat introuvable.');
end

if isempty(analysis_file)
    error('Fichier analysis_temp_results.mat introuvable.');
end

D = load(density_file);
A = load(analysis_file);

required_density = { ...
    'R','inc','N','phi_edges','phi_centers', ...
    'area_bins','lambda_theory_bins','lambda_empirical'};

for k = 1:numel(required_density)
    if ~isfield(D,required_density{k})
        error('densite_phi_results.mat doit contenir %s.', ...
            required_density{k});
    end
end

required_analysis = {'Positions','Adjacency','dmax'};

for k = 1:numel(required_analysis)
    if ~isfield(A,required_analysis{k})
        error('analysis_temp_results.mat doit contenir %s.', ...
            required_analysis{k});
    end
end

R = double(D.R);
inc = double(D.inc);
N_density = double(D.N);

phi_edges = double(D.phi_edges(:).');
phi_vals = double(D.phi_centers(:).');
area_bins = double(D.area_bins(:).');

lambda_phi_th = double(D.lambda_theory_bins(:).');
lambda_phi_emp_density = double(D.lambda_empirical(:).');

dmax = double(A.dmax);
Positions = A.Positions;
Adjacency = A.Adjacency;

Nt = min(numel(Positions),numel(Adjacency));
Nb = numel(phi_vals);

if numel(phi_edges) ~= Nb+1
    error('Dimensions incoherentes entre phi_edges et phi_centers.');
end

if numel(lambda_phi_th) ~= Nb || ...
   numel(lambda_phi_emp_density) ~= Nb
    error('Dimensions incoherentes dans densite_phi_results.mat.');
end

dphi_bins = diff(phi_edges);

if Nt < 2
    error('Il faut au moins deux instants temporels.');
end

% Verification des parametres communs avec l'analyse temporelle
N_analysis = size(Positions{1},1);

if abs(N_analysis-N_density) > 0
    warning(['N differe entre densite_phi_results.mat et ', ...
             'analysis_temp_results.mat : %d contre %d.'], ...
             round(N_density),N_analysis);
end

R_analysis = mean(vecnorm(double(Positions{1}),2,2));

if abs(R_analysis-R) > 1e-6*max(R,1)
    warning(['R differe entre les deux fichiers : %.6f contre ', ...
             '%.6f km.'],R,R_analysis);
end

%% ============================================================
%  3. Modeles de eta_sweep issus des densites deja validees
%% ============================================================

A_intersection = ...
    (2*pi/3-sqrt(3)/2)*dmax^2;

eta_sweep_phi_th = ...
    exp(-lambda_phi_th*A_intersection);

eta_sweep_phi_from_empirical_density = ...
    exp(-lambda_phi_emp_density*A_intersection);

%% ============================================================
%  4. Choix des instants
%% ============================================================

n_time_samples = min(n_time_samples,Nt-1);
time_indices = unique(round(linspace(1,Nt-1,n_time_samples)));

%% ============================================================
%  5. Accumulateurs pour la mesure directe
%% ============================================================

new_area_points_total = zeros(1,Nb);
new_area_points_uncovered = zeros(1,Nb);

eta_measurements = cell(1,Nb);
n_satellite_tests = zeros(1,Nb);

% Accumulateurs pour la vraie probabilite qu'un lien rompu soit un pont.
broken_links_total = zeros(1,Nb);
broken_bridges_total = zeros(1,Nb);
broken_no_common_neighbor_total = zeros(1,Nb);
broken_bridge_and_no_common_total = zeros(1,Nb);

alpha_max = 2*asin(min(dmax/(2*R),1));
cmax = cos(alpha_max);

%% ============================================================
%  6. Mesure empirique directe de eta_sweep(phi)
%% ============================================================

for it = 1:numel(time_indices)

    t = time_indices(it);

    X0 = double(Positions{t});
    X1 = double(Positions{t+1});

    U0 = X0./vecnorm(X0,2,2);
    U1 = X1./vecnorm(X1,2,2);

    adjacency_t = logical(spones(Adjacency{t}));
    adjacency_t = adjacency_t | adjacency_t.';
    adjacency_t(1:size(adjacency_t,1)+1:end) = false;

    adjacency_t1 = logical(spones(Adjacency{t+1}));
    adjacency_t1 = adjacency_t1 | adjacency_t1.';
    adjacency_t1(1:size(adjacency_t1,1)+1:end) = false;

    labels = conncomp(graph(adjacency_t));

    latitude = asin(max(min(U0(:,3),1),-1));

    % Vraie probabilite P(pont | rupture, phi)
    bridge_matrix_t = find_bridges_sparse(adjacency_t);
    broken_upper = triu(adjacency_t & ~adjacency_t1,1);
    [broken_i,broken_j] = find(broken_upper);

    if ~isempty(broken_i)
        link_midpoint = U0(broken_i,:) + U0(broken_j,:);
        midpoint_norm = vecnorm(link_midpoint,2,2);
        valid_midpoint = midpoint_norm > 1e-12;
        link_midpoint(valid_midpoint,:) = link_midpoint(valid_midpoint,:) ./ midpoint_norm(valid_midpoint);

        link_latitude = nan(numel(broken_i),1);
        link_latitude(valid_midpoint) = asin(max(min(link_midpoint(valid_midpoint,3),1),-1));
        link_bin = discretize(link_latitude,phi_edges);

        linear_idx = sub2ind(size(adjacency_t),broken_i,broken_j);
        is_bridge_broken = full(bridge_matrix_t(linear_idx));
        is_no_common_neighbor = false(numel(broken_i),1);

        for e = 1:numel(broken_i)
            u = broken_i(e);
            v = broken_j(e);
            is_no_common_neighbor(e) = ~any(adjacency_t(u,:) & adjacency_t(v,:));
        end

        for b_bridge = 1:Nb
            in_bin_bridge = link_bin == b_bridge;
            if ~any(in_bin_bridge), continue; end

            broken_links_total(b_bridge) = broken_links_total(b_bridge) + nnz(in_bin_bridge);
            broken_bridges_total(b_bridge) = broken_bridges_total(b_bridge) + nnz(is_bridge_broken(in_bin_bridge));
            broken_no_common_neighbor_total(b_bridge) = broken_no_common_neighbor_total(b_bridge) + nnz(is_no_common_neighbor(in_bin_bridge));
            broken_bridge_and_no_common_total(b_bridge) = broken_bridge_and_no_common_total(b_bridge) + nnz(is_bridge_broken(in_bin_bridge) & is_no_common_neighbor(in_bin_bridge));
        end
    end
    bin_id = discretize(latitude,phi_edges);

    for b = 1:Nb

        candidates = find(bin_id == b);

        if isempty(candidates)
            continue;
        end

        if numel(candidates) > n_sat_per_bin_per_time
            candidates = candidates(randperm( ...
                numel(candidates),n_sat_per_bin_per_time));
        end

        for idx = 1:numel(candidates)

            i = candidates(idx);

            angular_displacement = acos(max(min( ...
                dot(U0(i,:),U1(i,:)),1),-1));

            if angular_displacement < 1e-12
                continue;
            end

            % Points uniformes dans B_i(t+dt)
            cap_points = sample_spherical_cap( ...
                U1(i,:),alpha_max,n_points_cap);

            % Zone nouvelle : B_i(t+dt) \ B_i(t)
            in_old_cap = cap_points*U0(i,:).' >= cmax;
            crescent_points = cap_points(~in_old_cap,:);

            n_crescent = size(crescent_points,1);

            if n_crescent < min_crescent_points
                continue;
            end

            component_members = find(labels == labels(i));
            component_members(component_members == i) = [];

            if isempty(component_members)
                is_covered = false(n_crescent,1);
            else
                component_positions = U0(component_members,:);

                is_covered = any( ...
                    crescent_points*component_positions.' >= cmax, ...
                    2);
            end

            n_uncovered = nnz(~is_covered);
            eta_local = n_uncovered/n_crescent;

            new_area_points_total(b) = ...
                new_area_points_total(b)+n_crescent;

            new_area_points_uncovered(b) = ...
                new_area_points_uncovered(b)+n_uncovered;

            eta_measurements{b}(end+1,1) = eta_local; %#ok<SAGROW>
            n_satellite_tests(b) = n_satellite_tests(b)+1;
        end
    end
end

%% ============================================================
%  7. Moyennes empiriques directes
%% ============================================================

% Moyenne ponderee par le nombre de points du croissant
eta_sweep_phi_emp_direct = safe_divide( ...
    new_area_points_uncovered,new_area_points_total);

% Moyenne donnant le meme poids a chaque satellite teste
eta_sweep_phi_emp_mean_per_sat = nan(1,Nb);
eta_sweep_phi_emp_std = nan(1,Nb);
eta_sweep_phi_emp_sem = nan(1,Nb);

for b = 1:Nb
    values = eta_measurements{b};

    if isempty(values)
        continue;
    end

    eta_sweep_phi_emp_mean_per_sat(b) = mean(values);
    eta_sweep_phi_emp_std(b) = std(values);
    eta_sweep_phi_emp_sem(b) = ...
        eta_sweep_phi_emp_std(b)/sqrt(numel(values));
end

%% ============================================================
%  7.b Vraie probabilite empirique qu'un lien rompu soit un pont
%% ============================================================

p_bridge_bord_phi_emp_true = safe_divide(broken_bridges_total,broken_links_total);
p_no_common_neighbor_given_break_phi = safe_divide(broken_no_common_neighbor_total,broken_links_total);
p_true_bridge_given_no_common_break_phi = safe_divide(broken_bridge_and_no_common_total,broken_no_common_neighbor_total);

p_bridge_bord_emp_true_global = sum(broken_bridges_total)/max(sum(broken_links_total),1);
p_no_common_neighbor_given_break_global = sum(broken_no_common_neighbor_total)/max(sum(broken_links_total),1);
p_true_bridge_given_no_common_break_global = sum(broken_bridge_and_no_common_total)/max(sum(broken_no_common_neighbor_total),1);

%% ============================================================
%  8. Diagnostics
%% ============================================================

valid_direct = isfinite(eta_sweep_phi_emp_direct) ...
    & isfinite(eta_sweep_phi_th) ...
    & new_area_points_total > 0;

ratio_direct_th_phi = nan(1,Nb);
ratio_direct_th_phi(valid_direct) = ...
    eta_sweep_phi_emp_direct(valid_direct) ...
    ./ eta_sweep_phi_th(valid_direct);

ratio_density_emp_th_phi = ...
    eta_sweep_phi_from_empirical_density ...
    ./ eta_sweep_phi_th;

rmse_direct_vs_th = sqrt(mean( ...
    (eta_sweep_phi_emp_direct(valid_direct) ...
    -eta_sweep_phi_th(valid_direct)).^2));

rmse_direct_vs_emp_density = sqrt(mean( ...
    (eta_sweep_phi_emp_direct(valid_direct) ...
    -eta_sweep_phi_from_empirical_density(valid_direct)).^2));

mean_abs_error_direct_vs_th = mean(abs( ...
    eta_sweep_phi_emp_direct(valid_direct) ...
    -eta_sweep_phi_th(valid_direct)));

valid_bridge = isfinite(p_bridge_bord_phi_emp_true) & broken_links_total > 0;
ratio_true_bridge_to_eta_phi = nan(1,Nb);
ratio_true_bridge_to_eta_phi(valid_bridge) = p_bridge_bord_phi_emp_true(valid_bridge) ./ eta_sweep_phi_emp_direct(valid_bridge);

rmse_true_bridge_vs_eta = sqrt(mean((p_bridge_bord_phi_emp_true(valid_bridge)-eta_sweep_phi_emp_direct(valid_bridge)).^2,'omitnan'));
rmse_true_bridge_vs_th = sqrt(mean((p_bridge_bord_phi_emp_true(valid_bridge)-eta_sweep_phi_th(valid_bridge)).^2,'omitnan'));

%% ============================================================
%  9. Figures
%% ============================================================

figure;
plot(rad2deg(phi_vals),eta_sweep_phi_th,'LineWidth',2);
hold on;
plot(rad2deg(phi_vals), ...
    eta_sweep_phi_from_empirical_density,'--','LineWidth',1.8);
plot(rad2deg(phi_vals), ...
    eta_sweep_phi_emp_direct,':','LineWidth',1.8);
plot(rad2deg(phi_vals), ...
    p_bridge_bord_phi_emp_true,'-.','LineWidth',2.0);
grid on;
xlabel('Latitude \phi (deg)');
ylabel('\eta_{sweep}^\Delta(\phi)');
title('Facteur local de redondance spatiale');
legend('Modele avec \lambda theorique', ...
       'Modele avec \lambda empirique', ...
       'Mesure directe de \eta_{sweep}', ...
       'Vraie P(pont | rupture)', ...
       'Location','best');
ylim([0,1]);
hold off;

figure;
plot(rad2deg(phi_vals),eta_sweep_phi_emp_mean_per_sat, ...
    'LineWidth',1.8);
hold on;
plot(rad2deg(phi_vals),eta_sweep_phi_emp_direct, ...
    '--','LineWidth',1.8);
grid on;
xlabel('Latitude \phi (deg)');
ylabel('\eta_{sweep}^\Delta(\phi)');
title('Ponderation de la mesure directe');
legend('Moyenne par satellite', ...
       'Moyenne ponderee par aire echantillonnee', ...
       'Location','best');
ylim([0,1]);
hold off;

figure;
plot(rad2deg(phi_vals),ratio_direct_th_phi,'LineWidth',1.8);
hold on;
plot(rad2deg(phi_vals),ratio_density_emp_th_phi, ...
    '--','LineWidth',1.8);
yline(1,':','Accord parfait');
grid on;
xlabel('Latitude \phi (deg)');
ylabel('Rapport');
title('Qualite du modele local de \eta_{sweep}');
legend('\eta direct / \eta(\lambda theorique)', ...
       '\eta(\lambda empirique) / \eta(\lambda theorique)', ...
       'Location','best');
hold off;

figure;
hold on;
plot(rad2deg(phi_vals),eta_sweep_phi_emp_direct,'LineWidth',2,'DisplayName','\eta_{sweep}^{emp}');
plot(rad2deg(phi_vals),p_no_common_neighbor_given_break_phi,'--','LineWidth',1.8,'DisplayName','P(pas de voisin commun | rupture)');
plot(rad2deg(phi_vals),p_bridge_bord_phi_emp_true,'-.','LineWidth',2,'DisplayName','P(pont | rupture)');
grid on;
xlabel('Latitude \phi (deg)');
ylabel('Probabilite');
title('Approximation locale et vraie probabilite de pont');
legend('Location','best');
ylim([0,1]);
hold off;

figure;
plot(rad2deg(phi_vals),p_true_bridge_given_no_common_break_phi,'LineWidth',2);
grid on;
xlabel('Latitude \phi (deg)');
ylabel('P(pont | pas de voisin commun, rupture)');
title('Effet des chemins alternatifs plus longs');
ylim([0,1]);

figure;
yyaxis left
plot(rad2deg(phi_vals),broken_links_total,'LineWidth',1.8);
ylabel('Nombre de liens rompus');
yyaxis right
plot(rad2deg(phi_vals),ratio_true_bridge_to_eta_phi,'--','LineWidth',1.8);
ylabel('P(pont | rupture) / \eta_{sweep}^{emp}');
grid on;
xlabel('Latitude \phi (deg)');
title('Echantillonnage des ruptures et correction de pont');

figure;
plot(rad2deg(phi_vals),new_area_points_total,'LineWidth',1.8);
grid on;
xlabel('Latitude \phi (deg)');
ylabel('Nombre de points du croissant');
title('Echantillonnage direct selon la latitude');

%% ============================================================
%  10. Affichage
%% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' ETA_SWEEP(phi) - DENSITE REUTILISEE / MESURE DIRECTE\n');
fprintf('============================================================\n');
fprintf('Fichier de densite                  : %s\n',density_file);
fprintf('Fichier temporel                    : %s\n',analysis_file);
fprintf('N du fichier de densite             : %.8f\n',N_density);
fprintf('N du fichier temporel               : %d\n',N_analysis);
fprintf('Inclinaison                         : %.8f deg\n',rad2deg(inc));
fprintf('dmax                                : %.8f km\n',dmax);
fprintf('Nombre d''instants analyses          : %d\n', ...
    numel(time_indices));
fprintf('------------------------------------------------------------\n');
fprintf('A_intersection(dmax)                : %.10f km^2\n', ...
    A_intersection);
fprintf('RMSE direct / lambda theorique      : %.3e\n', ...
    rmse_direct_vs_th);
fprintf('RMSE direct / lambda empirique      : %.3e\n', ...
    rmse_direct_vs_emp_density);
fprintf('Erreur absolue moyenne directe/th   : %.3e\n', ...
    mean_abs_error_direct_vs_th);
fprintf('------------------------------------------------------------\n');
fprintf('Nombre total de liens rompus        : %d\n',round(sum(broken_links_total)));
fprintf('Nombre de ponts parmi les ruptures  : %d\n',round(sum(broken_bridges_total)));
fprintf('P vraie(pont | rupture) globale     : %.10f\n',p_bridge_bord_emp_true_global);
fprintf('P(pas voisin commun | rupture)      : %.10f\n',p_no_common_neighbor_given_break_global);
fprintf('P(pont | pas voisin commun, rupture): %.10f\n',p_true_bridge_given_no_common_break_global);
fprintf('RMSE vraie P(pont)/eta direct       : %.3e\n',rmse_true_bridge_vs_eta);
fprintf('RMSE vraie P(pont)/modele theorique : %.3e\n',rmse_true_bridge_vs_th);
fprintf('============================================================\n');

%% ============================================================
%  11. Sauvegarde
%% ============================================================

output_file = fullfile(script_dir,'eta_sweep_phi_results.mat');

save(output_file, ...
    'density_file','analysis_file', ...
    'R','inc','dmax','N_density','N_analysis','A_intersection', ...
    'phi_vals','phi_edges','dphi_bins','area_bins', ...
    'lambda_phi_th','lambda_phi_emp_density', ...
    'eta_sweep_phi_th', ...
    'eta_sweep_phi_from_empirical_density', ...
    'eta_sweep_phi_emp_direct', ...
    'eta_sweep_phi_emp_mean_per_sat', ...
    'eta_sweep_phi_emp_std','eta_sweep_phi_emp_sem', ...
    'eta_measurements','new_area_points_total', ...
    'new_area_points_uncovered','n_satellite_tests', ...
    'time_indices','n_time_samples','n_sat_per_bin_per_time', ...
    'n_points_cap','min_crescent_points','rng_seed', ...
    'ratio_direct_th_phi','ratio_density_emp_th_phi', ...
    'rmse_direct_vs_th','rmse_direct_vs_emp_density', ...
    'mean_abs_error_direct_vs_th', ...
    'broken_links_total','broken_bridges_total', ...
    'broken_no_common_neighbor_total','broken_bridge_and_no_common_total', ...
    'p_bridge_bord_phi_emp_true', ...
    'p_no_common_neighbor_given_break_phi', ...
    'p_true_bridge_given_no_common_break_phi', ...
    'p_bridge_bord_emp_true_global', ...
    'p_no_common_neighbor_given_break_global', ...
    'p_true_bridge_given_no_common_break_global', ...
    'ratio_true_bridge_to_eta_phi', ...
    'rmse_true_bridge_vs_eta','rmse_true_bridge_vs_th');

fprintf('Resultats sauvegardes dans %s\n',output_file);

%% ============================================================
%  Fonctions locales
%% ============================================================

function path_out = find_first_existing(candidates)
    path_out = '';

    for k = 1:numel(candidates)
        if isfile(candidates{k})
            path_out = candidates{k};
            return;
        end
    end
end

function points = sample_spherical_cap(center,alpha,n)
    center = center/norm(center);

    cos_gamma = 1-(1-cos(alpha))*rand(n,1);
    sin_gamma = sqrt(max(1-cos_gamma.^2,0));
    theta = 2*pi*rand(n,1);

    if abs(center(3)) < 0.9
        reference = [0,0,1];
    else
        reference = [1,0,0];
    end

    e1 = cross(reference,center);
    e1 = e1/norm(e1);

    e2 = cross(center,e1);
    e2 = e2/norm(e2);

    points = ...
        cos_gamma.*center ...
        + sin_gamma.*cos(theta).*e1 ...
        + sin_gamma.*sin(theta).*e2;

    points = points./vecnorm(points,2,2);
end

function bridge_matrix = find_bridges_sparse(adjacency)
% Algorithme de Tarjan : ponts exacts du graphe en O(N+E).
    adjacency = logical(spones(adjacency));
    adjacency = adjacency | adjacency.';
    n = size(adjacency,1);
    adjacency(1:n+1:end) = false;

    discovery = zeros(n,1);
    low = zeros(n,1);
    parent = zeros(n,1);
    visited = false(n,1);
    timer = 0;

    max_edges = nnz(triu(adjacency,1));
    bridge_i = zeros(max_edges,1);
    bridge_j = zeros(max_edges,1);
    n_bridges = 0;

    for root = 1:n
        if ~visited(root)
            dfs(root);
        end
    end

    bridge_i = bridge_i(1:n_bridges);
    bridge_j = bridge_j(1:n_bridges);
    bridge_matrix = sparse([bridge_i;bridge_j],[bridge_j;bridge_i],true,n,n);

    function dfs(u)
        visited(u) = true;
        timer = timer+1;
        discovery(u) = timer;
        low(u) = timer;
        neighbors = find(adjacency(u,:));

        for kk = 1:numel(neighbors)
            v = neighbors(kk);
            if ~visited(v)
                parent(v) = u;
                dfs(v);
                low(u) = min(low(u),low(v));
                if low(v) > discovery(u)
                    n_bridges = n_bridges+1;
                    bridge_i(n_bridges) = u;
                    bridge_j(n_bridges) = v;
                end
            elseif v ~= parent(u)
                low(u) = min(low(u),discovery(v));
            end
        end
    end
end

function ratio = safe_divide(num,den)
    ratio = nan(size(num));
    valid = den > 0;
    ratio(valid) = num(valid)./den(valid);
end
