%% plink_t_phi_compare_resimulation.m
% Comparaison theorie / empirique de
%
%   p_link(t,phi)
%     = P(d_12(t) <= dmax | Phi_1(t)=phi),
%
% en resimulant directement plusieurs realisations independantes du
% Walker Delta a uniformite spatiale initiale.
%
% Aucun fichier analysis_temp_results.mat n'est charge.
%
% Initialisation spatiale :
%   - sin(phi_0) uniforme sur [-sin(i),sin(i)] ;
%   - branche ascendante/descendante choisie equiprobablement ;
%   - RAAN uniforme sur [0,2pi).
%
% Cette construction est equivalente a la loi initiale de phase
%
%   f_{u_0}(u) = |cos(u)|/4.
%
% Dynamique :
%   u(t) = u_0 + omega*t.
%
% Position :
%
%   r(u,Omega) = R [
%       cos(Omega)cos(u)-sin(Omega)sin(u)cos(i)
%       sin(Omega)cos(u)+cos(Omega)sin(u)cos(i)
%       sin(u)sin(i)
%   ].
%
% Sortie :
%   plink_t_phi_results.mat

clear; clc; close all;

% Dossier contenant ce script. Cette variable est utilisee plus bas
% pour sauvegarder plink_t_phi_results.mat au meme endroit que le code.
script_dir = fileparts(mfilename('fullpath'));

% Cas de securite : si le script est lance d'une maniere pour laquelle
% mfilename('fullpath') est vide, on utilise le dossier courant.
if isempty(script_dir)
    script_dir = pwd;
end

%% ============================================================
%  1. Parametres physiques et numeriques
%% ============================================================

% Parametres de la constellation.
%
% N est recupere directement depuis analysis_temp_results.mat afin
% d'utiliser exactement la meme realisation que le code principal.
analysis_file = fullfile(script_dir,'analysis_temp_results.mat');

if ~isfile(analysis_file)
    analysis_file = fullfile(script_dir,'..','analysis_temp_results.mat');
end

if ~isfile(analysis_file)
    error('Fichier analysis_temp_results.mat introuvable.');
end

Sanalysis = load(analysis_file,'N');

if ~isfield(Sanalysis,'N')
    error('La variable N est absente de %s.',analysis_file);
end

N = double(Sanalysis.N);

if ~isscalar(N) || ~isfinite(N) || N < 2
    error('La valeur N chargee depuis analysis_temp_results.mat est invalide.');
end

N = round(N);

R = 6921;                % km
dmax = 1500;             % km
inc = deg2rad(58);       % rad

mu = 398600;             % km^3/s^2
omega = sqrt(mu/R^3);    % rad/s

% Grille temporelle.
dt = 20;                 % s
Tmax = 12000;            % s
time_values = 0:dt:Tmax;
Nt = numel(time_values);

% Nombre de realisations independantes.
n_iterations = 100;

% Grilles de latitude.
n_phi_bins_emp = 15;
n_phi_bins_th = 300;

% Quadrature theorique.
n_u_quad = 720;

% Seuil minimal total de satellites dans une tranche, toutes
% realisations confondues.
min_satellites_per_bin_total = 10;

% Nombre d'instants affiches.
n_selected_times = 5;

% Monte-Carlo conditionnel direct, utilise uniquement aux instants
% affiches pour verifier independamment la formule analytique.
n_mc_conditional = 10000;

% Graine reproductible.
rng_seed = 1;
rng(rng_seed);

%% ============================================================
%  2. Grilles de latitude
%% ============================================================

phi_edges_emp = linspace(-inc,inc,n_phi_bins_emp+1);
phi_vals_emp = ...
    0.5*(phi_edges_emp(1:end-1)+phi_edges_emp(2:end));
dphi_emp = diff(phi_edges_emp);

phi_edges_th = linspace(-inc,inc,n_phi_bins_th+1);
phi_vals_th = ...
    0.5*(phi_edges_th(1:end-1)+phi_edges_th(2:end));
dphi_th = diff(phi_edges_th);

%% ============================================================
%  3. Simulation de plusieurs realisations
%% ============================================================

p_link_emp_iterations = ...
    nan(n_iterations,Nt,n_phi_bins_emp);

satellite_count_emp_iterations = ...
    zeros(n_iterations,Nt,n_phi_bins_emp);

degree_sum_emp_iterations = ...
    zeros(n_iterations,Nt,n_phi_bins_emp);

p_link_emp_global_iterations = ...
    nan(n_iterations,Nt);

for r = 1:n_iterations

    %% 3.a Tirage initial uniforme en surface dans la bande
    sin_phi0 = ...
        -sin(inc) + 2*sin(inc)*rand(N,1);

    phi0 = asin(sin_phi0);

    % Relation orbitale :
    %
    %   sin(phi_0) = sin(i) sin(u_0).
    %
    % Il faut donc inverser cette relation avec
    %
    %   sin(u_0) = sin(phi_0)/sin(i).
    %
    % L'ancienne version utilisait directement
    % u_base = asin(sin_phi0), ce qui comprimait artificiellement
    % les latitudes vers l'equateur.
    sin_u0 = sin_phi0/sin(inc);
    sin_u0 = min(max(sin_u0,-1),1);

    % Une latitude correspond a deux branches orbitales.
    u_base = asin(sin_u0);

    ascending = rand(N,1) < 0.5;

    u0 = zeros(N,1);
    u0(ascending) = u_base(ascending);
    u0(~ascending) = pi-u_base(~ascending);
    u0 = mod(u0,2*pi);

    % Verification numerique de la latitude reconstruite.
    phi0_reconstructed = asin(max(min( ...
        sin(inc)*sin(u0),1),-1));

    reconstruction_error = max(abs( ...
        phi0_reconstructed-phi0));

    if reconstruction_error > 1e-12
        error(['Erreur dans l''inversion latitude-phase : ', ...
               'max |phi_reconstruite-phi0| = %.3e rad.'], ...
              reconstruction_error);
    end

    % RAAN independants et uniformes.
    Omega = 2*pi*rand(N,1);

    cos_Omega = cos(Omega);
    sin_Omega = sin(Omega);

    %% 3.b Dynamique temporelle
    for t_idx = 1:Nt

        t = time_values(t_idx);
        u = mod(u0+omega*t,2*pi);

        cos_u = cos(u);
        sin_u = sin(u);

        % Positions cartesiennes.
        X = R * [ ...
            cos_Omega.*cos_u ...
                - sin_Omega.*sin_u*cos(inc), ...
            sin_Omega.*cos_u ...
                + cos_Omega.*sin_u*cos(inc), ...
            sin_u*sin(inc)];

        % Latitude instantanee.
        latitude = asin(max(min(X(:,3)/R,1),-1));

        % Matrice des distances euclidiennes.
        gram = X*X.';
        squared_norm = sum(X.^2,2);

        D2 = ...
            squared_norm ...
            + squared_norm.' ...
            - 2*gram;

        D2 = max(D2,0);

        adjacency = ...
            D2 <= dmax^2;

        adjacency(1:N+1:end) = false;
        adjacency = triu(adjacency,1);
        adjacency = adjacency | adjacency.';

        degree = sum(adjacency,2);
        bin_id = discretize(latitude,phi_edges_emp);

        total_degree = 0;
        total_satellites = 0;

        for b = 1:n_phi_bins_emp
            mask = bin_id == b;
            n_b = nnz(mask);

            satellite_count_emp_iterations(r,t_idx,b) = n_b;

            if n_b == 0
                continue;
            end

            degree_sum_b = sum(degree(mask));

            degree_sum_emp_iterations(r,t_idx,b) = degree_sum_b;

            p_link_emp_iterations(r,t_idx,b) = ...
                degree_sum_b/(n_b*(N-1));

            total_degree = total_degree+degree_sum_b;
            total_satellites = total_satellites+n_b;
        end

        if total_satellites > 0
            p_link_emp_global_iterations(r,t_idx) = ...
                total_degree/(total_satellites*(N-1));
        end
    end

    fprintf('Realisation %d/%d terminee.\n',r,n_iterations);
end

%% ============================================================
%  4. Moyennes empiriques entre realisations
%% ============================================================

satellite_count_emp_total = ...
    squeeze(sum(satellite_count_emp_iterations,1));

degree_sum_emp_total = ...
    squeeze(sum(degree_sum_emp_iterations,1));

p_link_emp_mean = nan(Nt,n_phi_bins_emp);

valid_emp_total = ...
    satellite_count_emp_total >= min_satellites_per_bin_total;

p_link_emp_mean(valid_emp_total) = ...
    degree_sum_emp_total(valid_emp_total) ...
    ./ ((N-1)*satellite_count_emp_total(valid_emp_total));

% Moyenne simple, ecart-type et SEM entre realisations.
p_link_emp_iteration_mean = ...
    squeeze(mean(p_link_emp_iterations,1,'omitnan'));

p_link_emp_iteration_std = ...
    squeeze(std(p_link_emp_iterations,0,1,'omitnan'));

n_valid_iterations_per_bin = ...
    squeeze(sum(isfinite(p_link_emp_iterations),1));

p_link_emp_iteration_sem = ...
    p_link_emp_iteration_std ...
    ./ sqrt(max(n_valid_iterations_per_bin,1));

% Valeur globale empirique.
p_link_emp_global_mean = ...
    mean(p_link_emp_global_iterations,1,'omitnan').';

p_link_emp_global_std = ...
    std(p_link_emp_global_iterations,0,1,'omitnan').';

p_link_emp_global_sem = ...
    p_link_emp_global_std/sqrt(n_iterations);

%% ============================================================
%  5. Noyau geometrique theorique
%% ============================================================

alpha_max = 2*asin(min(dmax/(2*R),1));
cos_alpha_max = cos(alpha_max);

u2_grid = (0:n_u_quad-1)*(2*pi/n_u_quad);
du = 2*pi/n_u_quad;

s_th = sin(phi_vals_th)/sin(inc);
s_th = min(max(s_th,-1),1);

u_plus_th = mod(asin(s_th),2*pi);
u_minus_th = mod(pi-asin(s_th),2*pi);

G_plus_th = zeros(n_phi_bins_th,n_u_quad);
G_minus_th = zeros(n_phi_bins_th,n_u_quad);

for b = 1:n_phi_bins_th
    G_plus_th(b,:) = G_link_kernel( ...
        u_plus_th(b),u2_grid,inc,cos_alpha_max);

    G_minus_th(b,:) = G_link_kernel( ...
        u_minus_th(b),u2_grid,inc,cos_alpha_max);
end

%% ============================================================
%  6. Theorie fine de p_link(t,phi)
%% ============================================================

p_link_th_fine = nan(Nt,n_phi_bins_th);
branch_weight_plus_th = nan(Nt,n_phi_bins_th);
branch_weight_minus_th = nan(Nt,n_phi_bins_th);
f_phi_th_fine = nan(Nt,n_phi_bins_th);

for t_idx = 1:Nt

    t = time_values(t_idx);

    f_u2 = abs(cos(u2_grid-omega*t))/4;
    f_u2 = f_u2/(sum(f_u2)*du);

    f_plus = abs(cos(u_plus_th-omega*t))/4;
    f_minus = abs(cos(u_minus_th-omega*t))/4;

    denom_branch = f_plus+f_minus;

    w_plus = zeros(1,n_phi_bins_th);
    w_minus = zeros(1,n_phi_bins_th);

    valid_branch = denom_branch > 1e-14;

    w_plus(valid_branch) = ...
        f_plus(valid_branch)./denom_branch(valid_branch);

    w_minus(valid_branch) = ...
        f_minus(valid_branch)./denom_branch(valid_branch);

    w_plus(~valid_branch) = 0.5;
    w_minus(~valid_branch) = 0.5;

    branch_weight_plus_th(t_idx,:) = w_plus;
    branch_weight_minus_th(t_idx,:) = w_minus;

    p_plus = (G_plus_th*f_u2(:))*du;
    p_minus = (G_minus_th*f_u2(:))*du;

    p_link_th_fine(t_idx,:) = ...
        w_plus.*p_plus.' + w_minus.*p_minus.';

    jacobian_abs = ...
        sin(inc)*abs(cos(u_plus_th))./cos(phi_vals_th);

    valid_jacobian = jacobian_abs > 1e-14;

    f_phi = zeros(1,n_phi_bins_th);

    f_phi(valid_jacobian) = ...
        (f_plus(valid_jacobian)+f_minus(valid_jacobian)) ...
        ./ jacobian_abs(valid_jacobian);

    mass_phi = sum(f_phi.*dphi_th);

    if mass_phi > 0
        f_phi = f_phi/mass_phi;
    end

    f_phi_th_fine(t_idx,:) = f_phi;
end

p_link_th_fine = min(max(p_link_th_fine,0),1);

%% ============================================================
%  7. Theorie moyennee dans les tranches empiriques
%
% L'empirique correspond a une moyenne conditionnelle sur toute la
% tranche de latitude. La quantite theorique comparable est donc
%
%   p_link,b^th(t)
%     = int_b p_link^th(t,phi) f_Phi^th(t,phi) dphi
%       ---------------------------------------------
%              int_b f_Phi^th(t,phi) dphi.
%
% On ne compare plus l'empirique a la seule valeur theorique au centre.
%% ============================================================

p_link_th_on_emp = nan(Nt,n_phi_bins_emp);
f_phi_mass_on_emp = zeros(Nt,n_phi_bins_emp);
f_phi_th_on_emp = zeros(Nt,n_phi_bins_emp);
satellite_count_th_on_emp = zeros(Nt,n_phi_bins_emp);

for t_idx = 1:Nt
    for b = 1:n_phi_bins_emp

        in_bin = ...
            phi_vals_th >= phi_edges_emp(b) ...
            & phi_vals_th < phi_edges_emp(b+1);

        if b == n_phi_bins_emp
            in_bin = ...
                phi_vals_th >= phi_edges_emp(b) ...
                & phi_vals_th <= phi_edges_emp(b+1);
        end

        if ~any(in_bin)
            continue;
        end

        local_mass = ...
            f_phi_th_fine(t_idx,in_bin) ...
            .* dphi_th(in_bin);

        mass_bin = sum(local_mass);

        f_phi_mass_on_emp(t_idx,b) = mass_bin;
        f_phi_th_on_emp(t_idx,b) = mass_bin/dphi_emp(b);
        satellite_count_th_on_emp(t_idx,b) = N*mass_bin;

        if mass_bin > 0
            p_link_th_on_emp(t_idx,b) = ...
                sum( ...
                    p_link_th_fine(t_idx,in_bin) ...
                    .* local_mass) ...
                / mass_bin;
        end
    end
end

p_link_th_on_emp = min(max(p_link_th_on_emp,0),1);

% Valeur au centre conservee seulement comme diagnostic.
p_link_th_at_centers = nan(Nt,n_phi_bins_emp);

for t_idx = 1:Nt
    p_link_th_at_centers(t_idx,:) = interp1( ...
        phi_vals_th,p_link_th_fine(t_idx,:), ...
        phi_vals_emp,'pchip','extrap');
end

p_link_th_at_centers = min(max(p_link_th_at_centers,0),1);

%% ============================================================
%  8. Probabilite globale theorique
%% ============================================================

p_link_th_global = ...
    sum(p_link_th_fine.*f_phi_th_fine.*dphi_th,2);

p_link_th_global_direct = nan(Nt,1);

for t_idx = 1:Nt

    t = time_values(t_idx);

    f_u = abs(cos(u2_grid-omega*t))/4;
    f_u = f_u/(sum(f_u)*du);

    integral_value = 0;

    for j = 1:n_u_quad
        G_row = G_link_kernel( ...
            u2_grid(j),u2_grid,inc,cos_alpha_max);

        integral_value = integral_value ...
            + f_u(j)*sum(G_row.*f_u)*du^2;
    end

    p_link_th_global_direct(t_idx) = integral_value;
end


%% ============================================================
%  8.b Monte-Carlo conditionnel direct
%
% Ce calcul teste la formule locale sans simuler une constellation
% complete. Pour chaque tranche et chaque instant affiche :
%
%   1) on tire la latitude du premier satellite dans la tranche selon
%      f_Phi(phi,t) ;
%   2) on tire sa branche orbitale conditionnelle ;
%   3) on tire u_2 selon f_u(u_2,t) ;
%   4) on tire DeltaOmega uniformement ;
%   5) on teste directement la condition d <= dmax.
%
% Cette mesure doit coincider avec la theorie moyennee dans la tranche.
%% ============================================================

selected_indices = unique(round( ...
    linspace(1,Nt,n_selected_times)));

p_link_mc_conditional = nan(numel(selected_indices),n_phi_bins_emp);
p_link_mc_conditional_sem = nan(numel(selected_indices),n_phi_bins_emp);

for k_sel = 1:numel(selected_indices)

    t_idx = selected_indices(k_sel);
    fprintf('Monte-Carlo conditionnel : instant %d/%d\n', ...
        k_sel,numel(selected_indices));
    t = time_values(t_idx);

    for b = 1:n_phi_bins_emp

        % Tirage de phi dans la tranche selon f_Phi(phi,t), en utilisant
        % la grille theorique fine comme loi discrete.
        in_bin = ...
            phi_vals_th >= phi_edges_emp(b) ...
            & phi_vals_th < phi_edges_emp(b+1);

        if b == n_phi_bins_emp
            in_bin = ...
                phi_vals_th >= phi_edges_emp(b) ...
                & phi_vals_th <= phi_edges_emp(b+1);
        end

        idx_candidates = find(in_bin);

        if isempty(idx_candidates)
            continue;
        end

        weights_phi = ...
            f_phi_th_fine(t_idx,idx_candidates) ...
            .* dphi_th(idx_candidates);

        if sum(weights_phi) <= 0
            continue;
        end

        % Suppression des poids nuls afin d'obtenir une CDF strictement
        % croissante, compatible avec discretize.
        positive_weight = weights_phi > 0;
        idx_candidates = idx_candidates(positive_weight);
        weights_phi = weights_phi(positive_weight);

        if isempty(idx_candidates)
            continue;
        end

        weights_phi = weights_phi/sum(weights_phi);
        cdf_edges = [0,cumsum(weights_phi)];
        cdf_edges(end) = 1;

        draws_phi = rand(n_mc_conditional,1);

        % Tirage vectorise dans la loi discrete de latitude.
        local_index = discretize(draws_phi,cdf_edges);

        % IMPORTANT : forcer explicitement un vecteur colonne.
        % Sans cela, l'indexation d'un vecteur ligne peut produire un
        % vecteur ligne et les operations avec u2_draw (colonne) creent
        % par expansion implicite une matrice n_mc x n_mc.
        phi_draw = phi_vals_th(idx_candidates(local_index));
        phi_draw = phi_draw(:);

        % Deux phases compatibles avec chaque latitude tiree.
        s_draw = sin(phi_draw)/sin(inc);
        s_draw = min(max(s_draw,-1),1);

        u_plus_draw = mod(asin(s_draw),2*pi);
        u_minus_draw = mod(pi-asin(s_draw),2*pi);

        u_plus_draw = u_plus_draw(:);
        u_minus_draw = u_minus_draw(:);

        f_plus_draw = abs(cos(u_plus_draw-omega*t))/4;
        f_minus_draw = abs(cos(u_minus_draw-omega*t))/4;

        prob_plus = ...
            f_plus_draw./max(f_plus_draw+f_minus_draw,eps);

        choose_plus = rand(n_mc_conditional,1) < prob_plus;

        u1_draw = u_minus_draw;
        u1_draw(choose_plus) = u_plus_draw(choose_plus);

        % Tirage de u2 selon f_u(u2,t) par rejet.
        u2_draw = sample_phase_fu(n_mc_conditional,omega*t);

        % Difference de RAAN uniforme.
        delta_Omega = 2*pi*rand(n_mc_conditional,1);

        % Verification de dimensions avant les produits terme a terme.
        if ~isequal(size(u1_draw),size(u2_draw),size(delta_Omega))
            error(['Dimensions Monte-Carlo incompatibles : ', ...
                   'u1=%s, u2=%s, DeltaOmega=%s.'], ...
                  mat2str(size(u1_draw)), ...
                  mat2str(size(u2_draw)), ...
                  mat2str(size(delta_Omega)));
        end

        % Produit scalaire normalise.
        A_mc = ...
            cos(u1_draw).*cos(u2_draw) ...
            + cos(inc)^2 ...
            .* sin(u1_draw).*sin(u2_draw);

        B_mc = ...
            cos(inc) ...
            .* (sin(u1_draw).*cos(u2_draw) ...
            - cos(u1_draw).*sin(u2_draw));

        C_mc = ...
            sin(inc)^2 ...
            .* sin(u1_draw).*sin(u2_draw);

        dot_normalized = ...
            A_mc.*cos(delta_Omega) ...
            + B_mc.*sin(delta_Omega) ...
            + C_mc;

        is_link = dot_normalized >= cos_alpha_max;

        p_link_mc_conditional(k_sel,b) = mean(is_link);
        p_link_mc_conditional_sem(k_sel,b) = ...
            sqrt(mean(is_link)*(1-mean(is_link)) ...
            / n_mc_conditional);
    end
end

%% ============================================================
%  9. Diagnostics
%% ============================================================

valid_compare = ...
    isfinite(p_link_emp_mean) ...
    & isfinite(p_link_th_on_emp) ...
    & satellite_count_emp_total >= min_satellites_per_bin_total;

difference = nan(size(p_link_th_on_emp));

difference(valid_compare) = ...
    p_link_emp_mean(valid_compare) ...
    - p_link_th_on_emp(valid_compare);

rmse_global_grid = sqrt(mean( ...
    difference(valid_compare).^2));

mae_global_grid = mean(abs( ...
    difference(valid_compare)));

bias_global_grid = mean( ...
    difference(valid_compare));

weights_emp = satellite_count_emp_total;
valid_weighted = valid_compare & weights_emp > 0;

bias_weighted = ...
    sum(weights_emp(valid_weighted).*difference(valid_weighted)) ...
    / sum(weights_emp(valid_weighted));

rmse_weighted = sqrt( ...
    sum(weights_emp(valid_weighted).*difference(valid_weighted).^2) ...
    / sum(weights_emp(valid_weighted)));

rmse_by_phi = nan(1,n_phi_bins_emp);
mae_by_phi = nan(1,n_phi_bins_emp);
bias_by_phi = nan(1,n_phi_bins_emp);

for b = 1:n_phi_bins_emp
    valid_b = valid_compare(:,b);

    if ~any(valid_b)
        continue;
    end

    err_b = difference(valid_b,b);

    rmse_by_phi(b) = sqrt(mean(err_b.^2));
    mae_by_phi(b) = mean(abs(err_b));
    bias_by_phi(b) = mean(err_b);
end

global_theory_consistency_error = max(abs( ...
    p_link_th_global-p_link_th_global_direct));

% Accord entre Monte-Carlo conditionnel et theorie moyennee.
p_link_th_selected = p_link_th_on_emp(selected_indices,:);

valid_mc = ...
    isfinite(p_link_mc_conditional) ...
    & isfinite(p_link_th_selected);

rmse_mc_vs_theory = sqrt(mean( ...
    (p_link_mc_conditional(valid_mc) ...
    - p_link_th_selected(valid_mc)).^2));

bias_mc_vs_theory = mean( ...
    p_link_mc_conditional(valid_mc) ...
    - p_link_th_selected(valid_mc));

%% ============================================================
%  10. Figures
%% ============================================================

figure;
imagesc(time_values,rad2deg(phi_vals_th),p_link_th_fine.');
axis xy;
colorbar;
xlabel('Temps (s)');
ylabel('Latitude \phi (deg)');
title('p_{link}^{th}(t,\phi)');

figure;
imagesc(time_values,rad2deg(phi_vals_emp),p_link_emp_mean.');
axis xy;
colorbar;
xlabel('Temps (s)');
ylabel('Latitude \phi (deg)');
title(sprintf( ...
    'p_{link}^{emp}(t,\phi) moyen sur %d realisations', ...
    n_iterations));

figure;
imagesc(time_values,rad2deg(phi_vals_emp),difference.');
axis xy;
colorbar;
xlabel('Temps (s)');
ylabel('Latitude \phi (deg)');
title('Ecart p_{link}^{emp,moy}-p_{link}^{th}');

selected_indices = unique(round( ...
    linspace(1,Nt,n_selected_times)));

figure;
tiledlayout(numel(selected_indices),1, ...
    'TileSpacing','compact','Padding','compact');

for k = 1:numel(selected_indices)

    t_idx = selected_indices(k);

    nexttile;
    hold on;

    plot(rad2deg(phi_vals_th), ...
        p_link_th_fine(t_idx,:), ...
        'LineWidth',2, ...
        'DisplayName','Theorie fine');

    plot(rad2deg(phi_vals_emp), ...
        p_link_th_on_emp(t_idx,:), ...
        's','LineWidth',1.1, ...
        'DisplayName','Theorie moyennee par tranche');

    errorbar( ...
        rad2deg(phi_vals_emp), ...
        p_link_mc_conditional(k,:), ...
        p_link_mc_conditional_sem(k,:), ...
        'd--','LineWidth',1.1, ...
        'DisplayName','Monte-Carlo conditionnel');

    errorbar( ...
        rad2deg(phi_vals_emp), ...
        p_link_emp_mean(t_idx,:), ...
        p_link_emp_iteration_sem(t_idx,:), ...
        'o-','LineWidth',1.2, ...
        'DisplayName','Empirique moyen \pm SEM');

    grid on;
    ylabel('p_{link}');
    title(sprintf('t = %.1f s',time_values(t_idx)));

    if k == 1
        legend('Location','best');
    end

    if k == numel(selected_indices)
        xlabel('Latitude \phi (deg)');
    end

    hold off;
end

figure;
hold on;
plot(rad2deg(phi_vals_emp), ...
    mean(p_link_th_on_emp,1,'omitnan'), ...
    'LineWidth',2, ...
    'DisplayName','Theorie moyennee dans les tranches');
plot(rad2deg(phi_vals_emp), ...
    mean(p_link_th_at_centers,1,'omitnan'), ...
    '--','LineWidth',1.8, ...
    'DisplayName','Theorie aux centres');
grid on;
xlabel('Latitude \phi (deg)');
ylabel('Moyenne temporelle de p_{link}');
title('Effet de la moyenne theorique dans les tranches');
legend('Location','best');
hold off;

figure;
hold on;

plot(time_values,p_link_th_global, ...
    'LineWidth',2, ...
    'DisplayName','Theorie');

plot(time_values,p_link_emp_global_mean, ...
    'LineWidth',1.5, ...
    'DisplayName','Empirique moyen');

plot(time_values, ...
    p_link_emp_global_mean+p_link_emp_global_sem, ...
    ':','LineWidth',1, ...
    'DisplayName','Empirique + SEM');

plot(time_values, ...
    p_link_emp_global_mean-p_link_emp_global_sem, ...
    ':','LineWidth',1, ...
    'DisplayName','Empirique - SEM');

grid on;
xlabel('Temps (s)');
ylabel('p_{link}(t)');
title('Probabilite globale de liaison');
legend('Location','best');
hold off;

figure;
hold on;

plot(rad2deg(phi_vals_emp),rmse_by_phi, ...
    'LineWidth',2,'DisplayName','RMSE');

plot(rad2deg(phi_vals_emp),mae_by_phi, ...
    '--','LineWidth',1.8,'DisplayName','MAE');

plot(rad2deg(phi_vals_emp),bias_by_phi, ...
    ':','LineWidth',1.8,'DisplayName','Biais emp-th');

grid on;
xlabel('Latitude \phi (deg)');
ylabel('Erreur');
title('Erreur locale apres resimulation');
legend('Location','best');
hold off;

%% ============================================================
%  11. Affichage console
%% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' p_link(t,phi) - RESIMULATION MULTI-REALISATIONS\n');
fprintf('============================================================\n');
fprintf('Nombre de realisations              : %d\n',n_iterations);
fprintf('Graine aleatoire                    : %d\n',rng_seed);
fprintf('Initialisation                      : sin(u0)=sin(phi0)/sin(i)\n');
fprintf('N                                   : %d\n',N);
fprintf('Source de N                         : %s\n',analysis_file);
fprintf('R                                   : %.8f km\n',R);
fprintf('Inclinaison                         : %.8f deg\n',rad2deg(inc));
fprintf('dmax                                : %.8f km\n',dmax);
fprintf('dt                                  : %.8f s\n',dt);
fprintf('Tmax                                : %.8f s\n',Tmax);
fprintf('Tranches empiriques                 : %d\n',n_phi_bins_emp);
fprintf('Tranches theoriques                 : %d\n',n_phi_bins_th);
fprintf('Tirages MC conditionnels / tranche  : %d\n', ...
    n_mc_conditional);
fprintf('------------------------------------------------------------\n');
fprintf('RMSE non ponderee                   : %.10e\n', ...
    rmse_global_grid);
fprintf('MAE non ponderee                    : %.10e\n', ...
    mae_global_grid);
fprintf('Biais non pondere emp-theorie       : %.10e\n', ...
    bias_global_grid);
fprintf('RMSE ponderee                       : %.10e\n', ...
    rmse_weighted);
fprintf('Biais pondere                       : %.10e\n', ...
    bias_weighted);
fprintf('Coherence theorie globale max       : %.10e\n', ...
    global_theory_consistency_error);
fprintf('RMSE MC conditionnel / theorie      : %.10e\n', ...
    rmse_mc_vs_theory);
fprintf('Biais MC conditionnel - theorie     : %.10e\n', ...
    bias_mc_vs_theory);
fprintf('------------------------------------------------------------\n');
fprintf('Moyenne p_link globale empirique    : %.10e\n', ...
    mean(p_link_emp_global_mean,'omitnan'));
fprintf('Moyenne p_link globale theorique    : %.10e\n', ...
    mean(p_link_th_global,'omitnan'));
fprintf('============================================================\n');

%% ============================================================
%  12. Sauvegarde
%% ============================================================

output_file = fullfile(script_dir,'plink_t_phi_results.mat');

% Alias compatibles avec les codes aval.
p_link_th = p_link_th_on_emp;
f_phi_th = f_phi_th_on_emp;
phi_vals = phi_vals_emp;
dphi = dphi_emp;
phi_edges = phi_edges_emp;
n_phi_bins = n_phi_bins_emp;

p_link_emp = p_link_emp_mean;
p_link_emp_sampled = p_link_emp_mean;
satellite_count_emp = satellite_count_emp_total;
satellite_count_emp_sampled = satellite_count_emp_total;
degree_sum_emp = degree_sum_emp_total;

p_link_emp_global = p_link_emp_global_mean;
p_link_emp_global_sampled = p_link_emp_global_mean;

time_theory = time_values;
time_indices = 1:Nt;

save(output_file, ...
    'N','analysis_file','R','dmax','inc','mu','omega','alpha_max', ...
    'dt','Tmax','time_values','time_theory','time_indices', ...
    'n_iterations','rng_seed', ...
    'phi_edges_emp','phi_vals_emp','dphi_emp', ...
    'phi_edges_th','phi_vals_th','dphi_th', ...
    'n_phi_bins_emp','n_phi_bins_th', ...
    'n_u_quad','u2_grid','du', ...
    'u_plus_th','u_minus_th', ...
    'G_plus_th','G_minus_th', ...
    'branch_weight_plus_th','branch_weight_minus_th', ...
    'f_phi_th_fine','f_phi_th_on_emp', ...
    'f_phi_mass_on_emp','satellite_count_th_on_emp', ...
    'p_link_th_fine','p_link_th_on_emp', ...
    'p_link_th_at_centers', ...
    'selected_indices', ...
    'n_mc_conditional', ...
    'p_link_mc_conditional', ...
    'p_link_mc_conditional_sem', ...
    'rmse_mc_vs_theory','bias_mc_vs_theory', ...
    'p_link_emp_iterations', ...
    'p_link_emp_iteration_mean', ...
    'p_link_emp_iteration_std', ...
    'p_link_emp_iteration_sem', ...
    'n_valid_iterations_per_bin', ...
    'p_link_emp_mean', ...
    'satellite_count_emp_iterations', ...
    'satellite_count_emp_total', ...
    'degree_sum_emp_iterations','degree_sum_emp_total', ...
    'p_link_emp_global_iterations', ...
    'p_link_emp_global_mean', ...
    'p_link_emp_global_std', ...
    'p_link_emp_global_sem', ...
    'p_link_th_global','p_link_th_global_direct', ...
    'difference','valid_compare', ...
    'rmse_global_grid','mae_global_grid','bias_global_grid', ...
    'rmse_weighted','bias_weighted', ...
    'rmse_by_phi','mae_by_phi','bias_by_phi', ...
    'global_theory_consistency_error', ...
    'min_satellites_per_bin_total', ...
    ... % Alias
    'p_link_th','f_phi_th','phi_vals','dphi', ...
    'phi_edges','n_phi_bins', ...
    'p_link_emp','p_link_emp_sampled', ...
    'satellite_count_emp','satellite_count_emp_sampled', ...
    'degree_sum_emp', ...
    'p_link_emp_global','p_link_emp_global_sampled');

fprintf('Resultats sauvegardes dans %s\n',output_file);

%% ============================================================
%  Fonction locale
%% ============================================================

function u = sample_phase_fu(n,phase_shift)
% Tirage selon f_u(u,t)=|cos(u-phase_shift)|/4 sur [0,2pi).
% Rejet avec proposition uniforme et enveloppe 1/4.

    u = zeros(n,1);
    n_done = 0;

    while n_done < n
        n_batch = max(1000,2*(n-n_done));

        candidate = 2*pi*rand(n_batch,1);
        accept = rand(n_batch,1) ...
            <= abs(cos(candidate-phase_shift));

        accepted = candidate(accept);
        n_take = min(numel(accepted),n-n_done);

        if n_take > 0
            u(n_done+(1:n_take)) = accepted(1:n_take);
            n_done = n_done+n_take;
        end
    end
end

function G = G_link_kernel(u1,u2,inc,cos_alpha_max)

    u2 = double(u2);

    A = cos(u1).*cos(u2) ...
        + cos(inc)^2.*sin(u1).*sin(u2);

    B = cos(inc).*( ...
        sin(u1).*cos(u2) ...
        - cos(u1).*sin(u2));

    C = sin(inc)^2.*sin(u1).*sin(u2);

    rho = sqrt(A.^2+B.^2);
    G = zeros(size(u2));

    degenerate = rho < 1e-14;

    if any(degenerate)
        G(degenerate) = ...
            double(C(degenerate) >= cos_alpha_max);
    end

    regular = ~degenerate;

    if any(regular)
        q = (cos_alpha_max-C(regular))./rho(regular);

        G_regular = zeros(size(q));
        G_regular(q <= -1) = 1;

        interior = q > -1 & q < 1;
        G_regular(interior) = acos(q(interior))/pi;

        G_regular(q >= 1) = 0;
        G(regular) = G_regular;
    end
end
