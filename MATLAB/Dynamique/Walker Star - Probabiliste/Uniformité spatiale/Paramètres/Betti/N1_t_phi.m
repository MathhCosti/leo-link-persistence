%% satellites_isoles_t_phi_compare.m
% Comparaison theorie / empirique du nombre local de satellites isoles
%
%   N_1(t,phi)
%
% dans chaque tranche de latitude.
%
% Un satellite est isole lorsque son degre est nul.
%
% IMPORTANT :
% conditionnellement a la latitude phi, le premier satellite peut etre
% situe sur l'une des deux branches orbitales u_+(phi) ou u_-(phi).
% Il ne faut donc pas utiliser directement
%
%   (1-p_link(t,phi))^(N-1),
%
% car p_link(t,phi) est deja moyenne sur les deux branches.
%
% La probabilite correcte d'isolement est
%
%   p_iso(t,phi)
%     = w_+(t,phi) [1-p_+(t,phi)]^(N-1)
%       + w_-(t,phi) [1-p_-(t,phi)]^(N-1),
%
% ou p_+ et p_- sont les probabilites de lien conditionnelles a chaque
% branche du premier satellite.
%
% Le nombre theorique d'isoles dans une tranche b vaut alors
%
%   N_1,b^th(t)
%     = N int_b f_Phi(t,phi) p_iso(t,phi) dphi.
%
% Empiriquement :
%
%   N_1,b^{emp,(r)}(t)
%     = nombre de satellites de degre nul dans la tranche b
%       pour la simulation r.
%
% Entree :
%   plink_t_phi_results.mat
%
% Sortie :
%   satellites_isoles_t_phi_results.mat

clear; clc; close all;

%% ============================================================
%  1. Chargement des parametres
%% ============================================================

script_dir = fileparts(mfilename('fullpath'));

if isempty(script_dir)
    script_dir = pwd;
end

input_file = fullfile(script_dir, '..', 'plink_t_phi_results.mat');

if ~isfile(input_file)
    error('Fichier introuvable : %s',input_file);
end

S = load(input_file);

required_fields = { ...
    'N','R','dmax','inc','omega', ...
    'time_theory', ...
    'phi_edges_emp','phi_vals_emp','dphi_emp'};

for k = 1:numel(required_fields)
    if ~isfield(S,required_fields{k})
        error('Le fichier %s doit contenir %s.', ...
            input_file,required_fields{k});
    end
end

N = double(S.N);
R = double(S.R);
dmax = double(S.dmax);
inc = double(S.inc);
omega = double(S.omega);

time_values_full = double(S.time_theory(:).');

phi_edges_emp = double(S.phi_edges_emp(:).');
phi_vals_emp = double(S.phi_vals_emp(:).');
dphi_emp = double(S.dphi_emp(:).');

Nb = numel(phi_vals_emp);

%% ============================================================
%  2. Parametres numeriques
%% ============================================================

% Nombre de simulations empiriques.
n_iterations = 100;

% Sous-echantillonnage temporel eventuel.
% Mettre 1 pour utiliser tous les instants.
time_stride = 1;

time_indices = 1:time_stride:numel(time_values_full);
time_values = time_values_full(time_indices);
Nt = numel(time_values);

% Grille theorique fine et quadrature en phase.
n_phi_bins_th = 300;
n_u_quad = 720;

% Nombre de coupes affichees.
n_selected_times = 5;

% Graine reproductible.
rng_seed = 3;
rng(rng_seed);

%% ============================================================
%  3. Grille theorique fine
%% ============================================================

phi_edges_th = linspace(-inc,inc,n_phi_bins_th+1);
phi_vals_th = ...
    0.5*(phi_edges_th(1:end-1)+phi_edges_th(2:end));
dphi_th = diff(phi_edges_th);

u2_grid = (0:n_u_quad-1)*(2*pi/n_u_quad);
du = 2*pi/n_u_quad;

alpha_max = 2*asin(min(dmax/(2*R),1));
cos_alpha_max = cos(alpha_max);

s_th = sin(phi_vals_th)/sin(inc);
s_th = min(max(s_th,-1),1);

u_plus_th = mod(asin(s_th),2*pi);
u_minus_th = mod(pi-asin(s_th),2*pi);

%% ============================================================
%  4. Noyaux de probabilite de lien pour chaque branche
%% ============================================================

G_plus = zeros(n_phi_bins_th,n_u_quad);
G_minus = zeros(n_phi_bins_th,n_u_quad);

for b = 1:n_phi_bins_th
    G_plus(b,:) = G_link_kernel( ...
        u_plus_th(b),u2_grid,inc,cos_alpha_max);

    G_minus(b,:) = G_link_kernel( ...
        u_minus_th(b),u2_grid,inc,cos_alpha_max);
end

%% ============================================================
%  5. Theorie fine de la probabilite d'isolement
%% ============================================================

p_link_plus_th = nan(Nt,n_phi_bins_th);
p_link_minus_th = nan(Nt,n_phi_bins_th);
p_link_th_fine = nan(Nt,n_phi_bins_th);

p_iso_th_fine = nan(Nt,n_phi_bins_th);
f_phi_th_fine = nan(Nt,n_phi_bins_th);

branch_weight_plus = nan(Nt,n_phi_bins_th);
branch_weight_minus = nan(Nt,n_phi_bins_th);

for t_idx = 1:Nt

    t = time_values(t_idx);

    % Loi de phase du second satellite.
    f_u2 = abs(cos(u2_grid-omega*t))/4;
    f_u2 = f_u2/(sum(f_u2)*du);

    % Poids des deux branches du premier satellite.
    f_plus = abs(cos(u_plus_th-omega*t))/4;
    f_minus = abs(cos(u_minus_th-omega*t))/4;

    denominator = f_plus+f_minus;

    w_plus = zeros(1,n_phi_bins_th);
    w_minus = zeros(1,n_phi_bins_th);

    valid_branch = denominator > 1e-14;

    w_plus(valid_branch) = ...
        f_plus(valid_branch)./denominator(valid_branch);

    w_minus(valid_branch) = ...
        f_minus(valid_branch)./denominator(valid_branch);

    w_plus(~valid_branch) = 0.5;
    w_minus(~valid_branch) = 0.5;

    branch_weight_plus(t_idx,:) = w_plus;
    branch_weight_minus(t_idx,:) = w_minus;

    % Probabilite de lien conditionnelle a chaque branche.
    p_plus = (G_plus*f_u2(:))*du;
    p_minus = (G_minus*f_u2(:))*du;

    p_plus = min(max(p_plus.',0),1);
    p_minus = min(max(p_minus.',0),1);

    p_link_plus_th(t_idx,:) = p_plus;
    p_link_minus_th(t_idx,:) = p_minus;

    p_link_th_fine(t_idx,:) = ...
        w_plus.*p_plus+w_minus.*p_minus;

    % Probabilite correcte d'isolement, en conditionnant d'abord
    % sur la branche du premier satellite.
    p_iso_th_fine(t_idx,:) = ...
        w_plus.*(1-p_plus).^(N-1) ...
        + w_minus.*(1-p_minus).^(N-1);

    % Loi de latitude.
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

%% ============================================================
%  6. Nombre theorique d'isoles par tranche empirique
%% ============================================================

isolated_count_th = zeros(Nt,Nb);
p_iso_th_on_emp = nan(Nt,Nb);
satellite_count_th = zeros(Nt,Nb);

for t_idx = 1:Nt
    for b = 1:Nb

        in_bin = ...
            phi_vals_th >= phi_edges_emp(b) ...
            & phi_vals_th < phi_edges_emp(b+1);

        if b == Nb
            in_bin = ...
                phi_vals_th >= phi_edges_emp(b) ...
                & phi_vals_th <= phi_edges_emp(b+1);
        end

        if ~any(in_bin)
            continue;
        end

        satellite_mass = ...
            f_phi_th_fine(t_idx,in_bin) ...
            .* dphi_th(in_bin);

        mass_bin = sum(satellite_mass);

        satellite_count_th(t_idx,b) = N*mass_bin;

        isolated_mass = sum( ...
            satellite_mass ...
            .* p_iso_th_fine(t_idx,in_bin));

        isolated_count_th(t_idx,b) = N*isolated_mass;

        if mass_bin > 0
            p_iso_th_on_emp(t_idx,b) = ...
                isolated_mass/mass_bin;
        end
    end
end

%% ============================================================
%  7. Resimulation empirique
%% ============================================================

isolated_count_emp_iterations = ...
    zeros(n_iterations,Nt,Nb);

satellite_count_emp_iterations = ...
    zeros(n_iterations,Nt,Nb);

p_iso_emp_iterations = ...
    nan(n_iterations,Nt,Nb);

for r = 1:n_iterations

    % Uniformite spatiale initiale.
    sin_phi0 = ...
        -sin(inc)+2*sin(inc)*rand(N,1);

    phi0 = asin(sin_phi0);

    sin_u0 = sin_phi0/sin(inc);
    sin_u0 = min(max(sin_u0,-1),1);

    u_base = asin(sin_u0);

    ascending = rand(N,1) < 0.5;

    u0 = zeros(N,1);
    u0(ascending) = u_base(ascending);
    u0(~ascending) = pi-u_base(~ascending);
    u0 = mod(u0,2*pi);

    Omega = 2*pi*rand(N,1);

    cos_Omega = cos(Omega);
    sin_Omega = sin(Omega);

    for t_idx = 1:Nt

        t = time_values(t_idx);
        u = mod(u0+omega*t,2*pi);

        cos_u = cos(u);
        sin_u = sin(u);

        X = R * [ ...
            cos_Omega.*cos_u ...
                - sin_Omega.*sin_u*cos(inc), ...
            sin_Omega.*cos_u ...
                + cos_Omega.*sin_u*cos(inc), ...
            sin_u*sin(inc)];

        latitude = asin(max(min(X(:,3)/R,1),-1));

        gram = X*X.';
        squared_norm = sum(X.^2,2);

        D2 = ...
            squared_norm+squared_norm.'-2*gram;

        D2 = max(D2,0);

        adjacency = D2 <= dmax^2;
        adjacency(1:N+1:end) = false;
        adjacency = adjacency | adjacency.';

        degree = sum(adjacency,2);
        is_isolated = degree == 0;

        bin_id = discretize(latitude,phi_edges_emp);

        for b = 1:Nb

            mask = bin_id == b;
            n_b = nnz(mask);
            n_iso_b = nnz(mask & is_isolated);

            satellite_count_emp_iterations(r,t_idx,b) = n_b;
            isolated_count_emp_iterations(r,t_idx,b) = n_iso_b;

            if n_b > 0
                p_iso_emp_iterations(r,t_idx,b) = ...
                    n_iso_b/n_b;
            end
        end
    end

    fprintf('Realisation %d/%d terminee.\n',r,n_iterations);
end

%% ============================================================
%  8. Moyennes empiriques
%% ============================================================

isolated_count_emp = squeeze(mean( ...
    isolated_count_emp_iterations,1,'omitnan'));

isolated_count_emp_std = squeeze(std( ...
    isolated_count_emp_iterations,0,1,'omitnan'));

isolated_count_emp_sem = ...
    isolated_count_emp_std/sqrt(n_iterations);

satellite_count_emp_total = squeeze(sum( ...
    satellite_count_emp_iterations,1));

isolated_count_emp_total = squeeze(sum( ...
    isolated_count_emp_iterations,1));

p_iso_emp = nan(Nt,Nb);

valid_satellite_count = satellite_count_emp_total > 0;

p_iso_emp(valid_satellite_count) = ...
    isolated_count_emp_total(valid_satellite_count) ...
    ./ satellite_count_emp_total(valid_satellite_count);

p_iso_emp_iteration_mean = squeeze(mean( ...
    p_iso_emp_iterations,1,'omitnan'));

p_iso_emp_std = squeeze(std( ...
    p_iso_emp_iterations,0,1,'omitnan'));

n_valid_iterations = squeeze(sum( ...
    isfinite(p_iso_emp_iterations),1));

p_iso_emp_sem = ...
    p_iso_emp_std ...
    ./ sqrt(max(n_valid_iterations,1));

%% ============================================================
%  9. Totaux globaux
%% ============================================================

isolated_total_th = sum( ...
    isolated_count_th,2,'omitnan');

isolated_total_emp_iterations = squeeze(sum( ...
    isolated_count_emp_iterations,3,'omitnan'));

isolated_total_emp = mean( ...
    isolated_total_emp_iterations,1,'omitnan').';

isolated_total_emp_std = std( ...
    isolated_total_emp_iterations,0,1,'omitnan').';

isolated_total_emp_sem = ...
    isolated_total_emp_std/sqrt(n_iterations);

%% ============================================================
%  10. Diagnostics
%% ============================================================

valid_compare = ...
    isfinite(isolated_count_th) ...
    & isfinite(isolated_count_emp);

difference = nan(Nt,Nb);

difference(valid_compare) = ...
    isolated_count_emp(valid_compare) ...
    - isolated_count_th(valid_compare);

rmse_grid = sqrt(mean( ...
    difference(valid_compare).^2));

mae_grid = mean(abs( ...
    difference(valid_compare)));

bias_grid = mean( ...
    difference(valid_compare));

rmse_by_phi = nan(1,Nb);
mae_by_phi = nan(1,Nb);
bias_by_phi = nan(1,Nb);

for b = 1:Nb
    valid_b = valid_compare(:,b);

    if ~any(valid_b)
        continue;
    end

    err_b = difference(valid_b,b);

    rmse_by_phi(b) = sqrt(mean(err_b.^2));
    mae_by_phi(b) = mean(abs(err_b));
    bias_by_phi(b) = mean(err_b);
end

%% ============================================================
%  11. Figures
%% ============================================================

figure;
imagesc( ...
    time_values, ...
    rad2deg(phi_vals_emp), ...
    isolated_count_th.');
axis xy;
colorbar;
xlabel('Temps (s)');
ylabel('Latitude \phi (deg)');
title('Nombre theorique de satellites isoles par tranche');

figure;
imagesc( ...
    time_values, ...
    rad2deg(phi_vals_emp), ...
    isolated_count_emp.');
axis xy;
colorbar;
xlabel('Temps (s)');
ylabel('Latitude \phi (deg)');
title(sprintf( ...
    'Nombre empirique moyen d''isoles sur %d simulations', ...
    n_iterations));

figure;
imagesc( ...
    time_values, ...
    rad2deg(phi_vals_emp), ...
    difference.');
axis xy;
colorbar;
xlabel('Temps (s)');
ylabel('Latitude \phi (deg)');
title('Ecart N_1^{emp}(t,\phi)-N_1^{th}(t,\phi)');

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

    % Equivalent continu dans une tranche de largeur moyenne.
    isolated_density_phi = ...
        N*f_phi_th_fine(it,:) ...
        .* p_iso_th_fine(it,:);

    continuous_count_equivalent = ...
        isolated_density_phi*mean(dphi_emp);

    plot(rad2deg(phi_vals_th), ...
        continuous_count_equivalent, ...
        'LineWidth',2, ...
        'DisplayName','Theorie fine (equivalent tranche)');

    plot(rad2deg(phi_vals_emp), ...
        isolated_count_th(it,:), ...
        's', ...
        'LineWidth',1.1, ...
        'DisplayName','Theorie integree par tranche');

    errorbar(rad2deg(phi_vals_emp), ...
        isolated_count_emp(it,:), ...
        isolated_count_emp_sem(it,:), ...
        'o-', ...
        'LineWidth',1.2, ...
        'DisplayName','Empirique moyen \pm SEM');

    grid on;
    ylabel('Nombre d''isoles');
    title(sprintf('t = %.1f s',time_values(it)));

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

plot(time_values,isolated_total_th, ...
    'LineWidth',2, ...
    'DisplayName','Theorie');

plot(time_values,isolated_total_emp, ...
    'LineWidth',1.5, ...
    'DisplayName','Empirique moyen');

plot(time_values, ...
    isolated_total_emp+isolated_total_emp_sem, ...
    ':','LineWidth',1, ...
    'DisplayName','Empirique + SEM');

plot(time_values, ...
    isolated_total_emp-isolated_total_emp_sem, ...
    ':','LineWidth',1, ...
    'DisplayName','Empirique - SEM');

grid on;
xlabel('Temps (s)');
ylabel('Nombre total de satellites isoles');
title('Nombre total de satellites isoles');
legend('Location','best');
hold off;

figure;
hold on;

plot(rad2deg(phi_vals_emp),rmse_by_phi, ...
    'LineWidth',2, ...
    'DisplayName','RMSE');

plot(rad2deg(phi_vals_emp),mae_by_phi, ...
    '--','LineWidth',1.8, ...
    'DisplayName','MAE');

plot(rad2deg(phi_vals_emp),bias_by_phi, ...
    ':','LineWidth',1.8, ...
    'DisplayName','Biais emp-th');

grid on;
xlabel('Latitude \phi (deg)');
ylabel('Erreur en nombre d''isoles');
title('Erreur de N_1(t,\phi) selon la latitude');
legend('Location','best');
hold off;

%% ============================================================
%  12. Affichage console
%% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' SATELLITES ISOLES N_1(t,phi) - THEORIE / EMPIRIQUE\n');
fprintf('============================================================\n');
fprintf('Fichier charge                    : %s\n',input_file);
fprintf('N                                 : %d\n',round(N));
fprintf('Nombre de simulations             : %d\n',n_iterations);
fprintf('Nombre d''instants                : %d\n',Nt);
fprintf('Nombre de tranches empiriques     : %d\n',Nb);
fprintf('------------------------------------------------------------\n');
fprintf('RMSE locale                       : %.10e satellites\n', ...
    rmse_grid);
fprintf('MAE locale                        : %.10e satellites\n', ...
    mae_grid);
fprintf('Biais moyen emp-theorie           : %.10e satellites\n', ...
    bias_grid);
fprintf('------------------------------------------------------------\n');
fprintf('Moyenne totale empirique          : %.10f satellites\n', ...
    mean(isolated_total_emp,'omitnan'));
fprintf('Moyenne totale theorique          : %.10f satellites\n', ...
    mean(isolated_total_th,'omitnan'));
fprintf('============================================================\n');

%% ============================================================
%  13. Sauvegarde
%% ============================================================

output_file = fullfile( ...
    script_dir, ...
    'N1_t_phi_results.mat');

save(output_file, ...
    'input_file', ...
    'N','R','dmax','inc','omega', ...
    'n_iterations','rng_seed', ...
    'time_stride','time_indices','time_values','Nt', ...
    'phi_edges_emp','phi_vals_emp','dphi_emp','Nb', ...
    'phi_edges_th','phi_vals_th','dphi_th', ...
    'n_phi_bins_th','n_u_quad','u2_grid','du', ...
    'alpha_max','cos_alpha_max', ...
    'u_plus_th','u_minus_th', ...
    'G_plus','G_minus', ...
    'branch_weight_plus','branch_weight_minus', ...
    'p_link_plus_th','p_link_minus_th', ...
    'p_link_th_fine', ...
    'p_iso_th_fine','p_iso_th_on_emp', ...
    'f_phi_th_fine', ...
    'satellite_count_th','isolated_count_th', ...
    'satellite_count_emp_iterations', ...
    'isolated_count_emp_iterations', ...
    'p_iso_emp_iterations', ...
    'isolated_count_emp', ...
    'isolated_count_emp_std', ...
    'isolated_count_emp_sem', ...
    'satellite_count_emp_total', ...
    'isolated_count_emp_total', ...
    'p_iso_emp','p_iso_emp_iteration_mean', ...
    'p_iso_emp_std','p_iso_emp_sem', ...
    'n_valid_iterations', ...
    'isolated_total_th', ...
    'isolated_total_emp_iterations', ...
    'isolated_total_emp', ...
    'isolated_total_emp_std', ...
    'isolated_total_emp_sem', ...
    'difference','valid_compare', ...
    'rmse_grid','mae_grid','bias_grid', ...
    'rmse_by_phi','mae_by_phi','bias_by_phi', ...
    '-v7.3');

fprintf('Resultats sauvegardes dans %s\n',output_file);

%% ============================================================
%  Fonction locale
%% ============================================================

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
