%% vitesse_relative_t_phi_compare.m
% Comparaison theorie / empirique de la vitesse relative locale
% conditionnee a l'existence d'un lien :
%
%   v_rel^link(t,phi)
%     = E[ ||v_1-v_2|| | d_12 <= dmax, Phi_1 = phi ].
%
% La latitude phi est celle de la PREMIERE extremite du lien.
% Un lien non oriente est donc considere deux fois :
% une fois depuis chacune de ses extremites. Cette convention est
% exactement coherente avec la definition locale de p_link(t,phi).
%
% Theorie :
%
%   v_rel^th(t,phi)
%      = E[ ||v_1-v_2|| 1_{lien} | Phi_1=phi ]
%        ------------------------------------------------
%                    p_link(t,phi).
%
% Empirique, dans une tranche b et une simulation r :
%
%   v_rel,b^{emp,(r)}(t)
%      = somme des vitesses relatives des liens orientes
%        dont la premiere extremite appartient a b
%        -------------------------------------------------
%        nombre de liens orientes issus de la tranche b.
%
% Le script :
%   1) charge uniquement les parametres et grilles depuis
%      plink_t_phi_results.mat ;
%   2) resimule plusieurs realisations independantes ;
%   3) calcule la theorie par quadrature sur u_2 et DeltaOmega ;
%   4) compare theorie et moyenne empirique avec SEM.
%
% Entree :
%   plink_t_phi_results.mat
%
% Sortie :
%   vitesse_relative_t_phi_results.mat

clear; clc; close all;

%% ============================================================
%  1. Chargement des parametres
%% ============================================================

script_dir = fileparts(mfilename('fullpath'));

input_candidates = {
    fullfile(script_dir,'plink_t_phi_results.mat')
    fullfile(script_dir,'..','plink_t_phi_results.mat')
};

input_file = '';

for k = 1:numel(input_candidates)
    if isfile(input_candidates{k})
        input_file = input_candidates{k};
        break;
    end
end

if isempty(input_file)
    error('Fichier plink_t_phi_results.mat introuvable.');
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

v_orb = R*omega;

time_values_full = double(S.time_theory(:).');

phi_edges_emp = double(S.phi_edges_emp(:).');
phi_vals_emp = double(S.phi_vals_emp(:).');
dphi_emp = double(S.dphi_emp(:).');

Nb = numel(phi_vals_emp);

%% ============================================================
%  2. Parametres propres a ce calcul
%% ============================================================

% Nombre de simulations empiriques.
n_iterations = 50;

% Pour accelerer le calcul, on peut ne garder qu'un instant sur q.
% Mettre 1 pour conserver tous les instants de plink_t_phi_results.mat.
time_stride = 1;

time_indices = 1:time_stride:numel(time_values_full);
time_values = time_values_full(time_indices);
Nt = numel(time_values);

% Quadratures theoriques.
n_phi_bins_th = 150;
n_u_quad = 360;
n_delta_omega_quad = 360;

% Seuil minimal de liens orientes cumules dans une tranche.
min_oriented_links_total = 20;

% Nombre de coupes affichees.
n_selected_times = 5;

% Graine reproductible.
rng_seed = 2;
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

delta_omega_grid = ...
    (0:n_delta_omega_quad-1)*(2*pi/n_delta_omega_quad);

alpha_max = 2*asin(min(dmax/(2*R),1));
cos_alpha_max = cos(alpha_max);

s_th = sin(phi_vals_th)/sin(inc);
s_th = min(max(s_th,-1),1);

u_plus_th = mod(asin(s_th),2*pi);
u_minus_th = mod(pi-asin(s_th),2*pi);

%% ============================================================
%  4. Noyaux geometriques
%
% G(u1,u2) :
%   probabilite de lien apres moyenne sur DeltaOmega.
%
% H(u1,u2) :
%   E[ v_rel * 1_{lien} | u1,u2 ] apres moyenne sur DeltaOmega.
%% ============================================================

G_plus = zeros(n_phi_bins_th,n_u_quad);
G_minus = zeros(n_phi_bins_th,n_u_quad);

H_plus = zeros(n_phi_bins_th,n_u_quad);
H_minus = zeros(n_phi_bins_th,n_u_quad);

fprintf('Calcul des noyaux theoriques...\n');

for b = 1:n_phi_bins_th

    [G_plus(b,:),H_plus(b,:)] = ...
        link_velocity_kernels( ...
            u_plus_th(b),u2_grid,delta_omega_grid, ...
            inc,cos_alpha_max,v_orb);

    [G_minus(b,:),H_minus(b,:)] = ...
        link_velocity_kernels( ...
            u_minus_th(b),u2_grid,delta_omega_grid, ...
            inc,cos_alpha_max,v_orb);

    if mod(b,10) == 0 || b == n_phi_bins_th
        fprintf('  latitude theorique %d/%d\n', ...
            b,n_phi_bins_th);
    end
end

%% ============================================================
%  5. Theorie fine v_rel^link(t,phi)
%% ============================================================

p_link_th_fine = nan(Nt,n_phi_bins_th);
v_rel_th_fine = nan(Nt,n_phi_bins_th);
f_phi_th_fine = nan(Nt,n_phi_bins_th);

for t_idx = 1:Nt

    t = time_values(t_idx);

    % Loi de phase du second satellite.
    f_u2 = abs(cos(u2_grid-omega*t))/4;
    f_u2 = f_u2/(sum(f_u2)*du);

    % Poids conditionnels des deux branches du premier satellite.
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

    % Denominateur : probabilite de lien.
    p_plus = (G_plus*f_u2(:))*du;
    p_minus = (G_minus*f_u2(:))*du;

    p_local = ...
        w_plus.*p_plus.' ...
        + w_minus.*p_minus.';

    % Numerateur : E[v_rel 1_lien].
    h_plus = (H_plus*f_u2(:))*du;
    h_minus = (H_minus*f_u2(:))*du;

    h_local = ...
        w_plus.*h_plus.' ...
        + w_minus.*h_minus.';

    p_link_th_fine(t_idx,:) = p_local;

    valid_link = p_local > 1e-14;
    v_rel_th_fine(t_idx,valid_link) = ...
        h_local(valid_link)./p_local(valid_link);

    % Loi de latitude theorique.
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
%  6. Theorie moyennee dans les tranches empiriques
%
% La moyenne doit etre ponderee par le nombre attendu de liens
% orientes, donc par :
%
%   f_Phi(t,phi) p_link(t,phi).
%% ============================================================

v_rel_th_on_emp = nan(Nt,Nb);
p_link_th_on_emp = nan(Nt,Nb);
oriented_link_mass_th = zeros(Nt,Nb);

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

        link_weight = ...
            satellite_mass ...
            .* p_link_th_fine(t_idx,in_bin);

        mass_link = sum(link_weight);

        oriented_link_mass_th(t_idx,b) = mass_link;

        if mass_link > 0
            v_rel_th_on_emp(t_idx,b) = ...
                sum( ...
                    v_rel_th_fine(t_idx,in_bin) ...
                    .* link_weight) ...
                / mass_link;
        end

        mass_sat = sum(satellite_mass);

        if mass_sat > 0
            p_link_th_on_emp(t_idx,b) = ...
                sum( ...
                    p_link_th_fine(t_idx,in_bin) ...
                    .* satellite_mass) ...
                / mass_sat;
        end
    end
end

%% ============================================================
%  7. Resimulation empirique
%% ============================================================

oriented_link_count_iterations = ...
    zeros(n_iterations,Nt,Nb);

relative_speed_sum_iterations = ...
    zeros(n_iterations,Nt,Nb);

v_rel_emp_iterations = ...
    nan(n_iterations,Nt,Nb);

for r = 1:n_iterations

    % Uniformite spatiale initiale dans la bande.
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

        % Positions.
        X = R * [ ...
            cos_Omega.*cos_u ...
                - sin_Omega.*sin_u*cos(inc), ...
            sin_Omega.*cos_u ...
                + cos_Omega.*sin_u*cos(inc), ...
            sin_u*sin(inc)];

        % Vitesses.
        V = v_orb * [ ...
            -cos_Omega.*sin_u ...
                - sin_Omega.*cos_u*cos(inc), ...
            -sin_Omega.*sin_u ...
                + cos_Omega.*cos_u*cos(inc), ...
            cos_u*sin(inc)];

        latitude = asin(max(min(X(:,3)/R,1),-1));

        % Liens.
        gram_x = X*X.';
        norm_x2 = sum(X.^2,2);

        D2 = ...
            norm_x2+norm_x2.'-2*gram_x;

        D2 = max(D2,0);

        adjacency = D2 <= dmax^2;
        adjacency(1:N+1:end) = false;
        adjacency = adjacency | adjacency.';

        % Vitesses relatives de toutes les paires.
        gram_v = V*V.';
        norm_v2 = sum(V.^2,2);

        Vrel2 = ...
            norm_v2+norm_v2.'-2*gram_v;

        Vrel = sqrt(max(Vrel2,0));

        % Pour chaque satellite i :
        %   nombre de liens sortants = deg_i
        %   somme des vitesses relatives de ses liens.
        degree = sum(adjacency,2);
        speed_sum_by_satellite = sum(Vrel.*adjacency,2);

        bin_id = discretize(latitude,phi_edges_emp);

        for b = 1:Nb

            mask = bin_id == b;

            n_oriented_links = sum(degree(mask));
            speed_sum = sum(speed_sum_by_satellite(mask));

            oriented_link_count_iterations(r,t_idx,b) = ...
                n_oriented_links;

            relative_speed_sum_iterations(r,t_idx,b) = ...
                speed_sum;

            if n_oriented_links > 0
                v_rel_emp_iterations(r,t_idx,b) = ...
                    speed_sum/n_oriented_links;
            end
        end
    end

    fprintf('Realisation empirique %d/%d terminee.\n', ...
        r,n_iterations);
end

%% ============================================================
%  8. Moyenne empirique
%
% Estimateur principal : regroupement de tous les liens de toutes
% les simulations dans chaque cellule (t,b).
%% ============================================================

oriented_link_count_total = squeeze(sum( ...
    oriented_link_count_iterations,1));

relative_speed_sum_total = squeeze(sum( ...
    relative_speed_sum_iterations,1));

v_rel_emp = nan(Nt,Nb);

valid_total = ...
    oriented_link_count_total >= min_oriented_links_total;

v_rel_emp(valid_total) = ...
    relative_speed_sum_total(valid_total) ...
    ./ oriented_link_count_total(valid_total);

% Dispersion entre simulations.
v_rel_emp_iteration_mean = squeeze(mean( ...
    v_rel_emp_iterations,1,'omitnan'));

v_rel_emp_std = squeeze(std( ...
    v_rel_emp_iterations,0,1,'omitnan'));

n_valid_iterations = squeeze(sum( ...
    isfinite(v_rel_emp_iterations),1));

v_rel_emp_sem = ...
    v_rel_emp_std ...
    ./ sqrt(max(n_valid_iterations,1));

%% ============================================================
%  9. Vitesse relative globale conditionnee au lien
%% ============================================================

% Theorie globale.
global_link_weight = ...
    f_phi_th_fine ...
    .* p_link_th_fine ...
    .* dphi_th;

v_rel_th_global = ...
    sum(v_rel_th_fine.*global_link_weight,2,'omitnan') ...
    ./ sum(global_link_weight,2,'omitnan');

% Empirique globale par simulation.
v_rel_emp_global_iterations = nan(n_iterations,Nt);

for r = 1:n_iterations
    count_r = sum( ...
        oriented_link_count_iterations(r,:,:),3);

    speed_r = sum( ...
        relative_speed_sum_iterations(r,:,:),3);

    count_r = reshape(count_r,Nt,1);
    speed_r = reshape(speed_r,Nt,1);

    valid_r = count_r > 0;

    temp = nan(Nt,1);
    temp(valid_r) = speed_r(valid_r)./count_r(valid_r);

    v_rel_emp_global_iterations(r,:) = temp.';
end

v_rel_emp_global = mean( ...
    v_rel_emp_global_iterations,1,'omitnan').';

v_rel_emp_global_std = std( ...
    v_rel_emp_global_iterations,0,1,'omitnan').';

v_rel_emp_global_sem = ...
    v_rel_emp_global_std/sqrt(n_iterations);

%% ============================================================
%  10. Diagnostics
%% ============================================================

valid_compare = ...
    isfinite(v_rel_th_on_emp) ...
    & isfinite(v_rel_emp);

difference = nan(Nt,Nb);
difference(valid_compare) = ...
    v_rel_emp(valid_compare) ...
    - v_rel_th_on_emp(valid_compare);

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
    v_rel_th_on_emp.');
axis xy;
colorbar;
xlabel('Temps (s)');
ylabel('Latitude \phi (deg)');
title('v_{rel}^{th}(t,\phi) conditionnee au lien');

figure;
imagesc( ...
    time_values, ...
    rad2deg(phi_vals_emp), ...
    v_rel_emp.');
axis xy;
colorbar;
xlabel('Temps (s)');
ylabel('Latitude \phi (deg)');
title(sprintf( ...
    'v_{rel}^{emp}(t,\\phi) moyenne sur %d simulations', ...
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
title('Ecart v_{rel}^{emp}-v_{rel}^{th}');

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

    plot(rad2deg(phi_vals_th), ...
        v_rel_th_fine(it,:), ...
        'LineWidth',2, ...
        'DisplayName','Theorie fine');

    plot(rad2deg(phi_vals_emp), ...
        v_rel_th_on_emp(it,:), ...
        's', ...
        'LineWidth',1.1, ...
        'DisplayName','Theorie moyenne par tranche');

    errorbar(rad2deg(phi_vals_emp), ...
        v_rel_emp(it,:), ...
        v_rel_emp_sem(it,:), ...
        'o-', ...
        'LineWidth',1.2, ...
        'DisplayName','Empirique moyen \pm SEM');

    grid on;
    ylabel('v_{rel} (km/s)');
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

plot(time_values,v_rel_th_global, ...
    'LineWidth',2, ...
    'DisplayName','Theorie');

plot(time_values,v_rel_emp_global, ...
    'LineWidth',1.5, ...
    'DisplayName','Empirique moyen');

plot(time_values, ...
    v_rel_emp_global+v_rel_emp_global_sem, ...
    ':','LineWidth',1, ...
    'DisplayName','Empirique + SEM');

plot(time_values, ...
    v_rel_emp_global-v_rel_emp_global_sem, ...
    ':','LineWidth',1, ...
    'DisplayName','Empirique - SEM');

grid on;
xlabel('Temps (s)');
ylabel('v_{rel}^{link}(t) (km/s)');
title('Vitesse relative globale conditionnee au lien');
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
ylabel('Erreur de vitesse (km/s)');
title('Erreur de v_{rel}^{link}(t,\phi)');
legend('Location','best');
hold off;

%% ============================================================
%  12. Affichage console
%% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' VITESSE RELATIVE v_rel(t,phi) - THEORIE / EMPIRIQUE\n');
fprintf('============================================================\n');
fprintf('Fichier charge                    : %s\n',input_file);
fprintf('N                                 : %d\n',round(N));
fprintf('R                                 : %.8f km\n',R);
fprintf('v_orb                             : %.8f km/s\n',v_orb);
fprintf('Nombre de simulations             : %d\n',n_iterations);
fprintf('Nombre d''instants                : %d\n',Nt);
fprintf('Nombre de tranches empiriques     : %d\n',Nb);
fprintf('Quadrature u2                     : %d points\n',n_u_quad);
fprintf('Quadrature DeltaOmega             : %d points\n', ...
    n_delta_omega_quad);
fprintf('------------------------------------------------------------\n');
fprintf('RMSE locale                       : %.10e km/s\n', ...
    rmse_grid);
fprintf('MAE locale                        : %.10e km/s\n', ...
    mae_grid);
fprintf('Biais moyen emp-theorie           : %.10e km/s\n', ...
    bias_grid);
fprintf('------------------------------------------------------------\n');
fprintf('Moyenne globale empirique         : %.10f km/s\n', ...
    mean(v_rel_emp_global,'omitnan'));
fprintf('Moyenne globale theorique         : %.10f km/s\n', ...
    mean(v_rel_th_global,'omitnan'));
fprintf('============================================================\n');

%% ============================================================
%  13. Sauvegarde
%% ============================================================

output_file = fullfile( ...
    script_dir, ...
    'vrel_t_phi_results.mat');

save(output_file, ...
    'input_file', ...
    'N','R','dmax','inc','omega','v_orb', ...
    'n_iterations','rng_seed', ...
    'time_stride','time_indices','time_values','Nt', ...
    'phi_edges_emp','phi_vals_emp','dphi_emp','Nb', ...
    'phi_edges_th','phi_vals_th','dphi_th', ...
    'n_phi_bins_th','n_u_quad','n_delta_omega_quad', ...
    'u2_grid','delta_omega_grid','du', ...
    'alpha_max','cos_alpha_max', ...
    'u_plus_th','u_minus_th', ...
    'G_plus','G_minus','H_plus','H_minus', ...
    'f_phi_th_fine', ...
    'p_link_th_fine','p_link_th_on_emp', ...
    'v_rel_th_fine','v_rel_th_on_emp', ...
    'oriented_link_mass_th', ...
    'oriented_link_count_iterations', ...
    'relative_speed_sum_iterations', ...
    'v_rel_emp_iterations', ...
    'oriented_link_count_total', ...
    'relative_speed_sum_total', ...
    'v_rel_emp','v_rel_emp_iteration_mean', ...
    'v_rel_emp_std','v_rel_emp_sem', ...
    'n_valid_iterations', ...
    'v_rel_th_global', ...
    'v_rel_emp_global_iterations', ...
    'v_rel_emp_global', ...
    'v_rel_emp_global_std', ...
    'v_rel_emp_global_sem', ...
    'difference','valid_compare', ...
    'rmse_grid','mae_grid','bias_grid', ...
    'rmse_by_phi','mae_by_phi','bias_by_phi', ...
    '-v7.3');

fprintf('Resultats sauvegardes dans %s\n',output_file);

%% ============================================================
%  Fonction locale : noyaux de lien et de vitesse
%% ============================================================

function [G,H] = link_velocity_kernels( ...
    u1,u2_grid,delta_omega_grid,inc,cos_alpha_max,v_orb)

% Pour chaque u2 :
%
%   G(u1,u2) = moyenne_DeltaOmega 1_{lien}
%
%   H(u1,u2) = moyenne_DeltaOmega
%              [v_rel 1_{lien}].

    u2 = u2_grid(:);
    dOmega = delta_omega_grid(:).';

    cos_u1 = cos(u1);
    sin_u1 = sin(u1);

    cos_u2 = cos(u2);
    sin_u2 = sin(u2);

    cos_dO = cos(dOmega);
    sin_dO = sin(dOmega);

    % Produit scalaire normalise des positions.
    A_pos = ...
        cos_u1.*cos_u2 ...
        + cos(inc)^2.*sin_u1.*sin_u2;

    B_pos = ...
        cos(inc).*( ...
        sin_u1.*cos_u2 ...
        - cos_u1.*sin_u2);

    C_pos = ...
        sin(inc)^2.*sin_u1.*sin_u2;

    dot_position = ...
        A_pos.*cos_dO ...
        + B_pos.*sin_dO ...
        + C_pos;

    is_link = dot_position >= cos_alpha_max;

    % Vitesse normalisee du premier satellite pour Omega_1=0.
    v1x = -sin_u1;
    v1y = cos_u1*cos(inc);
    v1z = cos_u1*sin(inc);

    % Vitesse normalisee du second satellite.
    v2x = ...
        -cos_dO.*sin_u2 ...
        - sin_dO.*cos_u2*cos(inc);

    v2y = ...
        -sin_dO.*sin_u2 ...
        + cos_dO.*cos_u2*cos(inc);

    v2z = ...
        cos_u2*sin(inc);

    dot_velocity = ...
        v1x.*v2x ...
        + v1y.*v2y ...
        + v1z.*v2z;

    dot_velocity = min(max(dot_velocity,-1),1);

    relative_speed = ...
        v_orb*sqrt(max(2-2*dot_velocity,0));

    G = mean(is_link,2).';
    H = mean(relative_speed.*is_link,2).';
end
