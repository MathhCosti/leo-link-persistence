%% beta0_t_phi_th_delta_spatial.m
% beta_0 theorique local en fonction du temps et de la latitude
% Walker Delta a uniformite spatiale.
%
% Entree :
%   lambda_t_phi_results.mat
%   produit par lambda_t_phi.m
%
% Theorie locale :
%
%   c_s(t,phi)
%     = lambda(t,phi)^s / s!
%       * s^(s-2) * A_cap^(s-1)
%       * E_q[ exp(-lambda(t,phi) A_union(X))/tau_h(G(X)) ]
%
% puis :
%
%   beta0_density_t_phi(t,phi)
%     = sum_s c_s(t,phi)
%
% en composantes/km^2.
%
% Nombre moyen de composantes dans chaque tranche :
%
%   beta0_bin_th(t,b)
%     = beta0_density_t_phi(t,b) * A_b
%
% et nombre total de composantes :
%
%   beta0_total_th(t)
%     = sum_b beta0_bin_th(t,b).
%
% IMPORTANT :
% - Approximation locale : lambda(t,phi') ~ lambda(t,phi) a l'echelle
%   d'une composante.
% - Le Monte-Carlo sert seulement a evaluer l'integrale theorique.
% - Une meme banque de configurations geometriques est reutilisee
%   pour tous les couples (t,phi).

clear; clc; close all;
rng(1);

%% ============================================================
%  1. CHARGEMENT DE lambda(t,phi)
%% ============================================================

script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end

lambda_candidates = {
    fullfile(script_dir,'lambda_t_phi_results.mat')
    fullfile(script_dir,'..','lambda_t_phi_results.mat')
};

lambda_file = '';
for k = 1:numel(lambda_candidates)
    if isfile(lambda_candidates{k})
        lambda_file = lambda_candidates{k};
        break;
    end
end

if isempty(lambda_file)
    error(['Fichier lambda_t_phi_results.mat introuvable. ', ...
           'Executer lambda_t_phi.m auparavant.']);
end

S = load(lambda_file);

required_fields = { ...
    'N','R','time_values', ...
    'phi_edges_emp','phi_vals_emp','dphi_emp', ...
    'area_bin','lambda_bin_th'};

for k = 1:numel(required_fields)
    if ~isfield(S,required_fields{k})
        error('Le fichier %s doit contenir %s.', ...
            lambda_file,required_fields{k});
    end
end

N = double(S.N);
R = double(S.R);

time_values = double(S.time_values(:));
phi_edges = double(S.phi_edges_emp(:).');
phi_vals = double(S.phi_vals_emp(:).');
dphi_bins = double(S.dphi_emp(:).');
A_phi_bins = double(S.area_bin(:).');

lambda_t_phi = double(S.lambda_bin_th);

[Nt,Nb] = size(lambda_t_phi);

if numel(time_values) ~= Nt
    error('Dimension incompatible pour time_values.');
end

if numel(phi_vals) ~= Nb || ...
        numel(dphi_bins) ~= Nb || ...
        numel(A_phi_bins) ~= Nb
    error('Dimensions incompatibles pour la grille de latitude.');
end

%% ============================================================
%  2. PARAMETRES PHYSIQUES ET NUMERIQUES
%% ============================================================

% Recuperation de dmax, inc et omega.
% Pour reproduire exactement l'approche empirique de betti_t_phi.m,
% on cherche en priorite N1_t_phi_results.mat, qui contient ces parametres.
dmax = 1500;
inc = NaN;
omega = NaN;
mu = 398600; % km^3/s^2

param_candidates = {
    fullfile(script_dir,'N1_t_phi_results.mat')
    fullfile(script_dir,'..','N1_t_phi_results.mat')
};

param_file = '';
for k = 1:numel(param_candidates)
    if isfile(param_candidates{k})
        param_file = param_candidates{k};
        SP = load(param_file);

        if isfield(SP,'dmax'), dmax = double(SP.dmax); end
        if isfield(SP,'inc'), inc = double(SP.inc); end
        if isfield(SP,'omega'), omega = double(SP.omega); end
        break;
    end
end

% Repli sur le fichier source de lambda(t,phi) si disponible.
if (~isfinite(inc) || ~isfinite(omega)) && isfield(S,'input_file')
    plink_file = char(S.input_file);
    if isfile(plink_file)
        SPL = load(plink_file);
        if isfield(SPL,'dmax'), dmax = double(SPL.dmax); end
        if isfield(SPL,'inc'), inc = double(SPL.inc); end
        if isfield(SPL,'omega'), omega = double(SPL.omega); end
    end
end

% Dernier repli pour omega si l'inclinaison est connue.
if ~isfinite(omega)
    omega = sqrt(mu/R^3);
end

if ~isfinite(inc)
    error(['Inclinaison inc introuvable. Placer N1_t_phi_results.mat ', ...
           'dans le dossier du script ou verifier le fichier source ', ...
           'plink_t_phi_results.mat.']);
end

Smax = 25;

Nprobe = 4000;
Nsamp = 12000;

rho_tilt = 0.75;
beta_shape = 3.0;

% Comparaison empirique COMPLETE comme dans betti_t_phi.m :
% plusieurs constellations independantes, toutes les composantes comptees.
compare_empirical = true;
n_iterations_beta0_emp = 100;
rng_seed_beta0 = 31;

fprintf('============================================================\n');
fprintf('beta0(t,phi) theorique - Delta uniformite spatiale\n');
fprintf('N                  = %.3f\n',N);
fprintf('R                  = %.3f km\n',R);
fprintf('dmax               = %.3f km\n',dmax);
fprintf('Nt / Nb            = %d / %d\n',Nt,Nb);
fprintf('Smax               = %d\n',Smax);
fprintf('Nsamp / taille     = %d\n',Nsamp);
fprintf('Nprobe             = %d\n',Nprobe);
fprintf('Comparaison emp    = %d\n',compare_empirical);
fprintf('n simulations emp  = %d\n',n_iterations_beta0_emp);
fprintf('inclinaison        = %.3f deg\n',rad2deg(inc));
fprintf('omega              = %.8e rad/s\n',omega);
fprintf('============================================================\n\n');

%% ============================================================
%  3. GEOMETRIE DE CONNEXION
%% ============================================================

A_sphere = 4*pi*R^2;

alpha_max = 2*asin(min(1,dmax/(2*R)));
cos_alpha = cos(alpha_max);
one_minus_cos_alpha = 1-cos_alpha;

A_cap = 2*pi*R^2*one_minus_cos_alpha;

%% ============================================================
%  4. POINTS QMC POUR A_union
%% ============================================================

probe_unit = fibonacci_sphere(Nprobe);

%% ============================================================
%  5. TABLEAUX
%% ============================================================

% c_s(t,phi), stocke taille par taille pour diagnostics.
EC_density_t_phi = zeros(Smax,Nt,Nb);

% beta0 local surfacique [composantes/km^2]
beta0_density_t_phi = zeros(Nt,Nb);

% ESS local
ESS_t_phi = nan(Smax,Nt,Nb);

%% ============================================================
%  6. s = 1 EXACT
%% ============================================================

C1 = lambda_t_phi .* exp(-lambda_t_phi*A_cap);

EC_density_t_phi(1,:,:) = reshape(C1,1,Nt,Nb);
beta0_density_t_phi = beta0_density_t_phi + C1;
ESS_t_phi(1,:,:) = Inf;

fprintf('s =  1 : exact sur toute la grille (t,phi).\n');

%% ============================================================
%  7. s >= 2 : IMPORTANCE SAMPLING
%
% Pour chaque s :
% - on genere UNE banque de configurations geometriques ;
% - on calcule A_union et tau_h une seule fois ;
% - on repondere ensuite pour toutes les valeurs lambda(t,phi).
%% ============================================================

x1 = [0 0 1];

lambda_flat = lambda_t_phi(:);
n_grid = numel(lambda_flat);

for s = 2:Smax

    Aunion_all = zeros(Nsamp,1);
    logtau_all = zeros(Nsamp,1);

    for mm = 1:Nsamp

        % 1) Arbre de Cayley uniforme
        edges = random_labeled_tree_prufer(s);

        % 2) Orientation depuis la racine 1
        [parent,order] = root_tree(edges,s,1);

        % 3) Generation geometrique selon h(u)
        X = zeros(s,3);
        X(1,:) = x1;

        for kk = 2:s
            v = order(kk);
            p = parent(v);

            X(v,:) = sample_cap_around_tilted( ...
                X(p,:),alpha_max,rho_tilt,beta_shape);
        end

        % 4) Graphe geometrique induit
        Ddot = max(-1,min(1,X*X.'));
        Adj = (Ddot >= cos_alpha-1e-12);
        Adj(1:s+1:end) = false;
        Adj = Adj | Adj.';

        if ~is_connected_adj(Adj)
            error('Le sampler a produit un graphe non connexe.');
        end

        % 5) Poids h_ij
        U = (1-Ddot)/one_minus_cos_alpha;
        U = max(0,min(1,U));

        H = zeros(s,s);
        mask = Adj;

        H(mask) = ...
            (1-rho_tilt) ...
            + rho_tilt*beta_shape .* U(mask).^(beta_shape-1);

        H = (H+H.')/2;
        H(1:s+1:end) = 0;

        % 6) tau_h(G)
        logtau_all(mm) = ...
            log_weighted_spanning_tree_sum(H);

        % 7) Aire de l'union des calottes
        maxdot = max(probe_unit*X.',[],2);
        frac_union = mean(maxdot >= cos_alpha);

        Aunion_all(mm) = A_sphere*frac_union;
    end

    % Partie independante de lambda(t,phi)
    log_geom_prefactor = ...
        -gammaln(s+1) ...
        + (s-2)*log(s) ...
        + (s-1)*log(A_cap);

    Cs_flat = zeros(n_grid,1);
    ESS_flat = zeros(n_grid,1);

    for g = 1:n_grid

        lam = lambda_flat(g);

        if ~isfinite(lam) || lam <= 0
            Cs_flat(g) = 0;
            ESS_flat(g) = NaN;
            continue;
        end

        logw = -lam*Aunion_all - logtau_all;

        m = max(logw);
        log_mean_weight = ...
            m + log(mean(exp(logw-m)));

        log_cs = ...
            s*log(lam) ...
            + log_geom_prefactor ...
            + log_mean_weight;

        Cs_flat(g) = exp(log_cs);

        % ESS
        lsw = logsumexp_vec(logw);
        lsw2 = logsumexp_vec(2*logw);

        ESS_flat(g) = ...
            exp(2*lsw-lsw2);
    end

    Cs = reshape(Cs_flat,Nt,Nb);
    ESSs = reshape(ESS_flat,Nt,Nb);

    EC_density_t_phi(s,:,:) = reshape(Cs,1,Nt,Nb);
    ESS_t_phi(s,:,:) = reshape(ESSs,1,Nt,Nb);

    beta0_density_t_phi = ...
        beta0_density_t_phi + Cs;

    fprintf(['s = %2d : termine | ESS min/median = %.1f / %.1f ' ...
             '| max c_s = %.3e comp/km^2\n'], ...
        s, ...
        min(ESS_flat,[],'omitnan'), ...
        median(ESS_flat,'omitnan'), ...
        max(Cs_flat,[],'omitnan'));
end

%% ============================================================
%  8. NOMBRE DE COMPOSANTES PAR TRANCHE ET TOTAL
%% ============================================================

beta0_bin_th = ...
    beta0_density_t_phi .* A_phi_bins;

beta0_total_th_t = ...
    sum(beta0_bin_th,2,'omitnan');

% Densite par radian de latitude
beta0_per_rad_t_phi = ...
    beta0_bin_th ./ dphi_bins;

% Moyennes utiles
beta0_density_phi_mean_t = ...
    mean(beta0_density_t_phi,1,'omitnan');

beta0_per_rad_phi_mean_t = ...
    mean(beta0_per_rad_t_phi,1,'omitnan');

beta0_total_th_mean = ...
    mean(beta0_total_th_t,'omitnan');

%% ============================================================
%  9. COMPARAISON EMPIRIQUE COMPLETE
%
% Meme methode que dans betti_t_phi.m :
%
% - on genere n_iterations_beta0_emp constellations independantes ;
% - les positions initiales sont uniformes spatialement dans la bande ;
% - chaque satellite evolue ensuite sur son orbite ;
% - beta0 empirique est calcule directement sur le graphe complet ;
% - une composante de taille s contribue pour 1/s dans la tranche
%   de latitude de chacun de ses membres.
%
% Ainsi, pour chaque realisation et chaque instant :
%
%   sum_b beta0_emp(r,t,b) = beta0_emp,total(r,t).
%% ============================================================

beta0_emp_iterations = zeros( ...
    n_iterations_beta0_emp,Nt,Nb);

beta0_total_emp_iterations = zeros( ...
    n_iterations_beta0_emp,Nt);

if compare_empirical

    rng(rng_seed_beta0);

    for r = 1:n_iterations_beta0_emp

        % Uniformite SPATIALE initiale dans la bande |phi| <= inc.
        [u0,Omega] = sample_initial_orbits_spatial(N,inc);

        for it = 1:Nt

            X = positions_from_orbits( ...
                u0,Omega,time_values(it),R,inc,omega);

            latitude = asin(max(min(X(:,3)/R,1),-1));
            bin_id = discretize(latitude,phi_edges);

            A = adjacency_from_positions(X,dmax);
            G = graph(sparse(A));

            component_id = conncomp(G);
            component_sizes = accumarray(component_id(:),1);

            beta0_total_emp_iterations(r,it) = ...
                numel(component_sizes);

            % Attribution locale symetrique.
            for c = 1:numel(component_sizes)

                members = find(component_id==c);
                sc = component_sizes(c);

                for m = members(:).'
                    b = bin_id(m);

                    if ~isnan(b)
                        beta0_emp_iterations(r,it,b) = ...
                            beta0_emp_iterations(r,it,b)+1/sc;
                    end
                end
            end
        end

        fprintf('Beta0 empirique complet : realisation %d/%d\n', ...
            r,n_iterations_beta0_emp);
    end
end

% Moyenne empirique locale par tranche.
beta0_bin_emp = squeeze(mean( ...
    beta0_emp_iterations,1,'omitnan'));

beta0_bin_emp_std = squeeze(std( ...
    beta0_emp_iterations,0,1,'omitnan'));

beta0_bin_emp_sem = ...
    beta0_bin_emp_std/sqrt(n_iterations_beta0_emp);

% Valeur empirique globale.
beta0_total_emp_t = mean( ...
    beta0_total_emp_iterations,1,'omitnan').';

beta0_total_emp_std = std( ...
    beta0_total_emp_iterations,0,1,'omitnan').';

beta0_total_emp_sem = ...
    beta0_total_emp_std/sqrt(n_iterations_beta0_emp);

% Conversion en composantes/radian pour les cartes.
beta0_per_rad_emp = ...
    beta0_bin_emp ./ dphi_bins;

beta0_per_rad_emp_sem = ...
    beta0_bin_emp_sem ./ dphi_bins;

% Verification exacte de la repartition locale.
beta0_total_emp_from_local_iterations = squeeze(sum( ...
    beta0_emp_iterations,3,'omitnan'));

max_local_global_error_emp = max(abs( ...
    beta0_total_emp_from_local_iterations(:) ...
    - beta0_total_emp_iterations(:)));

% Diagnostics theorie / empirique sur le NOMBRE PAR TRANCHE,
% comme dans betti_t_phi.m.
difference_bin = beta0_bin_emp-beta0_bin_th;

valid = isfinite(difference_bin);

rmse_t_phi = sqrt(mean(difference_bin(valid).^2));
mae_t_phi = mean(abs(difference_bin(valid)));
bias_t_phi = mean(difference_bin(valid));

% Diagnostics par latitude.
rmse_by_phi = nan(1,Nb);
mae_by_phi = nan(1,Nb);
bias_by_phi = nan(1,Nb);

for b = 1:Nb
    valid_b = isfinite(difference_bin(:,b));

    if any(valid_b)
        err_b = difference_bin(valid_b,b);
        rmse_by_phi(b) = sqrt(mean(err_b.^2));
        mae_by_phi(b) = mean(abs(err_b));
        bias_by_phi(b) = mean(err_b);
    end
end

% Diagnostics globaux beta0(t).
difference_total = beta0_total_emp_t-beta0_total_th_t;

rmse_total_t = sqrt(mean(difference_total.^2,'omitnan'));
mae_total_t = mean(abs(difference_total),'omitnan');
bias_total_t = mean(difference_total,'omitnan');

%% ============================================================
%  10. AFFICHAGE CONSOLE
%% ============================================================

fprintf('\n============================================================\n');
fprintf('RESULTATS beta0(t,phi)\n');
fprintf('============================================================\n');
fprintf('beta0 theorique moyen temporel      : %.8f\n', ...
    beta0_total_th_mean);
fprintf('beta0 empirique moyen temporel      : %.8f\n', ...
    mean(beta0_total_emp_t,'omitnan'));
fprintf('Nombre de simulations empiriques    : %d\n', ...
    n_iterations_beta0_emp);
fprintf('RMSE locale par tranche             : %.6e\n', ...
    rmse_t_phi);
fprintf('MAE locale par tranche              : %.6e\n', ...
    mae_t_phi);
fprintf('Biais local emp-theorie             : %.6e\n', ...
    bias_t_phi);
fprintf('RMSE beta0 total(t)                  : %.6e\n', ...
    rmse_total_t);
fprintf('Erreur max somme locale/globale emp : %.6e\n', ...
    max_local_global_error_emp);

fprintf('============================================================\n');

%% ============================================================
%  11. FIGURES
%% ============================================================

% 1) Heatmap theorique : NOMBRE DE COMPOSANTES PAR TRANCHE
figure;
imagesc( ...
    time_values, ...
    rad2deg(phi_vals), ...
    beta0_bin_th.');
axis xy;
colorbar;
xlabel('Temps (s)');
ylabel('Latitude \phi (deg)');
title('\beta_0^{th}(t,\phi) par tranche - nouvelle theorie PPP');

% 2) Heatmap empirique complete
figure;
imagesc( ...
    time_values, ...
    rad2deg(phi_vals), ...
    beta0_bin_emp.');
axis xy;
colorbar;
xlabel('Temps (s)');
ylabel('Latitude \phi (deg)');
title(sprintf( ...
    '\\beta_0^{emp}(t,\\phi) - moyenne sur %d simulations', ...
    n_iterations_beta0_emp));

% 3) Heatmap de l'ecart emp - theorie
figure;
imagesc( ...
    time_values, ...
    rad2deg(phi_vals), ...
    difference_bin.');
axis xy;
colorbar;
xlabel('Temps (s)');
ylabel('Latitude \phi (deg)');
title('\beta_0^{emp}(t,\phi)-\beta_0^{th}(t,\phi) par tranche');

% 4) Coupes a plusieurs instants, avec SEM empirique.
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

    plot(rad2deg(phi_vals), ...
        beta0_bin_th(it,:), ...
        's-','LineWidth',1.8, ...
        'DisplayName','Nouvelle theorie PPP');

    errorbar(rad2deg(phi_vals), ...
        beta0_bin_emp(it,:), ...
        beta0_bin_emp_sem(it,:), ...
        'o-','LineWidth',1.2, ...
        'DisplayName','Empirique complet \pm SEM');

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

% 5) beta0 total(t) avec incertitude empirique
figure;
hold on;

plot(time_values,beta0_total_th_t, ...
    'LineWidth',2, ...
    'DisplayName','Nouvelle theorie PPP');

plot(time_values,beta0_total_emp_t, ...
    'LineWidth',1.5, ...
    'DisplayName','Empirique complet');

plot(time_values, ...
    beta0_total_emp_t+beta0_total_emp_sem, ...
    ':','LineWidth',1, ...
    'DisplayName','Empirique + SEM');

plot(time_values, ...
    beta0_total_emp_t-beta0_total_emp_sem, ...
    ':','LineWidth',1, ...
    'DisplayName','Empirique - SEM');

grid on;
xlabel('Temps (s)');
ylabel('\beta_0(t)');
title('Comparaison de \beta_0(t) theorique et empirique complet');
legend('Location','best');
hold off;

% 6) Moyenne temporelle selon la latitude
figure;
hold on;

plot(rad2deg(phi_vals), ...
    mean(beta0_bin_th,1,'omitnan'), ...
    'LineWidth',2, ...
    'DisplayName','Nouvelle theorie PPP');

errorbar(rad2deg(phi_vals), ...
    mean(beta0_bin_emp,1,'omitnan'), ...
    std(beta0_bin_emp,0,1,'omitnan')/sqrt(Nt), ...
    'o-','LineWidth',1.2, ...
    'DisplayName','Empirique moyen temporel');

grid on;
xlabel('Latitude \phi (deg)');
ylabel('\beta_0 moyen par tranche');
title('Moyenne temporelle de \beta_0(t,\phi)');
legend('Location','best');
hold off;

% 7) Erreurs selon la latitude
figure;
hold on;

plot(rad2deg(phi_vals),rmse_by_phi, ...
    'LineWidth',2, ...
    'DisplayName','RMSE');

plot(rad2deg(phi_vals),mae_by_phi, ...
    '--','LineWidth',1.8, ...
    'DisplayName','MAE');

plot(rad2deg(phi_vals),bias_by_phi, ...
    ':','LineWidth',1.8, ...
    'DisplayName','Biais emp-th');

grid on;
xlabel('Latitude \phi (deg)');
ylabel('Erreur sur \beta_0 par tranche');
title('Erreur de la nouvelle theorie selon la latitude');
legend('Location','best');
hold off;

% 8) ESS minimal par taille
ESS_min_by_s = nan(Smax,1);

for s = 2:Smax
    tmp = squeeze(ESS_t_phi(s,:,:));
    ESS_min_by_s(s) = min(tmp,[],'all','omitnan');
end

figure;
plot(2:Smax,ESS_min_by_s(2:end), ...
    'o-','LineWidth',1.4);
grid on;
xlabel('Taille s');
ylabel('ESS minimal sur (t,\phi)');
title('Qualite numerique de l''importance sampling');

%% ============================================================
%  12. SAUVEGARDE
%% ============================================================

output_file = fullfile(script_dir,'beta0_t_phi_results.mat');

save(output_file, ...
    'lambda_file','param_file', ...
    'N','R','dmax','inc','omega', ...
    'time_values','phi_edges','phi_vals','dphi_bins','A_phi_bins', ...
    'lambda_t_phi', ...
    'A_sphere','alpha_max','A_cap', ...
    'Smax','Nprobe','Nsamp','rho_tilt','beta_shape', ...
    'EC_density_t_phi','ESS_t_phi', ...
    'beta0_density_t_phi','beta0_bin_th','beta0_per_rad_t_phi', ...
    'beta0_total_th_t','beta0_total_th_mean', ...
    'beta0_density_phi_mean_t','beta0_per_rad_phi_mean_t', ...
    'compare_empirical','n_iterations_beta0_emp','rng_seed_beta0', ...
    'beta0_emp_iterations', ...
    'beta0_bin_emp','beta0_bin_emp_std','beta0_bin_emp_sem', ...
    'beta0_per_rad_emp','beta0_per_rad_emp_sem', ...
    'beta0_total_emp_iterations','beta0_total_emp_t', ...
    'beta0_total_emp_std','beta0_total_emp_sem', ...
    'beta0_total_emp_from_local_iterations', ...
    'max_local_global_error_emp', ...
    'difference_bin','difference_total', ...
    'rmse_t_phi','mae_t_phi','bias_t_phi', ...
    'rmse_by_phi','mae_by_phi','bias_by_phi', ...
    'rmse_total_t','mae_total_t','bias_total_t', ...
    '-v7.3');

fprintf('Resultats sauvegardes dans %s\n',output_file);

%% ============================================================
%  FONCTIONS LOCALES
%% ============================================================

function edges = random_labeled_tree_prufer(n)

    if n==2
        edges=[1 2];
        return;
    end

    P=randi(n,n-2,1);

    deg=ones(n,1);
    for k=1:numel(P)
        deg(P(k))=deg(P(k))+1;
    end

    edges=zeros(n-1,2);

    for k=1:n-2
        leaf=find(deg==1,1,'first');
        v=P(k);

        edges(k,:)=[leaf v];

        deg(leaf)=deg(leaf)-1;
        deg(v)=deg(v)-1;
    end

    last=find(deg==1);
    edges(n-1,:)=last(1:2).';
end

function [parent,order]=root_tree(edges,n,root)

    A=false(n,n);

    for k=1:size(edges,1)
        i=edges(k,1);
        j=edges(k,2);

        A(i,j)=true;
        A(j,i)=true;
    end

    parent=zeros(n,1);
    visited=false(n,1);
    order=zeros(n,1);

    queue=zeros(n,1);
    head=1;
    tail=1;

    queue(1)=root;
    visited(root)=true;

    no=0;

    while head<=tail

        v=queue(head);
        head=head+1;

        no=no+1;
        order(no)=v;

        neigh=find(A(v,:));

        for u=neigh
            if ~visited(u)
                visited(u)=true;
                parent(u)=v;

                tail=tail+1;
                queue(tail)=u;
            end
        end
    end
end

function x=sample_cap_around_tilted(c,alpha,rho,a)

    c=c(:).';
    c=c/norm(c);

    if rand<rho
        u=rand^(1/a);
    else
        u=rand;
    end

    cos_theta=1-u*(1-cos(alpha));
    sin_theta=sqrt(max(0,1-cos_theta^2));

    az=2*pi*rand;

    if abs(c(3))<0.9
        ref=[0 0 1];
    else
        ref=[1 0 0];
    end

    e1=cross(ref,c);
    e1=e1/norm(e1);

    e2=cross(c,e1);
    e2=e2/norm(e2);

    x=cos_theta*c ...
      + sin_theta*cos(az)*e1 ...
      + sin_theta*sin(az)*e2;

    x=x/norm(x);
end

function log_tau=log_weighted_spanning_tree_sum(W)

    n=size(W,1);

    if n==1
        log_tau=0;
        return;
    end

    d=sum(W,2);
    L=diag(d)-W;

    M=L(2:end,2:end);
    M=(M+M.')/2;

    [Rchol,p]=chol(M);

    if p==0
        log_tau=2*sum(log(diag(Rchol)));
    else
        ev=real(eig(M));

        tol=max(1e-14,1e-12*max(abs(ev)));
        ev(ev<tol)=tol;

        log_tau=sum(log(ev));
    end
end

function tf=is_connected_adj(A)

    n=size(A,1);
    seen=false(n,1);

    queue=zeros(n,1);
    head=1;
    tail=1;

    queue(1)=1;
    seen(1)=true;

    while head<=tail

        v=queue(head);
        head=head+1;

        neigh=find(A(v,:) & ~seen.');

        for u=neigh
            seen(u)=true;

            tail=tail+1;
            queue(tail)=u;
        end
    end

    tf=all(seen);
end

function X=fibonacci_sphere(n)

    k=(0:n-1).';
    z=1-2*(k+0.5)/n;

    golden_angle=pi*(3-sqrt(5));
    az=golden_angle*k;

    rxy=sqrt(max(0,1-z.^2));

    X=[rxy.*cos(az), ...
       rxy.*sin(az), ...
       z];
end

function y=logsumexp_vec(x)

    xmax=max(x);
    y=xmax+log(sum(exp(x-xmax)));
end


function [u0,Omega] = sample_initial_orbits_spatial(N,inc)
% Generation identique a betti_t_phi.m :
% sin(phi_0) uniforme sur [-sin(i),sin(i)] afin d'obtenir une
% distribution spatiale uniforme dans la bande au temps initial.

    sin_phi0 = -sin(inc)+2*sin(inc)*rand(N,1);

    sin_u0 = sin_phi0/sin(inc);
    sin_u0 = min(max(sin_u0,-1),1);

    u_base = asin(sin_u0);
    ascending = rand(N,1)<0.5;

    u0 = zeros(N,1);
    u0(ascending) = u_base(ascending);
    u0(~ascending) = pi-u_base(~ascending);

    u0 = mod(u0,2*pi);

    Omega = 2*pi*rand(N,1);
end

function X = positions_from_orbits(u0,Omega,t,R,inc,omega)

    u = mod(u0+omega*t,2*pi);

    cO = cos(Omega);
    sO = sin(Omega);
    cu = cos(u);
    su = sin(u);

    X = R*[ ...
        cO.*cu-sO.*su*cos(inc), ...
        sO.*cu+cO.*su*cos(inc), ...
        su*sin(inc)];
end

function A = adjacency_from_positions(X,dmax)

    gram = X*X.';
    norm2 = sum(X.^2,2);

    D2 = max(norm2+norm2.'-2*gram,0);

    A = D2<=dmax^2;
    A(1:size(A,1)+1:end) = false;

    A = triu(A,1);
    A = A|A.';
end
