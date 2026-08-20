%% beta0_phi_th_delta_orbital.m
% beta_0 theorique local en fonction de la latitude pour un Walker-Delta
% a uniformite orbitale.
%
% Approximation locale :
%   c_s(phi) = lambda(phi)^s/s! * s^(s-2) * A_cap^(s-1)
%              * E_q[ exp(-lambda(phi) A_union)/tau_h(G) ]
%
% avec c_s(phi) en composantes/km^2. Puis
%   beta0_density_phi = sum_s c_s(phi)
% et
%   dE[beta0]/dphi = 2*pi*R^2*cos(phi)*beta0_density_phi.
%
% Le Monte-Carlo sert uniquement a evaluer l'integrale theorique.

clear; clc; close all;
rng(1);

%% ===================== PARAMETRES WALKER DELTA =====================

% Le script cherche d'abord un fichier de densite deja produit.
density_candidates = {
    'densite_phi_results.mat'
    fullfile('..', '..', '..', 'Dynamique', 'Walker Delta - Probabiliste', 'Uniformité orbitale', 'Paramètres','densite_phi_results.mat')
};

density_file = '';
for k = 1:numel(density_candidates)
    if isfile(density_candidates{k})
        density_file = density_candidates{k};
        break;
    end
end

if ~isempty(density_file)
    D = load(density_file);

    if isfield(D,'N'), N = double(D.N); else, N = 204; end
    if isfield(D,'R'), R = double(D.R); else, R = 6371+550; end

    if isfield(D,'inc')
        inc = double(D.inc);
        inc_deg = rad2deg(inc);
    elseif isfield(D,'inc_deg')
        inc_deg = double(D.inc_deg);
        inc = deg2rad(inc_deg);
    else
        inc_deg = 58;
        inc = deg2rad(inc_deg);
    end
else
    R_earth = 6371;
    h = 550;
    R = R_earth+h;
    N = 204;
    inc_deg = 58;
    inc = deg2rad(inc_deg);
end

dmax = 1500;                  % km

%% ===================== PARAMETRES NUMERIQUES =====================

Smax = 50;
Nprobe = 4000;
Nsamp = 12000;

rho_tilt = 0.75;
beta_shape = 3.0;

% Grille de latitude sous forme de tranches.
% Les centres des tranches n'atteignent jamais exactement +/-inc :
% la singularite theorique de lambda(phi) est donc evitee sans supprimer
% de morceau de la bande orbitale.
n_phi = 120;
phi_edges = linspace(-inc,inc,n_phi+1);
phi = 0.5*(phi_edges(1:end-1)+phi_edges(2:end));
dphi_bins = diff(phi_edges);

% Aire exacte de chaque tranche sphérique [km^2].
A_phi_bins = ...
    2*pi*R^2 .* (sin(phi_edges(2:end))-sin(phi_edges(1:end-1)));

% Parametres de la comparaison empirique.
compare_empirical = true;
n_sim_emp = 500;
rng_seed_emp = 41;

%% ===================== DENSITE LOCALE =====================

% Walker Delta a uniformite orbitale :
den_phi = sqrt(max(sin(inc)^2-sin(phi).^2,eps));
lambda_phi = N ./ (2*pi^2*R^2.*den_phi);  % sat/km^2

surface_band = 4*pi*R^2*sin(inc);
lambda_band_mean = N/surface_band;

%% ===================== GEOMETRIE =====================

A_sphere = 4*pi*R^2;
alpha_max = 2*asin(min(1,dmax/(2*R)));
cos_alpha = cos(alpha_max);
one_minus_cos_alpha = 1-cos_alpha;
A_cap = 2*pi*R^2*one_minus_cos_alpha;

fprintf('============================================================\n');
fprintf('beta0(phi) theorique - Walker Delta uniformite orbitale\n');
fprintf('N                  = %d\n',N);
fprintf('R                  = %.3f km\n',R);
fprintf('inclinaison        = %.3f deg\n',inc_deg);
fprintf('dmax               = %.3f km\n',dmax);
fprintf('A_cap              = %.6e km^2\n',A_cap);
fprintf('Smax               = %d\n',Smax);
fprintf('Nsamp / taille     = %d\n',Nsamp);
fprintf('Nprobe             = %d\n',Nprobe);
fprintf('============================================================\n\n');

%% ===================== POINTS QMC POUR A_union =====================

probe_unit = fibonacci_sphere(Nprobe);

%% ===================== TABLEAUX =====================

% Densite locale de composantes de taille s [composantes/km^2]
EC_density_phi = zeros(Smax,n_phi);

% Contribution par radian de latitude [composantes/rad]
EC_per_rad_phi = zeros(Smax,n_phi);

% ESS local
ESS_phi = nan(Smax,n_phi);

%% ===================== s = 1 EXACT =====================

EC_density_phi(1,:) = lambda_phi .* exp(-lambda_phi*A_cap);
EC_per_rad_phi(1,:) = 2*pi*R^2*cos(phi).*EC_density_phi(1,:);
ESS_phi(1,:) = Inf;

fprintf('s =  1 : exact pour toutes les latitudes.\n');

%% ===================== s >= 2 =====================

x1 = [0 0 1];

for s = 2:Smax

    % Une seule banque de configurations geometriques pour toutes les phi.
    Aunion_all = zeros(Nsamp,1);
    logtau_all = zeros(Nsamp,1);

    for mm = 1:Nsamp

        % 1) Arbre etiquete uniforme de Cayley
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

        % 5) Matrice des poids h(u_ij)
        U = (1-Ddot)/one_minus_cos_alpha;
        U = max(0,min(1,U));

        H = zeros(s,s);
        mask = Adj;
        H(mask) = (1-rho_tilt) ...
            + rho_tilt*beta_shape.*U(mask).^(beta_shape-1);
        H = (H+H.')/2;
        H(1:s+1:end) = 0;

        % 6) Somme ponderee des arbres couvrants
        logtau_all(mm) = log_weighted_spanning_tree_sum(H);

        % 7) Aire de l'union des calottes
        maxdot = max(probe_unit*X.',[],2);
        frac_union = mean(maxdot >= cos_alpha);
        Aunion_all(mm) = A_sphere*frac_union;
    end

    % Partie du prefacteur independante de phi
    log_geom_prefactor = -gammaln(s+1) ...
        + (s-2)*log(s) ...
        + (s-1)*log(A_cap);

    for j = 1:n_phi
        lam = lambda_phi(j);

        % Poids local : exp[-lambda(phi) A_union]/tau_h(G)
        logw = -lam*Aunion_all-logtau_all;

        % Moyenne stable en log
        m = max(logw);
        log_mean_weight = m+log(mean(exp(logw-m)));

        log_cs = s*log(lam)+log_geom_prefactor+log_mean_weight;
        EC_density_phi(s,j) = exp(log_cs);

        % ESS local
        lsw = logsumexp_vec(logw);
        lsw2 = logsumexp_vec(2*logw);
        ESS_phi(s,j) = exp(2*lsw-lsw2);
    end

    EC_per_rad_phi(s,:) = 2*pi*R^2*cos(phi).*EC_density_phi(s,:);

    fprintf(['s = %2d : termine | ESS min/median = %.1f / %.1f ' ...
             '| max c_s = %.3e comp/km^2\n'], ...
        s,min(ESS_phi(s,:)),median(ESS_phi(s,:)), ...
        max(EC_density_phi(s,:)));
end

%% ===================== beta0(phi) =====================

% Densite locale du nombre de composantes
beta0_density_phi = sum(EC_density_phi,1);  % composantes/km^2

% Contribution au beta0 global par radian de latitude
dbeta0_dphi = sum(EC_per_rad_phi,1);        % composantes/rad

% Taille moyenne locale d'une composante (diagnostic)
mean_nodes_per_component_phi = ...
    lambda_phi ./ max(beta0_density_phi,eps);

%% ===================== INTEGRATION SUR LATITUDE =====================

% Integration par tranches sur toute la bande [-inc,+inc].
beta0_per_bin_th = beta0_density_phi .* A_phi_bins;
beta0_integrated_phi = sum(beta0_per_bin_th);

% Version "par radian", utile pour comparaison avec betti_phi.m.
beta0_density_per_rad_th = ...
    beta0_per_bin_th ./ dphi_bins;

N_per_bin_th = lambda_phi .* A_phi_bins;
N_reconstructed = sum(N_per_bin_th);

%% ===================== COMPARAISON EMPIRIQUE =====================
%
% Meme principe que dans betti_phi.m :
% chaque composante connexe de taille s contribue pour 1/s dans la
% tranche de latitude de chacun de ses s satellites.
%
% Ainsi, chaque composante contribue exactement pour 1 au total :
%
%   sum_b beta0_emp(b) = beta0
%
% pour chaque realisation.

beta0_per_bin_emp = nan(1,n_phi);
beta0_density_phi_emp = nan(1,n_phi);
beta0_density_per_rad_emp = nan(1,n_phi);
beta0_total_emp = NaN;
beta0_per_bin_emp_sem = nan(1,n_phi);
N_emp = NaN;

if compare_empirical

    rng(rng_seed_emp);

    beta0_per_bin_sim = zeros(n_sim_emp,n_phi);
    beta0_total_sim = zeros(n_sim_emp,1);
    N_draws_emp = zeros(n_sim_emp,1);

    cmax = 1-dmax^2/(2*R^2);

    for isim = 1:n_sim_emp

        Ns = max(draw_poisson(N),1);
        N_draws_emp(isim) = Ns;

        % Walker Delta a uniformite orbitale :
        % Omega et phase orbitale u uniformes.
        Omega = 2*pi*rand(Ns,1);
        uorb = 2*pi*rand(Ns,1);
        points = orbital_position(Omega,uorb,inc);

        gram = points*points.';
        adjacency = gram >= cmax;
        adjacency(1:Ns+1:end) = false;
        adjacency = adjacency | adjacency.';

        Gemp = graph(adjacency);
        labels = conncomp(Gemp).';
        component_sizes = accumarray(labels,1);

        beta0_total_sim(isim) = numel(component_sizes);

        latitudes = asin(max(min(points(:,3),1),-1));
        bin_id = discretize(latitudes,phi_edges);

        % Attribution symetrique des composantes aux latitudes.
        for c = 1:numel(component_sizes)

            members = find(labels==c);
            sc = numel(members);
            contribution = 1/sc;

            for r = 1:sc
                b = bin_id(members(r));

                if ~isnan(b)
                    beta0_per_bin_sim(isim,b) = ...
                        beta0_per_bin_sim(isim,b)+contribution;
                end
            end
        end
    end

    beta0_per_bin_emp = mean(beta0_per_bin_sim,1);
    beta0_per_bin_emp_std = std(beta0_per_bin_sim,0,1);
    beta0_per_bin_emp_sem = ...
        beta0_per_bin_emp_std/sqrt(n_sim_emp);

    % Deux normalisations possibles.
    beta0_density_phi_emp = ...
        beta0_per_bin_emp ./ A_phi_bins;           % composantes/km^2

    beta0_density_per_rad_emp = ...
        beta0_per_bin_emp ./ dphi_bins;            % composantes/rad

    beta0_total_emp = mean(beta0_total_sim);
    N_emp = mean(N_draws_emp);

    % Diagnostics.
    valid = isfinite(beta0_density_per_rad_th) ...
        & isfinite(beta0_density_per_rad_emp);

    rmse_per_rad = sqrt(mean( ...
        (beta0_density_per_rad_emp(valid) ...
        -beta0_density_per_rad_th(valid)).^2));

    relative_error_total = ...
        abs(beta0_total_emp-beta0_integrated_phi) ...
        / max(abs(beta0_total_emp),eps);

    ratio_emp_th_phi = nan(1,n_phi);
    mask_ratio = valid & beta0_density_per_rad_th>0;
    ratio_emp_th_phi(mask_ratio) = ...
        beta0_density_per_rad_emp(mask_ratio) ...
        ./ beta0_density_per_rad_th(mask_ratio);
else
    beta0_per_bin_sim = [];
    beta0_total_sim = [];
    N_draws_emp = [];
    beta0_per_bin_emp_std = nan(1,n_phi);
    rmse_per_rad = NaN;
    relative_error_total = NaN;
    ratio_emp_th_phi = nan(1,n_phi);
end

fprintf('\n============================================================\n');
fprintf('RESULTATS\n');
fprintf('beta0 theorique integre              : %.8f\n', ...
    beta0_integrated_phi);
fprintf('N reconstruit par lambda(phi)        : %.8f / %d\n', ...
    N_reconstructed,N);

if compare_empirical
    fprintf('beta0 empirique moyen                : %.8f\n', ...
        beta0_total_emp);
    fprintf('N empirique moyen                    : %.8f\n',N_emp);
    fprintf('Erreur relative totale beta0         : %.3f %%\n', ...
        100*relative_error_total);
    fprintf('RMSE locale (composantes/rad)        : %.6e\n', ...
        rmse_per_rad);
end
fprintf('============================================================\n');

%% ===================== GRAPHIQUES =====================

% 1) Comparaison principale en composantes par radian,
%    directement comparable au graphe de betti_phi.m.
figure;
plot(rad2deg(phi),beta0_density_per_rad_th,'LineWidth',2, ...
    'DisplayName','Theorie PPP par composantes');
hold on;

if compare_empirical
    plot(rad2deg(phi),beta0_density_per_rad_emp,'--','LineWidth',1.8, ...
        'DisplayName','Monte-Carlo direct');
end

grid on;
xlabel('Latitude \phi (deg)');
ylabel('Densite moyenne de composantes par radian');
title('\beta_0(\phi) : nouvelle theorie PPP et empirique');
legend('Location','best');
hold off;

% 2) Nombre moyen de composantes dans chaque tranche.
figure;
plot(rad2deg(phi),beta0_per_bin_th,'LineWidth',2, ...
    'DisplayName','Theorie PPP par composantes');
hold on;

if compare_empirical
    errorbar(rad2deg(phi),beta0_per_bin_emp,beta0_per_bin_emp_sem, ...
        '--','LineWidth',1.2,'DisplayName','Monte-Carlo direct');
end

grid on;
xlabel('Latitude \phi (deg)');
ylabel('Nombre moyen de composantes dans la tranche');
title('Nombre moyen local de composantes connexes');
legend('Location','best');
hold off;

% 3) Densite locale surfacique de composantes.
figure;
plot(rad2deg(phi),beta0_density_phi,'LineWidth',2, ...
    'DisplayName','Theorie PPP');
hold on;

if compare_empirical
    plot(rad2deg(phi),beta0_density_phi_emp,'--','LineWidth',1.8, ...
        'DisplayName','Monte-Carlo');
end

grid on;
xlabel('Latitude \phi (deg)');
ylabel('Densite locale de composantes (km^{-2})');
title('Densite surfacique locale de \beta_0');
legend('Location','best');
hold off;

% 4) Rapport empirique / theorie.
if compare_empirical
    figure;
    plot(rad2deg(phi),ratio_emp_th_phi,'LineWidth',1.8);
    hold on;
    yline(1,'--','Accord parfait');
    grid on;
    xlabel('Latitude \phi (deg)');
    ylabel('Empirique / theorie');
    title('Qualite locale de la nouvelle theorie de \beta_0');
    hold off;
end

% 5) lambda(phi)
figure;
plot(rad2deg(phi),lambda_phi,'LineWidth',2); hold on;
yline(lambda_band_mean,'--', ...
    sprintf('Moyenne bande = %.3e',lambda_band_mean), ...
    'LineWidth',1.3);
grid on;
xlabel('Latitude \phi (deg)');
ylabel('\lambda_{\Delta}(\phi) (satellites/km^2)');
title('Densite locale Walker Delta utilisee');
legend('\lambda_{\Delta}(\phi)','Moyenne bande','Location','best');
hold off;

% 6) Contributions des premieres tailles.
figure; hold on;
for s = 1:min(6,Smax)
    plot(rad2deg(phi),EC_density_phi(s,:).*A_phi_bins./dphi_bins, ...
        'LineWidth',1.4, ...
        'DisplayName',sprintf('C_%d',s));
end
grid on;
xlabel('Latitude \phi (deg)');
ylabel('Contribution theorique (composantes/rad)');
title('Contributions locales des premieres tailles de composantes');
legend('Location','best');
hold off;

% 7) ESS minimal selon s.
figure;
plot(2:Smax,min(ESS_phi(2:end,:),[],2),'o-','LineWidth',1.4);
grid on;
xlabel('Taille s');
ylabel('ESS minimal sur les latitudes');
title('Qualite numerique de l''importance sampling');

%% ===================== SAUVEGARDE =====================

save('beta0_phi_results.mat', ...
    'N','R','inc_deg','inc','dmax', ...
    'phi','phi_edges','dphi_bins','A_phi_bins', ...
    'lambda_phi','lambda_band_mean', ...
    'alpha_max','A_cap','A_sphere', ...
    'Smax','Nprobe','Nsamp','rho_tilt','beta_shape', ...
    'EC_density_phi','EC_per_rad_phi','ESS_phi', ...
    'beta0_density_phi','dbeta0_dphi', ...
    'beta0_per_bin_th','beta0_density_per_rad_th', ...
    'mean_nodes_per_component_phi', ...
    'beta0_integrated_phi','N_reconstructed', ...
    'compare_empirical','n_sim_emp','rng_seed_emp', ...
    'beta0_per_bin_emp','beta0_density_phi_emp', ...
    'beta0_density_per_rad_emp','beta0_per_bin_emp_sem', ...
    'beta0_per_bin_sim','beta0_total_emp','beta0_total_sim', ...
    'N_emp','N_draws_emp','rmse_per_rad','relative_error_total', ...
    'ratio_emp_th_phi');

fprintf('Resultats sauvegardes dans beta0_phi_th_results.mat\n');

%% ===================== FONCTIONS LOCALES =====================

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
        i=edges(k,1); j=edges(k,2);
        A(i,j)=true; A(j,i)=true;
    end

    parent=zeros(n,1);
    visited=false(n,1);
    order=zeros(n,1);
    queue=zeros(n,1);
    head=1; tail=1;
    queue(1)=root;
    visited(root)=true;
    no=0;

    while head<=tail
        v=queue(head); head=head+1;
        no=no+1; order(no)=v;
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

    e1=cross(ref,c); e1=e1/norm(e1);
    e2=cross(c,e1); e2=e2/norm(e2);

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
    head=1; tail=1;
    queue(1)=1; seen(1)=true;

    while head<=tail
        v=queue(head); head=head+1;
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
    X=[rxy.*cos(az),rxy.*sin(az),z];
end

function y=logsumexp_vec(x)
    xmax=max(x);
    y=xmax+log(sum(exp(x-xmax)));
end


function X = orbital_position(Omega,u,inc)
% Position unitaire d'un satellite Walker Delta a partir de
% l'ascension droite du noeud Omega et de la phase orbitale u.

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
% Tirage Poisson sans toolbox supplementaire.

    L = exp(-lambda);
    p = 1;
    k = 0;

    while p > L
        k = k+1;
        p = p*rand;
    end

    n = k-1;
end
