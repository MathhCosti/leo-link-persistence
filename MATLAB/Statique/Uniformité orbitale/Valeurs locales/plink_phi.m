%% plink_phi.m
% p_link(phi) theorique et empirique par Monte-Carlo
% Walker Delta a uniformite orbitale.
%
% Le nombre de satellites N est charge depuis analysis_temp_results.mat.
% La meme valeur fixe est utilisee dans toutes les simulations
% Monte-Carlo et dans toutes les formules theoriques.
%
% Sortie :
%   plink_phi_results.mat
%
% Les resultats sauvegardes servent ensuite a edges_phi.m et degree_phi.m.

clear; clc; close all;

%% ============================================================
%  1. Parametres
%% ============================================================

R = 6371 + 550;       % km
inc_deg = 90;
inc = deg2rad(inc_deg);
dmax = 1500;          % km

n_phi = 200;          % nombre de tranches/points en latitude
nQuad = 300;          % quadrature pour la theorie

% Monte-Carlo
n_sim = 500;
rng_seed = 1;

rng(rng_seed);

%% Chargement de N depuis analysis_temp_results.mat
script_dir = fileparts(mfilename('fullpath'));

analysis_candidates = {
    fullfile(script_dir,'..', '..', '..', 'Dynamique', 'Walker Delta - Probabiliste', 'Uniformité orbitale', 'analysis_temp_results.mat')
};

analysis_file = '';

for k = 1:numel(analysis_candidates)
    if isfile(analysis_candidates{k})
        analysis_file = analysis_candidates{k};
        break;
    end
end

if isempty(analysis_file)
    error(['Fichier analysis_temp_results.mat introuvable dans le ', ...
           'dossier du script ou dans son dossier parent.']);
end

analysis_data = load(analysis_file,'N');

if ~isfield(analysis_data,'N')
    error('analysis_temp_results.mat ne contient pas la variable N.');
end

N = double(analysis_data.N);

if ~isscalar(N) || ~isfinite(N) || N < 2
    error('La variable N chargee doit etre un scalaire >= 2.');
end

N = round(N);

%% ============================================================
%  2. Grille de latitude
%% ============================================================

eps_phi = 1e-6;
phi_vals = linspace(-inc+eps_phi,inc-eps_phi,n_phi);
phi_vals = phi_vals(:).';

% Bords des tranches associees aux centres phi_vals
phi_edges = zeros(1,n_phi+1);
phi_edges(1) = -inc;
phi_edges(end) = inc;
phi_edges(2:end-1) = ...
    (phi_vals(1:end-1)+phi_vals(2:end))/2;

dphi_bins = diff(phi_edges);

alpha_max = 2*asin(min(dmax/(2*R),1));
cmax = cos(alpha_max);

%% ============================================================
%  3. p_link(phi) theorique
%% ============================================================

p_link_phi = zeros(size(phi_vals));

[u2,w2] = gauss_legendre_interval(nQuad,0,2*pi);
u2 = u2(:);
w2 = w2(:);

c2 = cos(u2);
s2 = sin(u2);
ci = cos(inc);
si = sin(inc);

for k = 1:numel(phi_vals)

    phi = phi_vals(k);

    s1 = sin(phi)/si;
    s1 = max(min(s1,1),-1);

    % Deux branches : passage ascendant et passage descendant
    u1a = asin(s1);
    u1b = pi-u1a;

    p_branch = zeros(2,1);
    u1_list = [u1a,u1b];

    for b = 1:2
        u1 = u1_list(b);
        c1 = cos(u1);
        s1b = sin(u1);

        A = c1.*c2 + ci^2.*s1b.*s2;
        B = ci.*(s1b.*c2-c1.*s2);
        C = si^2.*s1b.*s2;

        rho = sqrt(A.^2+B.^2);

        g = zeros(size(rho));
        mask = rho > 1e-14;

        q = zeros(size(rho));
        q(mask) = (cmax-C(mask))./rho(mask);

        g(mask & q <= -1) = 1;

        middle = mask & q > -1 & q < 1;
        g(middle) = acos(q(middle))/pi;

        g(~mask) = double(C(~mask) >= cmax);

        p_branch(b) = sum(w2.*g)/(2*pi);
    end

    p_link_phi(k) = mean(p_branch);
end

% Densite theorique de latitude.
% Ne pas la renormaliser par trapz : elle possede des singularites
% integrables en +/-inc, que la grille uniforme surestime fortement.
f_phi = cos(phi_vals) ./ ...
    (pi*sqrt(max(si^2-sin(phi_vals).^2,eps)));

% Probabilite exacte de chaque tranche de latitude :
%
% P(phi in [a,b]) =
%   [asin(sin(phi)/sin(i))]_a^b / pi.
%
% Cette formule integre exactement la singularite de f_phi aux bords.
cdf_edges = asin(max(min(sin(phi_edges)/si,1),-1))/pi + 1/2;
phi_bin_probability = diff(cdf_edges);

% Reconstruction globale par somme sur les tranches.
% p_link_phi est evalue au centre de chaque tranche.
p_link_global_from_phi = ...
    sum(p_link_phi.*phi_bin_probability);

p_link_sphere = (1-cos(alpha_max))/2;

%% ============================================================
%  4. Monte-Carlo
%
%  Pour chaque simulation :
%   - utilisation du meme nombre fixe N charge depuis
%     analysis_temp_results.mat ;
%   - tirage de Omega et u uniformes ;
%   - construction du graphe de liaison ;
%   - mesure locale de p_link(phi), du degre et des aretes.
%% ============================================================

N_draws = zeros(n_sim,1);

% Sommes agregees, robustes lorsque certaines tranches sont vides
sum_degree_by_bin = zeros(1,n_phi);
sum_opportunities_by_bin = zeros(1,n_phi);
count_sat_by_bin = zeros(1,n_phi);

% Mesures par simulation pour ecarts-types et intervalles
p_link_phi_sim = nan(n_sim,n_phi);
degree_phi_sim = nan(n_sim,n_phi);
edges_per_bin_sim = zeros(n_sim,n_phi);
edges_total_sim = zeros(n_sim,1);

total_ordered_links = 0;
total_ordered_pairs = 0;

for s = 1:n_sim

    % Meme nombre de satellites que dans analysis_temp_results.mat
    Ns = N;
    N_draws(s) = Ns;

    Omega = 2*pi*rand(Ns,1);
    u = 2*pi*rand(Ns,1);

    cO = cos(Omega);
    sO = sin(Omega);
    cu = cos(u);
    su = sin(u);

    rhat = [ ...
        cO.*cu-sO.*su*cos(inc), ...
        sO.*cu+cO.*su*cos(inc), ...
        su*sin(inc)];

    % Produit scalaire = cos(angle geocentrique)
    gram = rhat*rhat.';
    adjacency = gram >= cmax;
    adjacency(1:Ns+1:end) = false;
    adjacency = adjacency | adjacency.';

    degree = sum(adjacency,2);
    latitude = asin(max(min(rhat(:,3),1),-1));
    bin_id = discretize(latitude,phi_edges);

    edges_total_sim(s) = nnz(triu(adjacency,1));

    total_ordered_links = total_ordered_links + sum(degree);
    total_ordered_pairs = total_ordered_pairs + Ns*(Ns-1);

    for b = 1:n_phi
        in_bin = bin_id == b;
        nb = nnz(in_bin);

        if nb == 0
            continue;
        end

        sum_deg = sum(degree(in_bin));

        count_sat_by_bin(b) = count_sat_by_bin(b)+nb;
        sum_degree_by_bin(b) = sum_degree_by_bin(b)+sum_deg;
        sum_opportunities_by_bin(b) = ...
            sum_opportunities_by_bin(b)+nb*(Ns-1);

        degree_phi_sim(s,b) = sum_deg/nb;
        p_link_phi_sim(s,b) = sum_deg/(nb*(Ns-1));

        % Chaque arete est partagee entre ses deux extremites.
        edges_per_bin_sim(s,b) = 0.5*sum_deg;
    end
end

%% ============================================================
%  5. Moyennes empiriques
%% ============================================================

% N reste exactement la valeur chargee depuis analysis_temp_results.mat.
N_std = std(N_draws);

p_link_phi_emp = safe_divide( ...
    sum_degree_by_bin,sum_opportunities_by_bin);

degree_phi_emp = safe_divide( ...
    sum_degree_by_bin,count_sat_by_bin);

% Nombre moyen d'aretes affecte a chaque tranche
edges_per_bin_emp = mean(edges_per_bin_sim,1);

% Densite empirique d'aretes par radian
edges_density_phi_emp = edges_per_bin_emp./dphi_bins;

p_link_global_emp = total_ordered_links/total_ordered_pairs;
degree_global_emp = total_ordered_links/sum(N_draws);
edges_total_emp = mean(edges_total_sim);

% Incertitudes Monte-Carlo sur les moyennes par simulation
p_link_phi_emp_std = std(p_link_phi_sim,0,1,'omitnan');
p_link_phi_emp_sem = p_link_phi_emp_std ./ ...
    sqrt(sum(isfinite(p_link_phi_sim),1));

degree_phi_emp_std = std(degree_phi_sim,0,1,'omitnan');
degree_phi_emp_sem = degree_phi_emp_std ./ ...
    sqrt(sum(isfinite(degree_phi_sim),1));

edges_per_bin_emp_std = std(edges_per_bin_sim,0,1);
edges_per_bin_emp_sem = edges_per_bin_emp_std/sqrt(n_sim);

%% ============================================================
%  6. Comparaison p_link theorique / empirique
%% ============================================================

valid = isfinite(p_link_phi_emp) & p_link_phi > 0;

rmse_p_link = sqrt(mean( ...
    (p_link_phi_emp(valid)-p_link_phi(valid)).^2));

relative_error_p_link_global = ...
    abs(p_link_global_emp-p_link_global_from_phi) ...
    / max(abs(p_link_global_from_phi),eps);

figure;
plot(rad2deg(phi_vals),p_link_phi,'LineWidth',2);
hold on;
plot(rad2deg(phi_vals),p_link_phi_emp,'--','LineWidth',1.8);
grid on;
xlabel('Latitude \phi (deg)');
ylabel('p_{link}(\phi)');
title('Probabilite locale de lien : theorie et Monte-Carlo');
legend('Theorie','Monte-Carlo','Location','best');
hold off;

figure;
yyaxis left
plot(rad2deg(phi_vals),p_link_phi,'LineWidth',2);
hold on;
plot(rad2deg(phi_vals),p_link_phi_emp,'--','LineWidth',1.8);
ylabel('p_{link}(\phi)');

yyaxis right
plot(rad2deg(phi_vals),f_phi,':','LineWidth',1.6);
ylabel('f_\Phi(\phi)');

grid on;
xlabel('Latitude \phi (deg)');
title('Probabilite locale de lien et densite de latitude');
legend('p_{link}^{th}','p_{link}^{emp}','f_\Phi', ...
    'Location','best');

fprintf('\n');
fprintf('============================================================\n');
fprintf(' p_link(phi) WALKER DELTA - THEORIE / MONTE-CARLO\n');
fprintf('============================================================\n');
fprintf('Nombre de simulations              : %d\n',n_sim);
fprintf('Fichier source de N                : %s\n',analysis_file);
fprintf('N fixe charge                      : %d\n',N);
fprintf('Ecart-type de N entre simulations  : %.8f\n',N_std);
fprintf('Somme des probabilites de tranches : %.12f\n', ...
    sum(phi_bin_probability));
fprintf('p_link global theorique            : %.10f\n', ...
    p_link_global_from_phi);
fprintf('p_link global empirique            : %.10f\n', ...
    p_link_global_emp);
fprintf('Erreur relative globale            : %.3e\n', ...
    relative_error_p_link_global);
fprintf('RMSE locale                        : %.3e\n',rmse_p_link);
fprintf('============================================================\n');

%% ============================================================
%  7. Sauvegarde
%% ============================================================

output_file = fullfile(fileparts(mfilename('fullpath')), ...
    'plink_phi_results.mat');

save(output_file, ...
    'R','inc_deg','inc','dmax','alpha_max','cmax', ...
    'n_phi','nQuad','phi_vals','phi_edges','dphi_bins', ...
    'phi_bin_probability', ...
    'p_link_phi','f_phi','p_link_global_from_phi','p_link_sphere', ...
    'analysis_file','n_sim','rng_seed','N_draws','N','N_std', ...
    'p_link_phi_emp','p_link_phi_sim', ...
    'p_link_phi_emp_std','p_link_phi_emp_sem', ...
    'degree_phi_emp','degree_phi_sim', ...
    'degree_phi_emp_std','degree_phi_emp_sem', ...
    'edges_per_bin_emp','edges_density_phi_emp', ...
    'edges_per_bin_sim','edges_per_bin_emp_std', ...
    'edges_per_bin_emp_sem','edges_total_sim', ...
    'p_link_global_emp','degree_global_emp','edges_total_emp', ...
    'count_sat_by_bin','sum_degree_by_bin', ...
    'sum_opportunities_by_bin', ...
    'rmse_p_link','relative_error_p_link_global');

fprintf('Resultats sauvegardes dans %s\n',output_file);

%% ============================================================
%  Fonctions locales
%% ============================================================

function [x,w] = gauss_legendre_interval(n,a,b)
    k = (1:n-1)';
    beta = k./sqrt(4*k.^2-1);
    J = diag(beta,1)+diag(beta,-1);
    [V,D] = eig(J);
    x0 = diag(D);
    [x0,idx] = sort(x0);
    V = V(:,idx);
    w0 = 2*(V(1,:).^2)';
    x = (b-a)/2*x0+(a+b)/2;
    w = (b-a)/2*w0;
end

function ratio = safe_divide(num,den)
    ratio = nan(size(num));
    valid = den > 0;
    ratio(valid) = num(valid)./den(valid);
end
