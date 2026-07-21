clear; clc; close all;

%% ============================================================
%  p_merge THEORIQUE - WALKER DELTA A UNIFORMITE ORBITALE
%
%  Modele local :
%
%  p_merge^Delta(phi)
%    = 1 - exp[-2 lambda^Delta(phi) d_max v_rel^Delta Delta_t chi^Delta]
%
%  avec
%
%  lambda^Delta(phi)
%    = N / [2*pi^2*R^2*sqrt(sin(i)^2 - sin(phi)^2)]
%
%  et
%
%  v_rel^Delta = v_orb/sqrt(2).
%
%  La probabilite globale est obtenue en moyennant en phase orbitale :
%
%  p_merge_bar^Delta
%    = (1/(2*pi)) int_0^(2*pi) p_merge^Delta(phi(u)) du,
%
%  phi(u) = asin(sin(i)*sin(u)).
%
%  L'integration en u evite la densite de latitude singuliere.
%% ============================================================

%% Parametres physiques
R_earth = 6371;          % km
h = 550;                 % km
R = R_earth + h;         % km

mu_earth = 398600;  % km^3/s^2
v_orb = sqrt(mu_earth/R);
v_rel = v_orb/sqrt(2);   % approximation retenue, km/s

%% Parametres de constellation
N = 204;
inc_deg = 90;
inc = deg2rad(inc_deg);

%% Parametres du modele de fusion
d_max = 1500;            % km
Delta_t = 20;             % s

%% Chargement des resultats analytiques externes
%
% Le facteur correctif est calcule a partir de :
%
%   chi_Delta = (N - E[beta_0^Delta]) / E[|E^Delta|]
%
% avec E[|E^Delta|] = C(N,2) p_link^Delta.

script_dir = fileparts(mfilename('fullpath'));
plink_results_file = fullfile(script_dir, '..', 'plink_walker_delta_semi_analytique_results.mat');
betti_results_file = fullfile(script_dir, '..', 'N1_N2_N3_walker_delta_results.mat');

if ~isfile(plink_results_file)
    error(['Fichier introuvable : %s\n' ...
        'Executer d''abord plink_quadrature.m.'],plink_results_file);
end
if ~isfile(betti_results_file)
    error(['Fichier introuvable : %s\n' ...
        'Executer d''abord betti_quadrature.m.'],betti_results_file);
end

plink_data = load(plink_results_file);
betti_data = load(betti_results_file);

if ~isfield(plink_data,'p_link_delta')
    error('La variable p_link_delta est absente de %s.',plink_results_file);
end
p_link_delta = plink_data.p_link_delta;

if isfield(plink_data,'E_delta')
    E_edges_delta = plink_data.E_delta;
else
    E_edges_delta = nchoosek(N,2)*p_link_delta;
end

if ~isfield(betti_data,'beta0_th_123')
    error('La variable beta0_th_123 est absente de %s.',betti_results_file);
end
beta0_delta = betti_data.beta0_th_123;

if isfield(plink_data,'N') && plink_data.N ~= N
    warning(['N differe entre pmerge_th (%d) et plink (%d). ' ...
        'Le nombre moyen d''aretes est recalcule avec N=%d.'], ...
        N,plink_data.N,N);
    E_edges_delta = nchoosek(N,2)*p_link_delta;
end
if isfield(betti_data,'N') && betti_data.N ~= N
    warning('N differe entre pmerge_th (%d) et betti (%d).',N,betti_data.N);
end
if isfield(plink_data,'dmax') && abs(plink_data.dmax-d_max) > 1e-12
    warning('d_max differe entre pmerge_th et plink_quadrature.');
end
if isfield(betti_data,'dmax') && abs(betti_data.dmax-d_max) > 1e-12
    warning('d_max differe entre pmerge_th et betti_quadrature.');
end
if isfield(plink_data,'inc_deg') && abs(plink_data.inc_deg-inc_deg) > 1e-12
    warning('Inclinaison differente dans plink_quadrature.');
end
if isfield(betti_data,'inc_deg') && abs(betti_data.inc_deg-inc_deg) > 1e-12
    warning('Inclinaison differente dans betti_quadrature.');
end

if E_edges_delta > 0
    chi_delta_raw = (N-beta0_delta)/E_edges_delta;
else
    chi_delta_raw = 0;
end
chi_delta = max(0,min(1,chi_delta_raw));
if abs(chi_delta-chi_delta_raw) > 1e-12
    warning(['chi_Delta brut = %.6f hors de [0,1]. ' ...
        'Valeur utilisee = %.6f.'],chi_delta_raw,chi_delta);
end

%% ============================================================
%  1. Densite locale et p_merge local pour le trace
%% ============================================================

n_phi = 2000;
eps_phi = 1e-6;

phi = linspace(-inc+eps_phi,inc-eps_phi,n_phi);

den_phi = sqrt(sin(inc)^2 - sin(phi).^2);

lambda_delta_phi = ...
    N ./ (2*pi^2*R^2.*den_phi);       % satellites/km^2

exponent_phi = ...
    2 .* lambda_delta_phi .* d_max .* v_rel .* Delta_t .* chi_delta;

p_merge_phi = 1-exp(-exponent_phi);

%% ============================================================
%  2. Moyenne orbitale theorique
%
%  Quadrature de Gauss-Legendre sur u dans [0,2*pi].
%% ============================================================

n_quad = 500;
[u_nodes,u_weights] = gauss_legendre(n_quad,0,2*pi);

phi_u = asin(sin(inc).*sin(u_nodes));

den_u_sq = sin(inc)^2 - sin(phi_u).^2;
den_u_sq = max(den_u_sq,0);
den_u = sqrt(den_u_sq);

% Aux points de retournement, lambda diverge mais p_merge tend vers 1.
lambda_u = zeros(size(u_nodes));
p_merge_u = zeros(size(u_nodes));

regular = den_u > 1e-13;

lambda_u(regular) = ...
    N ./ (2*pi^2*R^2.*den_u(regular));

p_merge_u(regular) = ...
    1-exp(-2.*lambda_u(regular).*d_max.*v_rel.*Delta_t.*chi_delta);

p_merge_u(~regular) = 1;

p_merge_mean = ...
    sum(u_weights.*p_merge_u)/(2*pi);

%% ============================================================
%  3. Approximation utilisant directement la densite moyenne
%% ============================================================

surface_band = 4*pi*R^2*sin(inc);
lambda_band_mean = N/surface_band;

p_merge_mean_density = ...
    1-exp(-2*lambda_band_mean*d_max*v_rel*Delta_t*chi_delta);

%% ============================================================
%  4. Traces
%% ============================================================

figure;
plot(rad2deg(phi),p_merge_phi,'LineWidth',2); hold on;

yline(p_merge_mean,'--', ...
    sprintf('Moyenne orbitale = %.4e',p_merge_mean), ...
    'LineWidth',1.5);

yline(p_merge_mean_density,':', ...
    sprintf('Avec densite moyenne = %.4e',p_merge_mean_density), ...
    'LineWidth',1.5);

grid on;
xlabel('Latitude \phi (deg)');
ylabel('p_{\mathrm{merge}}^{\Delta}(\phi)');
title(sprintf(['Probabilite theorique locale de fusion ' ...
    '(N=%d, i=%.1f deg, \\Delta t=%.2f s)'], ...
    N,inc_deg,Delta_t));

legend('p_{\mathrm{merge}}^{\Delta}(\phi)', ...
       'Moyenne orbitale', ...
       'Approximation par densite moyenne', ...
       'Location','best');

ylim([0, min(1,1.05*max([p_merge_phi,p_merge_mean_density]))]);
hold off;

figure;
plot(rad2deg(phi),lambda_delta_phi,'LineWidth',2);
grid on;
xlabel('Latitude \phi (deg)');
ylabel('\lambda_{\Delta}(\phi) (satellites/km^2)');
title('Densite locale utilisee dans p_{\mathrm{merge}}^{\Delta}');

%% ============================================================
%  5. Affichage console
%% ============================================================

fprintf('\n=== p_merge theorique Walker Delta ===\n');
fprintf('N                              : %d\n',N);
fprintf('Inclinaison                    : %.2f deg\n',inc_deg);
fprintf('Rayon orbital                  : %.2f km\n',R);
fprintf('Vitesse orbitale               : %.6f km/s\n',v_orb);
fprintf('Vitesse relative approximee    : %.6f km/s\n',v_rel);
fprintf('d_max                          : %.2f km\n',d_max);
fprintf('Delta_t                        : %.4f s\n',Delta_t);
fprintf('p_link_Delta                   : %.10f\n',p_link_delta);
fprintf('E[|E_Delta|]                   : %.6f\n',E_edges_delta);
fprintf('E[beta0_Delta] approxime       : %.6f\n',beta0_delta);
fprintf('chi_Delta brut                 : %.6f\n',chi_delta_raw);
fprintf('chi_Delta utilise              : %.6f\n',chi_delta);
fprintf('Densite moyenne sur la bande   : %.6e sat/km^2\n', ...
    lambda_band_mean);
fprintf('p_merge moyen, integration u   : %.10e\n',p_merge_mean);
fprintf('p_merge avec lambda moyenne    : %.10e\n',p_merge_mean_density);
fprintf('Ecart absolu                   : %.10e\n', ...
    p_merge_mean_density-p_merge_mean);
fprintf('Ecart relatif                  : %.4f %%\n', ...
    100*(p_merge_mean_density-p_merge_mean)/max(p_merge_mean,eps));

%% ============================================================
%  6. Sauvegarde
%% ============================================================

save('pmerge_theorique_walker_delta_results.mat', ...
    'R_earth','h','R','mu_earth','v_orb','v_rel', ...
    'N','inc_deg','inc','d_max','Delta_t', ...
    'p_link_delta','E_edges_delta','beta0_delta', ...
    'chi_delta_raw','chi_delta', ...
    'phi','lambda_delta_phi','p_merge_phi', ...
    'lambda_band_mean','p_merge_mean','p_merge_mean_density');

fprintf('\nResultats sauvegardes dans :\n');
fprintf('pmerge_theorique_walker_delta_results.mat\n');

%% ============================================================
%  FONCTION LOCALE : quadrature de Gauss-Legendre
%% ============================================================

function [x,w] = gauss_legendre(n,a,b)
    k = 1:n-1;
    beta = k ./ sqrt(4*k.^2-1);

    J = diag(beta,1) + diag(beta,-1);
    [V,D] = eig(J);

    x0 = diag(D);
    [x0,idx] = sort(x0);
    V = V(:,idx);

    w0 = 2*(V(1,:).^2)';

    x = (a+b)/2 + (b-a)/2*x0;
    w = (b-a)/2*w0;
end
