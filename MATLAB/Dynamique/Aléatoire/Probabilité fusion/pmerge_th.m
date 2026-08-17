clear; clc; close all;

%% ============================================================
%  CALCUL THEORIQUE DE p_merge
%  Modele aleatoire a vecteurs tangentiels
%
%  Corrections utilisees :
%
%    phi_sweep :
%      aire réellement nouvelle pour un satellite
%      -------------------------------------------
%      aire géométrique 2*dmax*v_rel*dt
%
%    eta_sweep :
%      fraction de cette aire nouvelle qui n'est pas déjà
%      couverte par un voisin de la composante.
%
%  Passage à l'échelle d'une composante :
%
%    Nbar_C = N / E[beta0]
%
%  avec beta0 approché par les composantes isolées,
%  les dimères et les trimères.
%
%  Entree :
%    analysis_temp_results.mat
%
%  Sortie :
%    pmerge_th_results.mat
%% ============================================================

%% Chargement des paramètres du modèle
script_dir = fileparts(mfilename('fullpath'));
input_file = fullfile(script_dir, '..', 'analysis_temp_results.mat');

if ~isfile(input_file)
    error('Fichier introuvable : %s', input_file);
end

S = load(input_file, ...
    'N', 'R', 'lambda', 'dmax', 'dt', 'beta0', ...
    'Positions', 'Adjacency');

N      = S.N;
R      = S.R;
lambda = S.lambda;
dmax   = S.dmax;
dt     = S.dt;

if ~isfield(S,'beta0')
    error('La variable beta0 est absente de %s.',input_file);
end

% Valeur empirique moyenne de beta0 sur toute la dynamique.
beta0_empirical = mean(double(S.beta0(:)),'omitnan');

% Chargement de la valeur empirique de eta_sweep.
eta_file_candidates = {
    fullfile(script_dir,'eta_sweep_emp_results.mat')
    fullfile(script_dir,'..','eta_sweep_emp_results.mat')
    fullfile(script_dir,'..','Paramètres','eta_sweep_emp_results.mat')
};

eta_file = '';

for k = 1:numel(eta_file_candidates)
    if isfile(eta_file_candidates{k})
        eta_file = eta_file_candidates{k};
        break;
    end
end

if isempty(eta_file)
    error('Fichier eta_sweep_emp_results.mat introuvable.');
end

Seta = load(eta_file);

if isfield(Seta,'eta_sweep_empi')
    eta_sweep_empirical = double(Seta.eta_sweep_empi);
elseif isfield(Seta,'eta_sweep_mean_t')
    eta_sweep_empirical = double(Seta.eta_sweep_mean_t);
elseif isfield(Seta,'eta_sweep_mean_satellite')
    eta_sweep_empirical = double(Seta.eta_sweep_mean_satellite);
else
    error(['Aucune variable empirique reconnue pour eta_sweep dans ', ...
           '%s.'],eta_file);
end

eta_sweep_empirical = min(max(eta_sweep_empirical,0),1);

%% ============================================================
%  1. Paramètres orbitaux et vitesse relative moyenne
%% ============================================================

mu = 398600;                    % km^3/s^2
omega = sqrt(mu / R^3);         % rad/s
v_orb = R * omega;              % km/s
v_rel = (4/pi) * v_orb;         % km/s

% Déplacement relatif moyen pendant un pas de temps
ell = v_rel * dt;

%% ============================================================
%  2. Aire géométrique balayée brute par un satellite
%% ============================================================

A_sweep_geom = 2 * dmax * ell;

p_merge_raw = 1 - exp(-lambda * A_sweep_geom);
p_merge_raw = min(max(p_merge_raw, 0), 1 - eps);

%% ============================================================
%  3. Probabilité de lien et nombre moyen théorique d'arêtes
%% ============================================================

alpha_max = 2 * asin(min(dmax / (2*R), 1));
p_link = (1 - cos(alpha_max)) / 2;

E_theory = N * (N - 1) / 2 * p_link;

%% ============================================================
%  4. Approximation théorique de beta0 jusqu'aux trimères
%
%  E[beta0] ≈ 1 + E[N1] + E[N2] + E[N3]
%% ============================================================

c2_union = 1 + 3*sqrt(3)/(4*pi);
c3_conn  = 1 + 3*sqrt(3)/(2*pi);
c3_union = 1.80;

q1_ext = max(1 - p_link, 0);
q2_ext = max(1 - c2_union*p_link, 0);
q3_ext = max(1 - c3_union*p_link, 0);

% Composantes isolées
N1_theory = N * q1_ext^(N - 1);

% Dimères
if N >= 2
    N2_theory = nchoosek(N, 2) ...
        * p_link ...
        * q2_ext^(N - 2);
else
    N2_theory = 0;
end

% Trimères
if N >= 3
    p_conn_3 = min(max(c3_conn * p_link^2, 0), 1);

    N3_theory = nchoosek(N, 3) ...
        * p_conn_3 ...
        * q3_ext^(N - 3);
else
    N3_theory = 0;
end

beta0_theory = ...
    1 + N1_theory + N2_theory + N3_theory;

beta0_theory = min(max(beta0_theory, 1), N);

%% ============================================================
%  5. Nombre moyen théorique de satellites par composante
%
%  Nbar_C = N / E[beta0]
%% ============================================================

if beta0_theory > 0
    mean_satellites_per_component = N / beta0_theory;
else
    mean_satellites_per_component = 0;
end

%% ============================================================
%  6. Facteur théorique phi_sweep
%
%  Aire nouvelle exacte entre deux disques plans de rayon dmax
%  séparés de ell :
%
%    A_new_sat = pi*dmax^2 - A_inter(ell)
%
%  puis :
%
%    phi_sweep = A_new_sat / (2*dmax*ell)
%% ============================================================

if ell <= 0
    phi_sweep_th = 1;
    A_new_sat_th = 0;

elseif ell < 2*dmax
    A_inter_ell = ...
        2*dmax^2 * acos(ell/(2*dmax)) ...
        - 0.5*ell*sqrt(max(4*dmax^2 - ell^2, 0));

    A_new_sat_th = pi*dmax^2 - A_inter_ell;

    phi_sweep_th = ...
        A_new_sat_th / A_sweep_geom;

else
    % Si le déplacement dépasse le diamètre, les deux disques
    % ne se recouvrent plus.
    A_inter_ell = 0;
    A_new_sat_th = pi*dmax^2;

    phi_sweep_th = ...
        A_new_sat_th / A_sweep_geom;
end

phi_sweep_th = min(max(phi_sweep_th, 0), 1);

%% ============================================================
%  7. Facteur théorique eta_sweep
%
%  Approximation de bord :
%
%    eta_sweep ≈ exp[-lambda*A_inter(dmax)]
%
%  où A_inter(dmax) est l'aire d'intersection de deux disques
%  de rayon dmax dont les centres sont séparés de dmax.
%
%  Cette approximation considère qu'un point nouvellement exploré
%  est utile lorsqu'aucun voisin direct du satellite ne le couvrait.
%% ============================================================

A_inter_bord = ...
    (2*pi/3 - sqrt(3)/2) * dmax^2;

eta_sweep_th = ...
    exp(-lambda * A_inter_bord);

eta_sweep_th = min(max(eta_sweep_th, 0), 1);

%% ============================================================
%  8. Aire balayée efficace à l'échelle d'une composante
%
%  A_eff,C =
%      (E/beta0) * 2*dmax*v_rel*dt * phi_sweep * eta_sweep
%% ============================================================

A_sweep_corrected = ...
    mean_satellites_per_component ...
    * A_sweep_geom ...
    * phi_sweep_th ...
    * eta_sweep_th;

%% ============================================================
%  9. Probabilité théorique de fusion
%% ============================================================

p_merge_th = ...
    1 - exp(-lambda * A_sweep_corrected);

p_merge_th = min(max(p_merge_th, 0), 1 - eps);

p_disp_fusion_th = 0.5 * p_merge_th;

%% ============================================================
%  9.b Correction empirique de beta0 et eta_sweep
%% ============================================================

correction_beta0 = beta0_theory / beta0_empirical;
correction_eta_sweep = eta_sweep_empirical / eta_sweep_th;
correction_total = correction_beta0 * correction_eta_sweep;

merge_exponent_th = lambda * A_sweep_corrected;
merge_exponent_corrected = merge_exponent_th * correction_total;

p_merge_th_corrected = 1-exp(-merge_exponent_corrected);
p_merge_th_corrected = min(max(p_merge_th_corrected,0),1-eps);

p_disp_fusion_th_corrected = 0.5*p_merge_th_corrected;

%% ============================================================
%  9.c Correction par le vrai flux empirique de nouveaux liens
%
% Un nouveau lien verifie :
%   A_ij(t) = 0, A_ij(t+dt) = 1.
%
% Le flux theorique brut de nouveaux liens est :
%
%   F_new^th = (N/2) * lambda * A_sweep_geom * phi_sweep_th
%
% Le facteur eta_sweep n'est PAS inclus dans ce diagnostic :
% il intervient ensuite comme correction topologique.
%
% Deux versions sont calculees :
%   - vrai flux + beta0 empirique + eta_sweep empirique ;
%   - vrai flux + beta0 empirique
%       + vraie P(C_i ~= C_j | nouveau lien).
%% ============================================================

if ~isfield(S,'Positions') || ~isfield(S,'Adjacency')
    error(['analysis_temp_results.mat doit contenir Positions et ', ...
           'Adjacency pour calculer le vrai flux de nouveaux liens.']);
end

[n_new_links_total, n_merge_links_total, n_transitions_flux] = ...
    empirical_new_link_flux(S.Positions,S.Adjacency);

new_links_per_step_emp = ...
    n_new_links_total / max(n_transitions_flux,1);

merge_links_per_step_emp = ...
    n_merge_links_total / max(n_transitions_flux,1);

new_links_per_step_th = ...
    0.5 * N * lambda * A_sweep_geom * phi_sweep_th;

new_link_flux_ratio = ...
    new_links_per_step_th / max(new_links_per_step_emp,eps);

if n_new_links_total > 0
    p_diffcomp_given_new_emp = ...
        n_merge_links_total / n_new_links_total;
else
    p_diffcomp_given_new_emp = NaN;
end

% Chaque nouveau lien est incident a deux composantes lorsque ses
% extremites appartiennent a deux composantes distinctes. Pour passer
% d'un nombre global de nouveaux liens a un taux vu par composante,
% on utilise donc 2*F_new/beta0.
mu_new_per_component_emp = ...
    2 * new_links_per_step_emp / beta0_empirical;

% Version "flux empirique + eta empirique"
merge_exponent_flux_eta_emp = ...
    mu_new_per_component_emp * eta_sweep_empirical;

p_merge_th_flux_eta_emp = ...
    1-exp(-merge_exponent_flux_eta_emp);

p_merge_th_flux_eta_emp = ...
    min(max(p_merge_th_flux_eta_emp,0),1-eps);

% Version la plus corrigee :
% vrai flux + vraie probabilite qu'un nouveau lien relie
% deux composantes distinctes.
if isfinite(p_diffcomp_given_new_emp)
    merge_exponent_flux_true_emp = ...
        mu_new_per_component_emp * p_diffcomp_given_new_emp;

    p_merge_th_flux_true_emp = ...
        1-exp(-merge_exponent_flux_true_emp);

    p_merge_th_flux_true_emp = ...
        min(max(p_merge_th_flux_true_emp,0),1-eps);
else
    merge_exponent_flux_true_emp = NaN;
    p_merge_th_flux_true_emp = NaN;
end

%% Approximation linéaire
p_merge_th_linear = ...
    lambda * A_sweep_corrected;

p_merge_th_linear = ...
    min(max(p_merge_th_linear, 0), 1 - eps);

%% ============================================================
%  10. Affichage
%% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' CALCUL THEORIQUE DE p_merge - MODELE ALEATOIRE\n');
fprintf('============================================================\n');

fprintf('N                                   : %d\n', N);
fprintf('R                                   : %.6f km\n', R);
fprintf('lambda                              : %.8e sat/km^2\n', lambda);
fprintf('dmax                                : %.6f km\n', dmax);
fprintf('dt                                  : %.6f s\n', dt);

fprintf('------------------------------------------------------------\n');

fprintf('v_orb                               : %.8f km/s\n', v_orb);
fprintf('v_rel = 4/pi v_orb                  : %.8f km/s\n', v_rel);
fprintf('ell = v_rel*dt                      : %.8f km\n', ell);
fprintf('alpha_max                           : %.8f rad\n', alpha_max);
fprintf('p_link                              : %.8e\n', p_link);
fprintf('E[|E|]                              : %.8f\n', E_theory);

fprintf('------------------------------------------------------------\n');

fprintf('E[N1]                               : %.8f\n', N1_theory);
fprintf('E[N2]                               : %.8f\n', N2_theory);
fprintf('E[N3]                               : %.8f\n', N3_theory);
fprintf('E[beta0]                            : %.8f\n', beta0_theory);
fprintf('N/E[beta0]                          : %.8f\n', ...
    mean_satellites_per_component);

fprintf('------------------------------------------------------------\n');

fprintf('Aire balayée géométrique/satellite  : %.8f km^2\n', ...
    A_sweep_geom);
fprintf('Aire nouvelle exacte/satellite      : %.8f km^2\n', ...
    A_new_sat_th);
fprintf('phi_sweep théorique                 : %.8f\n', ...
    phi_sweep_th);
fprintf('A_inter au bord                     : %.8f km^2\n', ...
    A_inter_bord);
fprintf('eta_sweep théorique                 : %.8f\n', ...
    eta_sweep_th);
fprintf('phi_sweep * eta_sweep               : %.8f\n', ...
    phi_sweep_th * eta_sweep_th);

fprintf('------------------------------------------------------------\n');

fprintf('p_merge brute par satellite         : %.8f\n', ...
    p_merge_raw);
fprintf('Aire balayée corrigée/composante    : %.8f km^2\n', ...
    A_sweep_corrected);
fprintf('p_merge théorique probabiliste      : %.8f\n', ...
    p_merge_th);
fprintf('p_disp par fusion théorique         : %.8f\n', ...
    p_disp_fusion_th);

fprintf('------------------------------------------------------------\n');
fprintf('beta0 empirique moyen               : %.8f\n', ...
    beta0_empirical);
fprintf('eta_sweep empirique                 : %.8f\n', ...
    eta_sweep_empirical);
fprintf('Facteur correctif beta0             : %.8f\n', ...
    correction_beta0);
fprintf('Facteur correctif eta_sweep         : %.8f\n', ...
    correction_eta_sweep);
fprintf('Facteur correctif total             : %.8f\n', ...
    correction_total);
fprintf('p_merge avec corrections empiriques : %.8f\n', ...
    p_merge_th_corrected);
fprintf('p_disp fusion corrige empirique     : %.8f\n', ...
    p_disp_fusion_th_corrected);

fprintf('------------------------------------------------------------\n');
fprintf(' DIAGNOSTIC DU FLUX DE NOUVEAUX LIENS\n');
fprintf('------------------------------------------------------------\n');
fprintf('Nouveaux liens theo / pas           : %.8f\n', ...
    new_links_per_step_th);
fprintf('Nouveaux liens emp / pas            : %.8f\n', ...
    new_links_per_step_emp);
fprintf('Rapport flux theo / emp             : %.8f\n', ...
    new_link_flux_ratio);
fprintf('P(C_i ~= C_j | nouveau lien) emp    : %.8f\n', ...
    p_diffcomp_given_new_emp);
fprintf('p_merge vrai flux + eta emp         : %.8f\n', ...
    p_merge_th_flux_eta_emp);
fprintf('p_merge vrai flux + vraie proba     : %.8f\n', ...
    p_merge_th_flux_true_emp);

fprintf('------------------------------------------------------------\n');
fprintf('p_merge théorique linéaire          : %.8f\n', ...
    p_merge_th_linear);

fprintf('============================================================\n');

%% ============================================================
%  11. Sauvegarde
%% ============================================================

output_file = 'pmerge_th_results.mat';

save(output_file, ...
    'N', 'R', 'lambda', 'dmax', 'dt', ...
    'mu', 'omega', 'v_orb', 'v_rel', 'ell', ...
    'alpha_max', 'p_link', 'E_theory', ...
    'c2_union', 'c3_conn', 'c3_union', ...
    'N1_theory', 'N2_theory', 'N3_theory', ...
    'beta0_theory', ...
    'mean_satellites_per_component', ...
    'A_sweep_geom', 'A_new_sat_th', ...
    'phi_sweep_th', ...
    'A_inter_bord', 'eta_sweep_th', ...
    'A_sweep_corrected', ...
    'p_merge_raw', ...
    'p_merge_th', 'p_disp_fusion_th', ...
    'beta0_empirical', 'eta_sweep_empirical', ...
    'correction_beta0', 'correction_eta_sweep', ...
    'correction_total', ...
    'merge_exponent_th', 'merge_exponent_corrected', ...
    'p_merge_th_corrected', ...
    'p_disp_fusion_th_corrected', ...
    'n_new_links_total','n_merge_links_total','n_transitions_flux', ...
    'new_links_per_step_emp','merge_links_per_step_emp', ...
    'new_links_per_step_th','new_link_flux_ratio', ...
    'p_diffcomp_given_new_emp','mu_new_per_component_emp', ...
    'merge_exponent_flux_eta_emp','p_merge_th_flux_eta_emp', ...
    'merge_exponent_flux_true_emp','p_merge_th_flux_true_emp', ...
    'p_merge_th_linear', ...
    'eta_file');

fprintf('Résultats sauvegardés dans %s\n', output_file);


%% ============================================================
%  Fonction locale : flux empirique de nouveaux liens
%% ============================================================
function [n_new_total,n_merge_total,n_transitions] = ...
    empirical_new_link_flux(Positions,Adjacency)

    Nt = min(numel(Positions),numel(Adjacency));
    n_new_total = 0;
    n_merge_total = 0;
    n_transitions = 0;

    for t = 1:Nt-1
        A0 = Adjacency{t};
        A1 = Adjacency{t+1};

        if isempty(A0) || isempty(A1)
            continue;
        end

        A0 = logical(A0);
        A1 = logical(A1);

        new_edges = triu(A1 & ~A0,1);
        [ii,jj] = find(new_edges);

        n_transitions = n_transitions + 1;
        n_new_total = n_new_total + numel(ii);

        if isempty(ii)
            continue;
        end

        G0 = graph(sparse(A0),'upper');
        component_id = conncomp(G0).';

        n_merge_total = n_merge_total + ...
            nnz(component_id(ii) ~= component_id(jj));
    end
end
