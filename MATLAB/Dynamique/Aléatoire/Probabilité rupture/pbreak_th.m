clear; clc; close all;

%% ============================================================
%  CALCUL THEORIQUE DE p_break
%  Modele aleatoire a vecteurs tangentiels
%
%  Entree :
%    analysis_temp_results.mat
%
%  Sortie :
%    pbreak_th_results.mat
%% ============================================================

%% Chargement des paramètres
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

% Valeur empirique moyenne du nombre de composantes.
beta0_empirical = mean(double(S.beta0(:)),'omitnan');

if ~isfinite(beta0_empirical) || beta0_empirical <= 0
    error('La moyenne empirique de beta0 est invalide.');
end

% Chargement de la probabilité empirique qu'un lien de bord soit un pont.
p_bridge_emp_candidates = {
    fullfile(script_dir,'p_bridge_emp_results.mat')
    fullfile(script_dir,'..','p_bridge_emp_results.mat')
    fullfile(script_dir,'..','Paramètres','p_bridge_emp_results.mat')
};

p_bridge_emp_file = '';

for k = 1:numel(p_bridge_emp_candidates)
    if isfile(p_bridge_emp_candidates{k})
        p_bridge_emp_file = p_bridge_emp_candidates{k};
        break;
    end
end

if isempty(p_bridge_emp_file)
    error('Fichier p_bridge_emp_results.mat introuvable.');
end

Sbridge = load(p_bridge_emp_file);

% p_bridge_emp est la moyenne pondérée globale.
% On conserve une compatibilité avec les autres noms possibles.
if isfield(Sbridge,'p_bridge_emp')
    p_bridge_empirical = double(Sbridge.p_bridge_emp);
elseif isfield(Sbridge,'p_bridge_bord_mean_t')
    p_bridge_empirical = double(Sbridge.p_bridge_bord_mean_t);
elseif isfield(Sbridge,'p_bridge_bord_t')
    p_bridge_empirical = ...
        mean(double(Sbridge.p_bridge_bord_t(:)),'omitnan');
else
    error(['Aucune variable empirique reconnue pour p_bridge ', ...
           'dans %s.'],p_bridge_emp_file);
end

p_bridge_empirical = min(max(p_bridge_empirical,0),1);

%% ============================================================
%  1. Vitesse relative moyenne
%% ============================================================

mu = 398600;                    % km^3/s^2
omega = sqrt(mu / R^3);         % rad/s
v_orb = R * omega;              % km/s

% Approximation du modèle aléatoire
v_rel = (4/pi) * v_orb;

%% ============================================================
%  2. Probabilité de rupture d'un lien individuel
%
%  p_break_link
%    = (2/pi) * v_rel*dt/dmax
%    = (8/pi^2) * v_orb*dt/dmax
%% ============================================================

p_break_link = (2/pi) * (v_rel * dt / dmax);
p_break_link = min(max(p_break_link, 0), 1 - eps);

%% ============================================================
%  3. Probabilité de lien et nombre moyen d'arêtes
%% ============================================================

alpha_max = 2 * asin(min(dmax / (2*R), 1));

p_link = (1 - cos(alpha_max)) / 2;

E_theory = N * (N - 1) / 2 * p_link;

%% ============================================================
%  4. DEUX approximations théoriques de beta0
%
%  A) Ancienne approximation :
%       beta0 ≈ 2 + N1 + N2 + N3
%
%  B) Nouvelle approximation PPP par tailles de composantes :
%       beta0 ≈ sum_{s=1}^{Smax} E[C_s]
%
%  Les deux sont conservées séparément pour comparer leur impact
%  sur p_break.
%% ============================================================

% ---------- A) Ancienne approximation isolés/dimères/trimères ----------
c2_union = 1 + 3*sqrt(3)/(4*pi);
c3_conn  = 1 + 3*sqrt(3)/(2*pi);
c3_union = 1.80;

q1_ext = max(1 - p_link, 0);
q2_ext = max(1 - c2_union*p_link, 0);
q3_ext = max(1 - c3_union*p_link, 0);

N1_theory_old = N * q1_ext^(N - 1);

if N >= 2
    N2_theory_old = nchoosek(N, 2) ...
        * p_link ...
        * q2_ext^(N - 2);
else
    N2_theory_old = 0;
end

if N >= 3
    p_conn_3 = min(max(c3_conn * p_link^2, 0), 1);

    N3_theory_old = nchoosek(N, 3) ...
        * p_conn_3 ...
        * q3_ext^(N - 3);
else
    N3_theory_old = 0;
end

beta0_theory_old = ...
    2 + N1_theory_old + N2_theory_old + N3_theory_old;

beta0_theory_old = min(max(beta0_theory_old,1),N);

% ---------- B) Nouvelle approximation PPP/cluster ----------
beta0_file_candidates = {
    fullfile(script_dir,'beta0_th_results.mat')
    fullfile(script_dir,'..','beta0_th_results.mat')
};

beta0_file = '';
for k = 1:numel(beta0_file_candidates)
    if isfile(beta0_file_candidates{k})
        beta0_file = beta0_file_candidates{k};
        break;
    end
end

if isempty(beta0_file)
    error(['Fichier beta0_th_results.mat introuvable. ', ...
           'Exécuter beta0_th.m avant pbreak_th.m.']);
end

Sbeta = load(beta0_file);

if ~isfield(Sbeta,'beta0_th_trunc') || ~isfield(Sbeta,'EC')
    error(['beta0_th_results.mat doit contenir ', ...
           'beta0_th_trunc et EC.']);
end

beta0_theory_cluster = double(Sbeta.beta0_th_trunc);
EC_cluster = double(Sbeta.EC(:));

N1_theory_cluster = EC_cluster(1);

if isfield(Sbeta,'Smax')
    Smax_cluster = double(Sbeta.Smax);
else
    Smax_cluster = numel(EC_cluster);
end

beta0_theory_cluster = ...
    min(max(beta0_theory_cluster,1),N);

% Alias conservés uniquement pour les anciennes sections empiriques.
% Les calculs principaux ci-dessous utilisent explicitement _old et _cluster.
N1_theory = N1_theory_old;
N2_theory = N2_theory_old;
N3_theory = N3_theory_old;
beta0_theory = beta0_theory_old;

%% ============================================================
%  5. Probabilité théorique qu'un lien de la couronne soit un pont
%
%  Modèle 2 : moyenne conditionnelle sur la couronne de rupture
%
%  D appartient à [dmax-l_out_eff, dmax]
%
%  p_bridge_bord =
%    int exp[-lambda A_inter(D)] f_D(D) dD
%    ------------------------------------------------
%              int f_D(D) dD
%% ============================================================

% Largeur radiale sortante moyenne
l_out_eff = v_rel * dt / pi;
l_out_eff = min(max(l_out_eff, 0), dmax);

D_min = max(dmax - l_out_eff, 0);
D_max = dmax;

% Aire d'intersection de deux disques de rayon dmax
A_inter = @(D) ...
    2 * dmax^2 .* ...
    acos(min(max(D ./ (2*dmax), -1), 1)) ...
    - 0.5 .* D .* ...
    sqrt(max(4*dmax^2 - D.^2, 0));

% Densité radiale conditionnelle à l'existence d'un lien
f_D_given_link = @(D) 2 .* D ./ dmax^2;

if D_max > D_min
    numerator = integral( ...
        @(D) exp(-lambda .* A_inter(D)) .* f_D_given_link(D), ...
        D_min, D_max, ...
        'RelTol', 1e-8, ...
        'AbsTol', 1e-11);

    denominator = integral( ...
        f_D_given_link, ...
        D_min, D_max, ...
        'RelTol', 1e-10, ...
        'AbsTol', 1e-13);

    if denominator > 0
        p_bridge_bord = numerator / denominator;
    else
        p_bridge_bord = exp(-lambda * A_inter(dmax));
    end
else
    p_bridge_bord = exp(-lambda * A_inter(dmax));
end

p_bridge_bord = min(max(p_bridge_bord, 0), 1);

%% ============================================================
%  6. Nombre moyen de liens par composante non isolée
%     pour les DEUX approximations de beta0
%% ============================================================

% ---------- Ancienne approximation ----------
n_nonisolated_components_old = ...
    beta0_theory_old - N1_theory_old;

if n_nonisolated_components_old > 0
    mean_links_per_component_old = ...
        E_theory / n_nonisolated_components_old;
else
    mean_links_per_component_old = 0;
end

mean_breaking_bridges_per_component_old = ...
    mean_links_per_component_old * p_bridge_bord;

% ---------- Nouvelle approximation PPP/cluster ----------
n_nonisolated_components_cluster = ...
    beta0_theory_cluster - N1_theory_cluster;

if n_nonisolated_components_cluster > 0
    mean_links_per_component_cluster = ...
        E_theory / n_nonisolated_components_cluster;
else
    mean_links_per_component_cluster = 0;
end

mean_breaking_bridges_per_component_cluster = ...
    mean_links_per_component_cluster * p_bridge_bord;

%% ============================================================
%  7. Probabilité théorique de rupture d'une composante
%
%  Pour chaque approximation :
%
%    mu_break = Lbar_C * p_break_link * p_bridge_bord
%
%    p_break|non-isolee = 1-exp(-mu_break)
%
%  puis déconditionnement par la fraction de composantes non isolées.
%% ============================================================

% ---------- A) Ancienne approximation ----------
mu_break_old = ...
    mean_links_per_component_old ...
    * p_break_link ...
    * p_bridge_bord;

p_break_th_conditional_old = ...
    1 - exp(-mu_break_old);

fraction_nonisolated_components_old = ...
    n_nonisolated_components_old / max(beta0_theory_old,eps);

fraction_nonisolated_components_old = ...
    min(max(fraction_nonisolated_components_old,0),1);

p_break_th_old = ...
    fraction_nonisolated_components_old ...
    * p_break_th_conditional_old;

p_break_th_old = min(max(p_break_th_old,0),1-eps);

p_break_th_linear_old = ...
    fraction_nonisolated_components_old * mu_break_old;

p_break_th_linear_old = ...
    min(max(p_break_th_linear_old,0),1-eps);

% ---------- B) Nouvelle approximation PPP/cluster ----------
mu_break_cluster = ...
    mean_links_per_component_cluster ...
    * p_break_link ...
    * p_bridge_bord;

p_break_th_conditional_cluster = ...
    1 - exp(-mu_break_cluster);

fraction_nonisolated_components_cluster = ...
    n_nonisolated_components_cluster ...
    / max(beta0_theory_cluster,eps);

fraction_nonisolated_components_cluster = ...
    min(max(fraction_nonisolated_components_cluster,0),1);

p_break_th_cluster = ...
    fraction_nonisolated_components_cluster ...
    * p_break_th_conditional_cluster;

p_break_th_cluster = ...
    min(max(p_break_th_cluster,0),1-eps);

p_break_th_linear_cluster = ...
    fraction_nonisolated_components_cluster * mu_break_cluster;

p_break_th_linear_cluster = ...
    min(max(p_break_th_linear_cluster,0),1-eps);

% Aliases pour compatibilité avec les sections empiriques suivantes :
% on garde l'ancienne approximation comme référence historique.
n_nonisolated_components = n_nonisolated_components_old;
mean_links_per_component = mean_links_per_component_old;
mean_breaking_bridges_per_component = ...
    mean_breaking_bridges_per_component_old;
mu_break = mu_break_old;
fraction_nonisolated_components = ...
    fraction_nonisolated_components_old;
p_break_th_conditional = p_break_th_conditional_old;
p_break_th = p_break_th_old;
p_break_th_linear_conditional = mu_break_old;
p_break_th_linear = p_break_th_linear_old;

%% ============================================================
%  7.b Corrections par beta0 empirique et p_bridge empirique
%
% beta0 agit à deux endroits :
%   - dans le nombre de composantes non isolées ;
%   - dans le facteur de déconditionnement.
%
% Le remplacement exact est donc effectué dans la formule complète,
% plutôt que de multiplier directement p_break par un facteur unique.
%% ============================================================

% N1 reste ici théorique : seule l'approximation de beta0 est corrigée.
n_nonisolated_components_corrected = ...
    beta0_empirical - N1_theory;

n_nonisolated_components_corrected = ...
    max(n_nonisolated_components_corrected,0);

if n_nonisolated_components_corrected > 0
    mean_links_per_component_corrected = ...
        E_theory / n_nonisolated_components_corrected;
else
    mean_links_per_component_corrected = 0;
end

if beta0_empirical > 0
    fraction_nonisolated_components_corrected = ...
        n_nonisolated_components_corrected / beta0_empirical;
else
    fraction_nonisolated_components_corrected = 0;
end

fraction_nonisolated_components_corrected = ...
    min(max(fraction_nonisolated_components_corrected,0),1);

% Facteur correctif de beta0 dans l'exposant conditionnel :
% [E/(beta0_emp-N1)] / [E/(beta0_th-N1)].
if n_nonisolated_components_corrected > 0
    correction_beta0_exponent = ...
        n_nonisolated_components ...
        / n_nonisolated_components_corrected;
else
    correction_beta0_exponent = NaN;
end

if p_bridge_bord > 0
    correction_p_bridge = ...
        p_bridge_empirical / p_bridge_bord;
else
    correction_p_bridge = NaN;
end

correction_mu_break = ...
    correction_beta0_exponent * correction_p_bridge;

mu_break_corrected = ...
    mean_links_per_component_corrected ...
    * p_break_link ...
    * p_bridge_empirical;

p_break_th_conditional_corrected = ...
    1-exp(-mu_break_corrected);

p_break_th_conditional_corrected = ...
    min(max(p_break_th_conditional_corrected,0),1-eps);

p_break_th_corrected = ...
    fraction_nonisolated_components_corrected ...
    * p_break_th_conditional_corrected;

p_break_th_corrected = ...
    min(max(p_break_th_corrected,0),1-eps);

% Version linéaire corrigée, comparable au comptage moyen de ruptures.
p_break_th_linear_corrected = ...
    fraction_nonisolated_components_corrected ...
    * mu_break_corrected;

p_break_th_linear_corrected = ...
    min(max(p_break_th_linear_corrected,0),1-eps);

%% ============================================================
%  7.c Correction par le vrai flux empirique de ruptures
%
% Un lien rompu verifie :
%   A_ij(t) = 1, A_ij(t+dt) = 0.
%
% Flux theorique brut :
%
%   F_break^th = E_theory * p_break_link.
%
% Flux empirique :
%
%   F_break^emp = nombre moyen reel de liens rompus par pas.
%
% La version corrigee conserve beta0 empirique et p_bridge empirique.
%% ============================================================

if ~isfield(S,'Positions') || ~isfield(S,'Adjacency')
    error(['analysis_temp_results.mat doit contenir Positions et ', ...
           'Adjacency pour calculer le vrai flux de ruptures.']);
end

[n_broken_links_total,n_transitions_break_flux] = ...
    empirical_break_link_flux(S.Adjacency);

broken_links_per_step_emp = ...
    n_broken_links_total / max(n_transitions_break_flux,1);

broken_links_per_step_th = ...
    E_theory * p_break_link;

break_flux_ratio = ...
    broken_links_per_step_th / max(broken_links_per_step_emp,eps);

% Chaque arete appartient a une seule composante : il n'y a PAS
% de facteur 2 ici.
if n_nonisolated_components_corrected > 0
    mean_broken_links_per_component_emp = ...
        broken_links_per_step_emp ...
        / n_nonisolated_components_corrected;
else
    mean_broken_links_per_component_emp = 0;
end

mu_break_flux_emp = ...
    mean_broken_links_per_component_emp ...
    * p_bridge_empirical;

p_break_th_conditional_flux_emp = ...
    1-exp(-mu_break_flux_emp);

p_break_th_conditional_flux_emp = ...
    min(max(p_break_th_conditional_flux_emp,0),1-eps);

p_break_th_flux_emp = ...
    fraction_nonisolated_components_corrected ...
    * p_break_th_conditional_flux_emp;

p_break_th_flux_emp = ...
    min(max(p_break_th_flux_emp,0),1-eps);

%% ============================================================
%  8. Affichage
%% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' CALCUL THEORIQUE DE p_break - MODELE p_bridge_bord MOYENNE SUR LA COURONNE\n');
fprintf('============================================================\n');

fprintf('N                                   : %d\n', N);
fprintf('R                                   : %.6f km\n', R);
fprintf('lambda                              : %.8e sat/km^2\n', lambda);
fprintf('dmax                                : %.6f km\n', dmax);
fprintf('dt                                  : %.6f s\n', dt);

fprintf('------------------------------------------------------------\n');

fprintf('v_orb                               : %.8f km/s\n', v_orb);
fprintf('v_rel = 4/pi v_orb                  : %.8f km/s\n', v_rel);
fprintf('alpha_max                           : %.8f rad\n', alpha_max);
fprintf('p_link                              : %.8e\n', p_link);
fprintf('E[|E|]                              : %.8f\n', E_theory);

fprintf('------------------------------------------------------------\n');

fprintf('E[N1]                               : %.8f\n', N1_theory);
fprintf('E[N2]                               : %.8f\n', N2_theory);
fprintf('E[N3]                               : %.8f\n', N3_theory);
fprintf('E[beta0]                            : %.8f\n', beta0_theory);
fprintf('E[beta0]-E[N1]                      : %.8f\n', ...
    n_nonisolated_components);

fprintf('------------------------------------------------------------\n');

fprintf('l_out_eff                           : %.8f km\n', l_out_eff);
fprintf('Couronne [D_min,D_max]              : [%.8f, %.8f] km\n', D_min, D_max);
fprintf('p_bridge_bord                       : %.8f\n', p_bridge_bord);
fprintf('Liens moyens par composante         : %.8f\n', ...
    mean_links_per_component);
fprintf('Ponts de bord moyens/composante     : %.8f\n', ...
    mean_breaking_bridges_per_component);
fprintf('p_break d''un lien                   : %.8f\n', ...
    p_break_link);
fprintf('mu_break                            : %.8f\n', ...
    mu_break);

fprintf('------------------------------------------------------------\n');

fprintf('Fraction de composantes non isolées : %.8f\n', ...
    fraction_nonisolated_components);
fprintf('p_break conditionnel non isolé      : %.8f\n', ...
    p_break_th_conditional);
fprintf('p_break probabiliste global         : %.8f\n', ...
    p_break_th);
fprintf('p_break linéaire conditionnel       : %.8f\n', ...
    p_break_th_linear_conditional);
fprintf('p_break linéaire global             : %.8f\n', ...
    p_break_th_linear);

fprintf('------------------------------------------------------------\n');
fprintf(' COMPARAISON DES DEUX APPROXIMATIONS DE beta0\n');
fprintf('------------------------------------------------------------\n');
fprintf('Ancien beta0 (N1+N2+N3+2)           : %.8f\n', beta0_theory_old);
fprintf('Nouveau beta0 PPP (s <= %d)          : %.8f\n', ...
    Smax_cluster, beta0_theory_cluster);
fprintf('p_break ancien beta0                 : %.8f\n', p_break_th_old);
fprintf('p_break nouveau beta0 PPP            : %.8f\n', p_break_th_cluster);
fprintf('Ecart nouveau - ancien               : %+.8f\n', ...
    p_break_th_cluster-p_break_th_old);
fprintf('p_break lineaire ancien              : %.8f\n', ...
    p_break_th_linear_old);
fprintf('p_break lineaire nouveau             : %.8f\n', ...
    p_break_th_linear_cluster);

fprintf('------------------------------------------------------------\n');
fprintf(' CORRECTIONS AVEC VALEURS EMPIRIQUES\n');
fprintf('------------------------------------------------------------\n');
fprintf('beta0 empirique moyen               : %.8f\n', ...
    beta0_empirical);
fprintf('p_bridge empirique                  : %.8f\n', ...
    p_bridge_empirical);
fprintf('Composantes non isolées corrigées   : %.8f\n', ...
    n_nonisolated_components_corrected);
fprintf('Fraction non isolée corrigée        : %.8f\n', ...
    fraction_nonisolated_components_corrected);
fprintf('Facteur beta0 dans l''exposant       : %.8f\n', ...
    correction_beta0_exponent);
fprintf('Facteur correctif p_bridge          : %.8f\n', ...
    correction_p_bridge);
fprintf('Facteur total sur mu_break          : %.8f\n', ...
    correction_mu_break);
fprintf('mu_break corrigé                    : %.8f\n', ...
    mu_break_corrected);
fprintf('p_break conditionnel corrigé        : %.8f\n', ...
    p_break_th_conditional_corrected);
fprintf('p_break global corrigé              : %.8f\n', ...
    p_break_th_corrected);
fprintf('p_break linéaire global corrigé     : %.8f\n', ...
    p_break_th_linear_corrected);

fprintf('------------------------------------------------------------\n');
fprintf(' DIAGNOSTIC DU FLUX DE RUPTURES\n');
fprintf('------------------------------------------------------------\n');
fprintf('Liens rompus theo / pas             : %.8f\n', ...
    broken_links_per_step_th);
fprintf('Liens rompus emp / pas              : %.8f\n', ...
    broken_links_per_step_emp);
fprintf('Rapport flux theo / emp             : %.8f\n', ...
    break_flux_ratio);
fprintf('p_break vrai flux + beta0 + pbridge : %.8f\n', ...
    p_break_th_flux_emp);

fprintf('============================================================\n');

%% ============================================================
%  9. Sauvegarde
%% ============================================================

output_file = 'pbreak_th_results.mat';

save(output_file, ...
    'N', 'R', 'lambda', 'dmax', 'dt', ...
    'mu', 'omega', 'v_orb', 'v_rel', ...
    'alpha_max', 'p_link', 'E_theory', ...
    'c2_union', 'c3_conn', 'c3_union', ...
    'N1_theory', 'N2_theory', 'N3_theory', ...
    'beta0_theory', ...
    'beta0_file', 'Smax_cluster', 'EC_cluster', ...
    'beta0_theory_old', 'beta0_theory_cluster', ...
    'N1_theory_old', 'N2_theory_old', 'N3_theory_old', ...
    'N1_theory_cluster', ...
    'n_nonisolated_components_old', ...
    'n_nonisolated_components_cluster', ...
    'mean_links_per_component_old', ...
    'mean_links_per_component_cluster', ...
    'mu_break_old', 'mu_break_cluster', ...
    'fraction_nonisolated_components_old', ...
    'fraction_nonisolated_components_cluster', ...
    'p_break_th_conditional_old', ...
    'p_break_th_conditional_cluster', ...
    'p_break_th_old', 'p_break_th_cluster', ...
    'p_break_th_linear_old', ...
    'p_break_th_linear_cluster', ...
    'n_nonisolated_components', ...
    'l_out_eff', 'D_min', 'D_max', ...
    'p_bridge_bord', ...
    'mean_links_per_component', ...
    'mean_breaking_bridges_per_component', ...
    'p_break_link', ...
    'mu_break', ...
    'fraction_nonisolated_components', ...
    'p_break_th_conditional', ...
    'p_break_th', ...
    'p_break_th_linear_conditional', ...
    'p_break_th_linear', ...
    'beta0_empirical', ...
    'p_bridge_empirical', ...
    'p_bridge_emp_file', ...
    'n_nonisolated_components_corrected', ...
    'mean_links_per_component_corrected', ...
    'fraction_nonisolated_components_corrected', ...
    'correction_beta0_exponent', ...
    'correction_p_bridge', ...
    'correction_mu_break', ...
    'mu_break_corrected', ...
    'p_break_th_conditional_corrected', ...
    'p_break_th_corrected', ...
    'p_break_th_linear_corrected', ...
    'n_broken_links_total','n_transitions_break_flux', ...
    'broken_links_per_step_emp','broken_links_per_step_th', ...
    'break_flux_ratio', ...
    'mean_broken_links_per_component_emp', ...
    'mu_break_flux_emp', ...
    'p_break_th_conditional_flux_emp', ...
    'p_break_th_flux_emp');

fprintf('Résultats sauvegardés dans %s\n', output_file);

%% ============================================================
%  Fonction locale : flux empirique de ruptures de liens
%% ============================================================
function [n_broken_total,n_transitions] = ...
    empirical_break_link_flux(Adjacency)

    Nt = numel(Adjacency);
    n_broken_total = 0;
    n_transitions = 0;

    for t = 1:Nt-1
        A0 = Adjacency{t};
        A1 = Adjacency{t+1};

        if isempty(A0) || isempty(A1)
            continue;
        end

        A0 = logical(A0);
        A1 = logical(A1);

        broken_edges = triu(A0 & ~A1,1);

        n_broken_total = ...
            n_broken_total + nnz(broken_edges);

        n_transitions = n_transitions + 1;
    end
end
