clear; clc; close all;

%% ============================================================
%  CALCUL THEORIQUE DE p_break
%  Modele aleatoire a vecteurs tangentiels
%
%  Entree :
%    leo_zigzag_analysis_random_vectors_results.mat
%
%  Sortie :
%    pbreak_th_random_results.mat
%% ============================================================

%% Chargement des paramètres
input_file = 'analysis_temp_results.mat';

if ~isfile(input_file)
    error('Fichier introuvable : %s', input_file);
end

S = load(input_file, 'N', 'R', 'lambda', 'dmax', 'dt');

N      = S.N;
R      = S.R;
lambda = S.lambda;
dmax   = S.dmax;
dt     = S.dt;

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
%  4. Approximation théorique de beta0
%
%  E[beta0] ≈ 1 + E[N1] + E[N2] + E[N3]
%% ============================================================

c2_union = 1 + 3*sqrt(3)/(4*pi);
c3_conn  = 1 + 3*sqrt(3)/(2*pi);
c3_union = 1.80;

q1_ext = max(1 - p_link, 0);
q2_ext = max(1 - c2_union*p_link, 0);
q3_ext = max(1 - c3_union*p_link, 0);

N1_theory = N * q1_ext^(N - 1);

if N >= 2
    N2_theory = nchoosek(N, 2) ...
        * p_link ...
        * q2_ext^(N - 2);
else
    N2_theory = 0;
end

if N >= 3
    p_conn_3 = c3_conn * p_link^2;
    p_conn_3 = min(max(p_conn_3, 0), 1);

    N3_theory = nchoosek(N, 3) ...
        * p_conn_3 ...
        * q3_ext^(N - 3);
else
    N3_theory = 0;
end

beta0_theory = ...
    2 + N1_theory + N2_theory + N3_theory;

beta0_theory = min(max(beta0_theory, 1), N);

%% ============================================================
%  5. Fraction analytique de liens critiques
%
%  chi_bridge =
%    integral exp[-lambda A_inter(r(alpha))]
%             f_{alpha|link}(alpha) d alpha
%% ============================================================

if alpha_max > 0 && N >= 2

    % Distance corde entre les deux satellites
    r_of_alpha = @(alpha) ...
        2 * R .* sin(alpha ./ 2);

    % Aire d'intersection de deux disques de rayon dmax
    A_inter = @(r) ...
        2 * dmax^2 .* ...
        acos(min(max(r ./ (2*dmax), -1), 1)) ...
        - 0.5 .* r .* ...
        sqrt(max(4*dmax^2 - r.^2, 0));

    % Densité de alpha conditionnellement à l'existence d'un lien
    f_alpha_given_link = @(alpha) ...
        sin(alpha) ./ (1 - cos(alpha_max));

    % Approximation : un lien est critique lorsqu'il ne possède
    % aucun voisin commun
    p_no_common_neighbor = @(alpha) ...
        exp(-lambda .* A_inter(r_of_alpha(alpha)));

    integrand = @(alpha) ...
        p_no_common_neighbor(alpha) ...
        .* f_alpha_given_link(alpha);

    chi_bridge = integral( ...
        integrand, ...
        0, alpha_max, ...
        'RelTol', 1e-8, ...
        'AbsTol', 1e-11);

else
    chi_bridge = 0;
end

chi_bridge = min(max(chi_bridge, 0), 1);

%% ============================================================
%  6. Nombre moyen de liens critiques par composante
%
%  Lbar_C = E[|E|] / E[beta0]
%  Bbar_C = Lbar_C * chi_bridge
%% ============================================================

if beta0_theory > 0
    mean_links_per_component = ...
        E_theory / beta0_theory;
else
    mean_links_per_component = 0;
end

mean_bridges_per_component = ...
    mean_links_per_component * chi_bridge;

%% ============================================================
%  7. Probabilité théorique de rupture d'une composante
%
%  Forme probabiliste :
%
%  p_break_th = 1-(1-p_break_link)^Bbar_C
%% ============================================================

p_break_th = ...
    1 - (1 - p_break_link)^mean_bridges_per_component;

p_break_th = min(max(p_break_th, 0), 1 - eps);

%% Approximation linéaire du passage LaTeX
p_break_th_linear = ...
    p_break_link * mean_bridges_per_component;

p_break_th_linear = ...
    min(max(p_break_th_linear, 0), 1 - eps);

%% ============================================================
%  8. Affichage
%% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' CALCUL THEORIQUE DE p_break - MODELE ALEATOIRE\n');
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

fprintf('------------------------------------------------------------\n');

fprintf('chi_bridge                          : %.8f\n', chi_bridge);
fprintf('Liens moyens par composante         : %.8f\n', ...
    mean_links_per_component);
fprintf('Liens critiques par composante      : %.8f\n', ...
    mean_bridges_per_component);
fprintf('p_break d''un lien                   : %.8f\n', ...
    p_break_link);

fprintf('------------------------------------------------------------\n');

fprintf('p_break théorique probabiliste      : %.8f\n', ...
    p_break_th);
fprintf('p_break théorique linéaire          : %.8f\n', ...
    p_break_th_linear);

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
    'chi_bridge', ...
    'mean_links_per_component', ...
    'mean_bridges_per_component', ...
    'p_break_link', ...
    'p_break_th', ...
    'p_break_th_linear');

fprintf('Résultats sauvegardés dans %s\n', output_file);