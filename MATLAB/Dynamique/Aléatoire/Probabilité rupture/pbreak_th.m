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
%
%  Les composantes isolées, comptées par N1, ne possèdent aucun
%  lien et ne peuvent donc pas être fragmentées par une rupture.
%
%  Lbar_C =
%      E[|E|] / (E[beta0] - E[N1])
%
%  p_bridge_bord représente directement la fraction des liens
%  susceptibles de se rompre qui sont des ponts.
%% ============================================================

n_nonisolated_components = beta0_theory - N1_theory;

if n_nonisolated_components > 0
    mean_links_per_component = ...
        E_theory / n_nonisolated_components;
else
    mean_links_per_component = 0;
end

mean_breaking_bridges_per_component = ...
    mean_links_per_component * p_bridge_bord;

%% ============================================================
%  7. Probabilité théorique de rupture d'une composante
%
%  Le nombre moyen de ruptures critiques par composante non isolée
%  pendant un pas de temps vaut
%
%    mu_break =
%      Lbar_C * p_break_link * p_bridge_bord.
%
%  Dans l'approximation de Poisson :
%
%    p_break_th = 1 - exp(-mu_break).
%% ============================================================

mu_break = ...
    mean_links_per_component ...
    * p_break_link ...
    * p_bridge_bord;

p_break_th = 1 - exp(-mu_break);
p_break_th = min(max(p_break_th, 0), 1 - eps);

%% Approximation linéaire pour mu_break << 1
p_break_th_linear = mu_break;

p_break_th_linear = ...
    min(max(p_break_th_linear, 0), 1 - eps);

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
fprintf('Couronne [D_min,D_max]             : [%.8f, %.8f] km\n', D_min, D_max);
fprintf('p_bridge_bord                       : %.8f\n', p_bridge_bord);
fprintf('Liens moyens par composante         : %.8f\n', ...
    mean_links_per_component);
fprintf('Ponts de bord moyens/composante     : %.8f\n', ...
    mean_breaking_bridges_per_component);
fprintf('p_break d''un lien                  : %.8f\n', ...
    p_break_link);
fprintf('mu_break                            : %.8f\n', ...
    mu_break);

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
    'n_nonisolated_components', ...
    'l_out_eff', 'D_min', 'D_max', ...
    'p_bridge_bord', ...
    'mean_links_per_component', ...
    'mean_breaking_bridges_per_component', ...
    'p_break_link', ...
    'mu_break', ...
    'p_break_th', ...
    'p_break_th_linear');

fprintf('Résultats sauvegardés dans %s\n', output_file);