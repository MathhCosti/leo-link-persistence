clear; clc; close all;

%% ============================================================
%  CALCUL THEORIQUE DE p_merge
%  Modele aleatoire a vecteurs tangentiels
%
%  Ce script reprend uniquement la partie theorique de p_merge
%  presente dans barcodes.m.
%
%  Entree :
%    analysis_temp_results.mat
%
%  Sortie :
%    pmerge_th_random_results.mat
%% ============================================================

%% Chargement des parametres du modele
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

%% Parametres orbitaux et vitesse relative moyenne
mu = 398600;                    % km^3/s^2
omega = sqrt(mu / R^3);         % rad/s
v_orb = R * omega;              % km/s
v_rel = (4/pi) * v_orb;         % km/s

%% Probabilite de fusion brute
A_sweep = 2 * dmax * v_rel * dt;

p_merge_raw = 1 - exp(-lambda * A_sweep);
p_merge_raw = min(max(p_merge_raw, 0), 1 - eps);

%% Probabilite de lien et nombre moyen d'aretes
alpha_max = 2 * asin(min(dmax / (2*R), 1));
p_link = (1 - cos(alpha_max)) / 2;
E_theory = N * (N - 1) / 2 * p_link;

%% Approximation theorique de beta0
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
    p_conn_3 = min(max(c3_conn * p_link^2, 0), 1);
    N3_theory = nchoosek(N, 3) ...
        * p_conn_3 ...
        * q3_ext^(N - 3);
else
    N3_theory = 0;
end

% Formule reprise telle quelle de barcodes.m
beta0_theory = ...
    1 + N1_theory + N2_theory + N3_theory;

beta0_theory = min(max(beta0_theory, 1), N);

%% Facteur topologique de fusion
if E_theory > 0
    chi_merge = (N - beta0_theory) / E_theory;
else
    chi_merge = 0;
end

chi_merge = min(max(chi_merge, 0), 1);

%% Probabilite theorique de fusion corrigee
A_sweep_corrected = ...
    2 * dmax * v_rel * dt * chi_merge;

p_merge_th = ...
    1 - exp(-lambda * A_sweep_corrected);

p_merge_th = min(max(p_merge_th, 0), 1 - eps);

%% Approximation lineaire
p_merge_th_linear = ...
    lambda * A_sweep_corrected;

p_merge_th_linear = ...
    min(max(p_merge_th_linear, 0), 1 - eps);

%% Affichage
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
fprintf('alpha_max                           : %.8f rad\n', alpha_max);
fprintf('p_link                              : %.8e\n', p_link);
fprintf('E[|E|]                              : %.8f\n', E_theory);
fprintf('------------------------------------------------------------\n');
fprintf('E[N1]                               : %.8f\n', N1_theory);
fprintf('E[N2]                               : %.8f\n', N2_theory);
fprintf('E[N3]                               : %.8f\n', N3_theory);
fprintf('E[beta0]                            : %.8f\n', beta0_theory);
fprintf('------------------------------------------------------------\n');
fprintf('Aire balayee brute                  : %.8f km^2\n', A_sweep);
fprintf('p_merge brute                       : %.8f\n', p_merge_raw);
fprintf('chi_merge                           : %.8f\n', chi_merge);
fprintf('Aire balayee corrigee               : %.8f km^2\n', A_sweep_corrected);
fprintf('------------------------------------------------------------\n');
fprintf('p_merge theorique probabiliste      : %.8f\n', p_merge_th);
fprintf('p_merge theorique lineaire          : %.8f\n', p_merge_th_linear);
fprintf('============================================================\n');

%% Sauvegarde
output_file = 'pmerge_th_results.mat';

save(output_file, ...
    'N', 'R', 'lambda', 'dmax', 'dt', ...
    'mu', 'omega', 'v_orb', 'v_rel', ...
    'alpha_max', 'p_link', 'E_theory', ...
    'c2_union', 'c3_conn', 'c3_union', ...
    'N1_theory', 'N2_theory', 'N3_theory', ...
    'beta0_theory', ...
    'A_sweep', 'p_merge_raw', ...
    'chi_merge', 'A_sweep_corrected', ...
    'p_merge_th', 'p_merge_th_linear');

fprintf('Resultats sauvegardes dans %s\n', output_file);
