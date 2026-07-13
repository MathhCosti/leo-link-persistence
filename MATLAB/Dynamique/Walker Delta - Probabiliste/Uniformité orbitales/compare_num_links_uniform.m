function results = compare_num_links_uniform(mat_file)
% Compare le nombre de liens temporel empirique a la theorie uniforme.
% Utilisation :
%   compare_num_links_uniform
%   compare_num_links_uniform('leo_zigzag_analysis_results_delta.mat')
%
% Seules les variables classiques num_edges et time_values sont utilisees.
% Les donnees zigzag ne sont pas utilisees.

close all;

%% Parametres orbitaux aleatoires a inclinaison deterministe
inc_deg = 53;                  % inclinaison commune imposee, en degres
inc = deg2rad(inc_deg);        % radians


if nargin < 1 || isempty(mat_file)
    mat_file = 'leo_zigzag_analysis_results_delta.mat';
end

if ~isfile(mat_file)
    error('Fichier introuvable : %s', mat_file);
end

S = load(mat_file, 'N', 'R', 'dmax', 'lambda', ...
    'time_values', 'num_edges', 'inc_deg');

required = {'N','R','dmax','time_values','num_edges'};
for k = 1:numel(required)
    if ~isfield(S, required{k})
        error('Variable manquante : %s', required{k});
    end
end

N = S.N;
R = S.R;
dmax = S.dmax;
time_values = S.time_values(:);
num_edges = S.num_edges(:);

if numel(time_values) ~= numel(num_edges)
    error('time_values et num_edges doivent avoir la meme longueur.');
end

%% Theorie uniforme exacte sur la sphere
if dmax >= 2*R
    theta_max = pi;
    p_link_uniform = 1;
else
    theta_max = 2*asin(dmax/(2*R));
    p_link_uniform = (1-cos(theta_max))/(2*sin(inc));
end

% Pour une distance euclidienne chordale :
p_link_chord = min(1, dmax^2/(4*R^2));

num_pairs = N*(N-1)/2;
E_uniform = num_pairs * p_link_uniform;
degree_uniform = (N-1)*p_link_uniform;

%% Approximation plane locale
surface_sphere = 4*pi*R^2;
lambda_realized = N/surface_sphere;
degree_plane = lambda_realized*pi*dmax^2;
E_plane = N*degree_plane/2;

%% Statistiques empiriques
E_emp_mean = mean(num_edges);
E_emp_std = std(num_edges);
E_emp_min = min(num_edges);
E_emp_max = max(num_edges);

relative_error = (E_emp_mean-E_uniform)/E_uniform;
relative_time_error = (num_edges-E_uniform)/E_uniform;

%% Figure principale
figure;
hold on;
grid on;

plot(time_values, num_edges, 'LineWidth', 1.3, ...
    'DisplayName', 'Nombre de liens empirique');

yline(E_uniform, '--', 'LineWidth', 2, ...
    'DisplayName', sprintf('Theorie uniforme exacte : %.2f', E_uniform));

yline(E_emp_mean, ':', 'LineWidth', 2, ...
    'DisplayName', sprintf('Moyenne empirique : %.2f', E_emp_mean));

yline(E_plane, '-.', 'LineWidth', 1.5, ...
    'DisplayName', sprintf('Approximation plane : %.2f', E_plane));

xlabel('Temps (s)');
ylabel('Nombre de liens');
title('Nombre de liens : simulation orbitale vs modele uniforme');
legend('Location', 'best');
hold off;

%% Figure de l'ecart relatif
figure;
hold on;
grid on;

plot(time_values, relative_time_error, 'LineWidth', 1.3, ...
    'DisplayName', 'Ecart relatif temporel');

yline(0, '--', 'LineWidth', 1.5, ...
    'DisplayName', 'Accord parfait');

yline(mean(relative_time_error), ':', 'LineWidth', 2, ...
    'DisplayName', sprintf('Ecart moyen : %.2f %%', ...
    100*mean(relative_time_error)));

xlabel('Temps (s)');
ylabel('(E_{emp}(t)-E_{th})/E_{th}');
title('Ecart relatif au modele uniforme');
legend('Location', 'best');
hold off;

%% Console
fprintf('\n=== Comparaison du nombre de liens ===\n');
fprintf('N                                  : %d\n', N);
fprintf('R                                  : %.3f km\n', R);
fprintf('dmax                               : %.3f km\n', dmax);
fprintf('Angle central maximal              : %.6f rad\n', theta_max);
fprintf('p_link uniforme exact              : %.8f\n', p_link_uniform);
fprintf('Verification dmax^2/(4R^2)         : %.8f\n', p_link_chord);
fprintf('Degre moyen theorique uniforme     : %.4f\n', degree_uniform);
fprintf('Nombre de liens theorique uniforme : %.4f\n', E_uniform);
fprintf('Approximation plane                : %.4f\n', E_plane);
fprintf('Moyenne empirique                  : %.4f\n', E_emp_mean);
fprintf('Ecart-type empirique               : %.4f\n', E_emp_std);
fprintf('Minimum / maximum                  : %.0f / %.0f\n', ...
    E_emp_min, E_emp_max);
fprintf('Ecart relatif moyen                : %.2f %%\n', ...
    100*relative_error);

if isfield(S, 'inc_deg')
    fprintf('Inclinaison simulee                : %.2f deg\n', S.inc_deg);
    fprintf(['Attention : la reference theorique suppose une repartition ' ...
        'uniforme sur toute la sphere.\n']);
end

if isfield(S, 'lambda')
    fprintf('lambda impose                      : %.8e km^-2\n', S.lambda);
end
fprintf('lambda realise N/(4*pi*R^2)        : %.8e km^-2\n', ...
    lambda_realized);

%% Structure de sortie
results.N = N;
results.R = R;
results.dmax = dmax;
results.theta_max = theta_max;
results.p_link_uniform = p_link_uniform;
results.degree_uniform = degree_uniform;
results.E_uniform = E_uniform;
results.E_plane = E_plane;
results.E_emp_mean = E_emp_mean;
results.E_emp_std = E_emp_std;
results.E_emp_min = E_emp_min;
results.E_emp_max = E_emp_max;
results.relative_error = relative_error;
results.relative_time_error = relative_time_error;
end
