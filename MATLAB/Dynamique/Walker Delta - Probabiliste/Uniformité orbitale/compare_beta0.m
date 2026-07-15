function results = compare_beta0(mat_file)
%COMPARE_BETA0_GEOMETRIC_COMPONENTS
% Compare beta_0(t) empirique (sans zigzag) aux approximations
% geometriques fondees sur les composantes isolees, dimeres et trimeres.
%
% Utilisation :
%   compare_beta0_geometric_components
%   results = compare_beta0_geometric_components( ...
%       'leo_zigzag_analysis_results_delta.mat');
%
% Le modele utilise :
%   p_link = (1-cos(alpha_max))/2
%
%   E[N1] = N (1-p_link)^(N-1)
%
%   c2 = 1 + 3sqrt(3)/(4pi)
%   E[N2] = C(N,2) p_link (1-c2 p_link)^(N-2)
%
%   c3,conn = 1 + 3sqrt(3)/(2pi)
%   c3,union ~= 1.8
%   E[N3] = C(N,3)c3,conn p_link^2
%           (1-c3,union p_link)^(N-3)
%
% Pour convertir ces nombres de petites composantes en beta_0, on ajoute
% une composante residuelle connectee (+1), sauf dans le modele
% entierement connecte beta_0 = 1.

%% Parametres orbitaux aleatoires a inclinaison deterministe
inc_deg = 53;                  % inclinaison commune imposee, en degres
inc = deg2rad(inc_deg);        % radians


    if nargin < 1 || isempty(mat_file)
        mat_file = 'leo_zigzag_analysis_results_delta.mat';
    end

    if ~isfile(mat_file)
        error('Fichier introuvable : %s', mat_file);
    end

    S = load(mat_file, ...
        'N', 'R', 'dmax', 'time_values', 'beta0', 'inc_deg');

    required = {'N','R','dmax','time_values','beta0'};
    for k = 1:numel(required)
        if ~isfield(S, required{k})
            error('Variable manquante dans le fichier MAT : %s', required{k});
        end
    end

    N = S.N;
    R = S.R;
    dmax = S.dmax;
    time_values = S.time_values(:);
    beta0_emp = S.beta0(:);

    if numel(time_values) ~= numel(beta0_emp)
        error('time_values et beta0 doivent avoir la meme longueur.');
    end

    %% Probabilite uniforme de lien
    if dmax >= 2*R
        alpha_max = pi;
        p_link = 1;
    else
        alpha_max = 2*asin(dmax/(2*R));
        p_link = (1-cos(alpha_max))/(2*sin(inc));
    end

    %% Coefficients geometriques
    c2 = 1 + 3*sqrt(3)/(4*pi);
    c3_conn = 1 + 3*sqrt(3)/(2*pi);
    c3_union = 1.8;

    %% Nombres moyens de petites composantes
    EN1 = N*(1-p_link)^(N-1);

    if N >= 2
        base2 = max(0, 1-c2*p_link);
        EN2 = nchoosek(N,2)*p_link*base2^(N-2);
    else
        EN2 = 0;
    end

    if N >= 3
        base3 = max(0, 1-c3_union*p_link);
        EN3 = nchoosek(N,3)*c3_conn*p_link^2*base3^(N-3);
    else
        EN3 = 0;
    end

    %% Approximations de beta_0
    % Une composante residuelle geante est ajoutee dans les trois
    % approximations tronquees.
    beta0_connected = 1;
    beta0_isolated = 1 + EN1;
    beta0_isolated_dimers = 1 + EN1 + EN2;
    beta0_isolated_dimers_trimers = 1 + EN1 + EN2 + EN3;

    theories = [ ...
        beta0_connected, ...
        beta0_isolated, ...
        beta0_isolated_dimers, ...
        beta0_isolated_dimers_trimers];

    theory_names = { ...
        'Connecte', ...
        'Isoles + reste connecte', ...
        'Isoles + dimeres + reste connecte', ...
        'Isoles + dimeres + trimeres + reste connecte'};

    %% Statistiques empiriques
    beta0_mean = mean(beta0_emp);
    beta0_std = std(beta0_emp);
    beta0_min = min(beta0_emp);
    beta0_max = max(beta0_emp);

    abs_errors = abs(theories-beta0_mean);
    rel_errors = abs_errors/max(beta0_mean, eps);
    [~, best_idx] = min(abs_errors);

    %% Figure temporelle
    figure;
    hold on;
    grid on;

    plot(time_values, beta0_emp, 'LineWidth', 1.4, ...
        'DisplayName', '\beta_0 empirique');

    yline(beta0_mean, ':', 'LineWidth', 2, ...
        'DisplayName', sprintf('Moyenne empirique : %.3f', beta0_mean));

    yline(beta0_connected, '-.', 'LineWidth', 1.8, ...
        'DisplayName', sprintf('Connecte : %.3f', beta0_connected));

    yline(beta0_isolated, '--', 'LineWidth', 1.8, ...
        'DisplayName', sprintf('Isoles : %.3f', beta0_isolated));

    yline(beta0_isolated_dimers, '--', 'LineWidth', 1.8, ...
        'DisplayName', sprintf('Isoles + dimeres : %.3f', ...
        beta0_isolated_dimers));

    yline(beta0_isolated_dimers_trimers, '--', 'LineWidth', 2.2, ...
        'DisplayName', sprintf('Isoles + dimeres + trimeres : %.3f', ...
        beta0_isolated_dimers_trimers));

    xlabel('Temps (s)');
    ylabel('\beta_0');
    title('\beta_0 sans zigzag : simulation vs approximations geometriques');
    legend('Location', 'best');
    hold off;

    %% Figure des erreurs relatives
    figure;
    bar(100*rel_errors);
    grid on;
    xticks(1:numel(theory_names));
    xticklabels(theory_names);
    xtickangle(20);
    ylabel('Erreur relative sur la moyenne empirique (%)');
    title('Erreur des approximations geometriques de \beta_0');

    %% Console
    fprintf('\n=== Approximations geometriques de beta_0 ===\n');
    fprintf('N                                      : %d\n', N);
    fprintf('R                                      : %.3f km\n', R);
    fprintf('dmax                                   : %.3f km\n', dmax);
    fprintf('alpha_max                              : %.8f rad\n', alpha_max);
    fprintf('p_link                                 : %.8f\n', p_link);
    fprintf('\n');
    fprintf('c2                                     : %.8f\n', c2);
    fprintf('c3_conn                                : %.8f\n', c3_conn);
    fprintf('c3_union                               : %.8f\n', c3_union);
    fprintf('\n');
    fprintf('E[N1] satellites isoles               : %.8f\n', EN1);
    fprintf('E[N2] dimeres                          : %.8f\n', EN2);
    fprintf('E[N3] trimeres                         : %.8f\n', EN3);
    fprintf('\n');
    fprintf('beta0 connecte                         : %.8f\n', beta0_connected);
    fprintf('beta0 isoles                           : %.8f\n', beta0_isolated);
    fprintf('beta0 isoles + dimeres                 : %.8f\n', ...
        beta0_isolated_dimers);
    fprintf('beta0 isoles + dimeres + trimeres     : %.8f\n', ...
        beta0_isolated_dimers_trimers);
    fprintf('\n');
    fprintf('beta0 empirique moyen                  : %.8f\n', beta0_mean);
    fprintf('Ecart-type empirique                   : %.8f\n', beta0_std);
    fprintf('Minimum / maximum empiriques           : %.0f / %.0f\n', ...
        beta0_min, beta0_max);
    fprintf('Meilleure approximation                : %s\n', ...
        theory_names{best_idx});
    fprintf('Erreur relative correspondante         : %.2f %%\n', ...
        100*rel_errors(best_idx));

    if isfield(S, 'inc_deg')
        fprintf('Inclinaison de la simulation           : %.2f deg\n', ...
            S.inc_deg);
    end

    %% Sortie
    results = struct();
    results.N = N;
    results.R = R;
    results.dmax = dmax;
    results.alpha_max = alpha_max;
    results.p_link = p_link;

    results.c2 = c2;
    results.c3_conn = c3_conn;
    results.c3_union = c3_union;

    results.EN1 = EN1;
    results.EN2 = EN2;
    results.EN3 = EN3;

    results.beta0_connected = beta0_connected;
    results.beta0_isolated = beta0_isolated;
    results.beta0_isolated_dimers = beta0_isolated_dimers;
    results.beta0_isolated_dimers_trimers = ...
        beta0_isolated_dimers_trimers;

    results.beta0_empirical_mean = beta0_mean;
    results.beta0_empirical_std = beta0_std;
    results.beta0_empirical_min = beta0_min;
    results.beta0_empirical_max = beta0_max;

    results.absolute_errors = abs_errors;
    results.relative_errors = rel_errors;
    results.best_theory = theory_names{best_idx};
end
