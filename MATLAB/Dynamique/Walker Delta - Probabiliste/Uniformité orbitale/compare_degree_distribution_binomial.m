function results = compare_degree_distribution_binomial(mat_file)
%COMPARE_DEGREE_DISTRIBUTION_BINOMIAL
% Compare la distribution empirique des degres, moyennee dans le temps,
% a une loi binomiale Bin(N-1,p).
%
% Utilisation :
%   compare_degree_distribution_binomial
%   results = compare_degree_distribution_binomial( ...
%       'leo_zigzag_analysis_results_delta.mat');
%
% Le fichier MAT doit contenir :
%   N, R, dmax, time_values, Adjacency, num_edges
%
% Deux probabilites sont comparees :
%   1) p_link empirique = 2*mean(num_edges)/(N*(N-1))
%   2) p_link uniforme sur la sphere
%      = (1-cos(alpha_max))/2

    if nargin < 1 || isempty(mat_file)
        mat_file = 'leo_zigzag_analysis_results_delta.mat';
    end

    if ~isfile(mat_file)
        error('Fichier introuvable : %s', mat_file);
    end

    S = load(mat_file, ...
        'N', 'R', 'dmax', 'time_values', 'Adjacency', ...
        'num_edges', 'inc_deg');

    required = {'N','R','dmax','time_values','Adjacency','num_edges'};
    for k = 1:numel(required)
        if ~isfield(S, required{k})
            error('Variable manquante dans le fichier MAT : %s', required{k});
        end
    end

    N = S.N;
    R = S.R;
    dmax = S.dmax;
    time_values = S.time_values(:);
    Adjacency = S.Adjacency;
    num_edges = S.num_edges(:);

    Nt = numel(time_values);

    if numel(Adjacency) ~= Nt
        error('Le nombre de matrices Adjacency ne correspond pas a time_values.');
    end

    if numel(num_edges) ~= Nt
        error('Le nombre de valeurs num_edges ne correspond pas a time_values.');
    end

    %% ========================================================
    %  Distribution empirique temporelle des degres
    %% ========================================================

    degree_counts_time = zeros(Nt, N);   % colonnes k=0,...,N-1
    degree_values_all = zeros(Nt*N,1);

    cursor = 1;

    for t_idx = 1:Nt
        A = Adjacency{t_idx};

        deg = full(sum(A,2));

        counts = accumarray(deg+1, 1, [N,1]);
        degree_counts_time(t_idx,:) = counts.';

        degree_values_all(cursor:cursor+N-1) = deg;
        cursor = cursor + N;
    end

    mean_degree_count = mean(degree_counts_time,1);
    empirical_degree_pmf = mean_degree_count / N;

    degree_axis = 0:N-1;

    empirical_mean_degree = mean(degree_values_all);
    empirical_var_degree = var(degree_values_all,1);

    empirical_isolated_prob = empirical_degree_pmf(1);
    empirical_isolated_count = N*empirical_isolated_prob;

    %% ========================================================
    %  Probabilites de lien
    %% ========================================================

    p_emp = 2*mean(num_edges)/(N*(N-1));

    if dmax >= 2*R
        alpha_max = pi;
        p_uniform = 1;
    else
        alpha_max = 2*asin(dmax/(2*R));
        p_uniform = (1-cos(alpha_max))/2;
    end

    %% ========================================================
    %  Lois binomiales
    %% ========================================================

    pmf_bin_emp = binopdf(degree_axis, N-1, p_emp);
    pmf_bin_uniform = binopdf(degree_axis, N-1, p_uniform);

    mean_bin_emp = (N-1)*p_emp;
    var_bin_emp = (N-1)*p_emp*(1-p_emp);

    isolated_bin_emp = N*(1-p_emp)^(N-1);
    isolated_bin_uniform = N*(1-p_uniform)^(N-1);

    %% ========================================================
    %  Distance entre distributions
    %% ========================================================

    tv_emp = 0.5*sum(abs(empirical_degree_pmf-pmf_bin_emp));
    tv_uniform = 0.5*sum(abs(empirical_degree_pmf-pmf_bin_uniform));

    rmse_emp = sqrt(mean((empirical_degree_pmf-pmf_bin_emp).^2));
    rmse_uniform = sqrt(mean((empirical_degree_pmf-pmf_bin_uniform).^2));

    overdispersion_ratio = empirical_var_degree/max(var_bin_emp,eps);

    %% ========================================================
    %  Limite d'affichage
    %% ========================================================

    max_degree_observed = find(mean_degree_count>0,1,'last')-1;
    max_degree_plot = min(N-1, max(15, max_degree_observed+3));
    idx_plot = 1:(max_degree_plot+1);

    %% ========================================================
    %  Figure principale : PMF
    %% ========================================================

    figure;
    hold on;
    grid on;

    bar(degree_axis(idx_plot), empirical_degree_pmf(idx_plot), 1, ...
        'FaceAlpha', 0.45, ...
        'DisplayName', 'Distribution empirique');

    plot(degree_axis(idx_plot), pmf_bin_emp(idx_plot), 'o-', ...
        'LineWidth', 1.8, ...
        'DisplayName', sprintf('Binomiale avec p_{emp}=%.5f', p_emp));

    plot(degree_axis(idx_plot), pmf_bin_uniform(idx_plot), '--', ...
        'LineWidth', 1.8, ...
        'DisplayName', sprintf('Binomiale uniforme p=%.5f', p_uniform));

    xlabel('Degre k');
    ylabel('Probabilite P(D=k)');
    title('Distribution empirique des degres vs loi binomiale');
    legend('Location','best');
    hold off;

    %% Figure semi-log pour les queues
    positive_idx = idx_plot(empirical_degree_pmf(idx_plot)>0 | ...
        pmf_bin_emp(idx_plot)>0 | pmf_bin_uniform(idx_plot)>0);

    figure;
    hold on;
    grid on;

    semilogy(degree_axis(positive_idx), ...
        max(empirical_degree_pmf(positive_idx), eps), 'o-', ...
        'LineWidth', 1.5, ...
        'DisplayName', 'Empirique');

    semilogy(degree_axis(positive_idx), ...
        max(pmf_bin_emp(positive_idx), eps), '--', ...
        'LineWidth', 1.8, ...
        'DisplayName', 'Binomiale avec p_{emp}');

    semilogy(degree_axis(positive_idx), ...
        max(pmf_bin_uniform(positive_idx), eps), ':', ...
        'LineWidth', 1.8, ...
        'DisplayName', 'Binomiale uniforme');

    xlabel('Degre k');
    ylabel('Probabilite (echelle logarithmique)');
    title('Queues de distribution des degres');
    legend('Location','best');
    hold off;

    %% Figure temporelle : nombre de satellites isoles
    isolated_time = degree_counts_time(:,1);

    figure;
    hold on;
    grid on;

    plot(time_values, isolated_time, 'LineWidth', 1.3, ...
        'DisplayName', 'Nombre empirique de satellites isoles');

    yline(isolated_bin_emp, '--', 'LineWidth', 2, ...
        'DisplayName', sprintf('Binomiale avec p_{emp}: %.2f', ...
        isolated_bin_emp));

    yline(isolated_bin_uniform, ':', 'LineWidth', 2, ...
        'DisplayName', sprintf('Binomiale uniforme: %.2f', ...
        isolated_bin_uniform));

    yline(mean(isolated_time), '-.', 'LineWidth', 2, ...
        'DisplayName', sprintf('Moyenne empirique: %.2f', ...
        mean(isolated_time)));

    xlabel('Temps (s)');
    ylabel('Nombre de satellites de degre nul');
    title('Satellites isoles : empirique vs loi binomiale');
    legend('Location','best');
    hold off;

    %% ========================================================
    %  Console
    %% ========================================================

    fprintf('\n=== Distribution des degres ===\n');
    fprintf('N                                  : %d\n', N);
    fprintf('Nombre d''instants                 : %d\n', Nt);
    fprintf('p_link empirique                   : %.8f\n', p_emp);
    fprintf('p_link uniforme                    : %.8f\n', p_uniform);
    fprintf('\n');
    fprintf('Degre moyen empirique              : %.6f\n', ...
        empirical_mean_degree);
    fprintf('Degre moyen binomial               : %.6f\n', ...
        mean_bin_emp);
    fprintf('Variance empirique des degres      : %.6f\n', ...
        empirical_var_degree);
    fprintf('Variance binomiale                 : %.6f\n', ...
        var_bin_emp);
    fprintf('Rapport de surdispersion           : %.6f\n', ...
        overdispersion_ratio);
    fprintf('\n');
    fprintf('N1 empirique moyen                 : %.6f\n', ...
        empirical_isolated_count);
    fprintf('N1 binomial avec p_emp             : %.6f\n', ...
        isolated_bin_emp);
    fprintf('N1 binomial uniforme               : %.6f\n', ...
        isolated_bin_uniform);
    fprintf('\n');
    fprintf('Distance variation totale p_emp    : %.6f\n', tv_emp);
    fprintf('Distance variation totale uniforme : %.6f\n', tv_uniform);
    fprintf('RMSE p_emp                         : %.6e\n', rmse_emp);
    fprintf('RMSE uniforme                      : %.6e\n', rmse_uniform);

    if overdispersion_ratio > 1
        fprintf(['Conclusion : la distribution empirique est surdispersee ' ...
            'par rapport a la loi binomiale.\n']);
    else
        fprintf(['Conclusion : la distribution empirique n''est pas ' ...
            'surdispersee par rapport a la loi binomiale.\n']);
    end

    %% ========================================================
    %  Sauvegarde
    %% ========================================================

    save('degree_distribution_comparison.mat', ...
        'degree_axis', 'degree_counts_time', ...
        'mean_degree_count', 'empirical_degree_pmf', ...
        'pmf_bin_emp', 'pmf_bin_uniform', ...
        'p_emp', 'p_uniform', 'alpha_max', ...
        'empirical_mean_degree', 'empirical_var_degree', ...
        'mean_bin_emp', 'var_bin_emp', ...
        'overdispersion_ratio', ...
        'empirical_isolated_count', ...
        'isolated_bin_emp', 'isolated_bin_uniform', ...
        'tv_emp', 'tv_uniform', 'rmse_emp', 'rmse_uniform');

    %% Structure de sortie
    results = struct();
    results.degree_axis = degree_axis;
    results.empirical_degree_pmf = empirical_degree_pmf;
    results.pmf_bin_emp = pmf_bin_emp;
    results.pmf_bin_uniform = pmf_bin_uniform;
    results.p_emp = p_emp;
    results.p_uniform = p_uniform;
    results.empirical_mean_degree = empirical_mean_degree;
    results.empirical_var_degree = empirical_var_degree;
    results.binomial_mean_degree = mean_bin_emp;
    results.binomial_var_degree = var_bin_emp;
    results.overdispersion_ratio = overdispersion_ratio;
    results.empirical_isolated_count = empirical_isolated_count;
    results.isolated_bin_emp = isolated_bin_emp;
    results.isolated_bin_uniform = isolated_bin_uniform;
    results.tv_emp = tv_emp;
    results.tv_uniform = tv_uniform;
end
