function res = percolation_N_sweep(N_vals, alpha_max, numTests, eta, delta)
    % Calcule les courbes quand alpha_max est fixe et N varie.
    %
    % Bornes renvoyees :
    % - Bound_upper_nonisolated : borne sup par satellites non isoles
    % - Bound_upper_edges       : borne sup par nombre d'aretes
    % - Bound_upper_tree        : borne sup par union bound sur arbres couvrants
    % - Bound_upper_math        : minimum des bornes superieures precedentes
    % - Bound_lower_star        : borne inf conservative par etoile

    P_perc_hat = zeros(size(N_vals));
    Hoeffding_low = zeros(size(N_vals));
    Hoeffding_up  = zeros(size(N_vals));

    Bound_upper_nonisolated = zeros(size(N_vals));
    Bound_upper_edges = zeros(size(N_vals));
    Bound_upper_tree = zeros(size(N_vals));
    Bound_upper_math = zeros(size(N_vals));
    Bound_lower_star = zeros(size(N_vals));

    eps_H = sqrt(log(2/delta) / (2*numTests));
    p_link = (1 - cos(alpha_max)) / 2;

    for k = 1:length(N_vals)

        N = N_vals(k);
        m = ceil(eta * N);

        perc_events = zeros(numTests, 1);

        for t = 1:numTests

            positions = randn(N, 3);
            positions = positions ./ vecnorm(positions, 2, 2);

            cosAlpha = positions * positions';
            cosAlpha = max(min(cosAlpha, 1), -1);

            alpha = acos(cosAlpha);

            A = (alpha <= alpha_max);
            A = A & ~eye(N);

            G = graph(A);

            comp = conncomp(G);
            comp_sizes = accumarray(comp', 1);
            Cmax = max(comp_sizes);

            perc_events(t) = (Cmax >= m);
        end

        P_perc_hat(k) = mean(perc_events);

        Hoeffding_low(k) = max(0, P_perc_hat(k) - eps_H);
        Hoeffding_up(k)  = min(1, P_perc_hat(k) + eps_H);

        %% Borne superieure 1 : satellites non isoles
        p_isolated = (1 - p_link)^(N-1);
        E_nonisolated = N * (1 - p_isolated);
        Bound_upper_nonisolated(k) = min(1, E_nonisolated / m);

        %% Borne superieure 2 : nombre d'aretes
        if m <= 1
            Bound_upper_edges(k) = 1;
        else
            E_edges = (N * (N-1) / 2) * p_link;
            Bound_upper_edges(k) = min(1, E_edges / (m - 1));
        end

        %% Borne superieure 3 : union bound sur arbres couvrants
        if p_link == 0
            Bound_upper_tree(k) = 0;
        elseif m <= 1
            Bound_upper_tree(k) = 1;
        else
            logBound = gammaln(N+1) ...
                - gammaln(m+1) ...
                - gammaln(N-m+1) ...
                + (m-2)*log(m) ...
                + (m-1)*log(p_link);

            Bound_upper_tree(k) = min(1, exp(logBound));
        end

        %% Meilleure borne superieure disponible
        Bound_upper_math(k) = min([Bound_upper_nonisolated(k), ...
                                   Bound_upper_edges(k), ...
                                   Bound_upper_tree(k)]);

        %% Borne inferieure conservative par etoile
        Bound_lower_star(k) = 1 - binocdf(m-2, N-1, p_link);
    end

    res.P_perc_hat = P_perc_hat;
    res.Hoeffding_low = Hoeffding_low;
    res.Hoeffding_up = Hoeffding_up;

    res.Bound_upper_nonisolated = Bound_upper_nonisolated;
    res.Bound_upper_edges = Bound_upper_edges;
    res.Bound_upper_tree = Bound_upper_tree;
    res.Bound_upper_math = Bound_upper_math;
    res.Bound_lower_star = Bound_lower_star;

    %% Approximation par seuils, sur l'axe N
    % Seuil de percolation : point où P ~= 1/2.
    % Seuil de connexité   : point où P ~= 1-epsilon.
    deg_perc = 4.512;
    N_perc = 1 + deg_perc / p_link;

    f_conn = @(x) 1 - x.^(-1./(x-1)) - p_link;
    N_conn = NaN;

    N_min_search = max(2.0001, min(N_vals)/10);
    N_max_search = max(N_vals)*10;

    try
        if f_conn(N_min_search) * f_conn(N_max_search) <= 0
            N_conn = fzero(f_conn, [N_min_search, N_max_search]);
        end
    catch
        N_conn = NaN;
    end

    epsilon_threshold = 0.01;

    if isfinite(N_perc) && isfinite(N_conn) && abs(N_conn-N_perc) > eps
        k_sig = log((1-epsilon_threshold)/epsilon_threshold) / abs(N_conn-N_perc);
        s_sig = sign(N_conn-N_perc);
        res.P_threshold_approx = 1 ./ (1 + exp(-s_sig*k_sig*(N_vals-N_perc)));
    else
        res.P_threshold_approx = NaN(size(N_vals));
    end

    res.threshold_percolation = N_perc;
    res.threshold_connectivity = N_conn;
end
