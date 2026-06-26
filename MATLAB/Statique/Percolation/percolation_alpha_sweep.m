function res = percolation_alpha_sweep(N, alpha_vals, numTests, eta, delta)
    % Calcule les courbes quand N est fixe et alpha_max varie.
    %
    % Bornes renvoyees :
    % - Bound_upper_nonisolated : borne sup par satellites non isoles
    % - Bound_upper_edges       : borne sup par nombre d'aretes
    % - Bound_upper_tree        : borne sup par union bound sur arbres couvrants
    % - Bound_upper_math        : minimum des bornes superieures precedentes
    % - Bound_lower_star        : borne inf conservative par etoile

    P_perc_hat = zeros(size(alpha_vals));
    Hoeffding_low = zeros(size(alpha_vals));
    Hoeffding_up  = zeros(size(alpha_vals));

    Bound_upper_nonisolated = zeros(size(alpha_vals));
    Bound_upper_edges = zeros(size(alpha_vals));
    Bound_upper_tree = zeros(size(alpha_vals));
    Bound_upper_math = zeros(size(alpha_vals));
    Bound_lower_star = zeros(size(alpha_vals));

    eps_H = sqrt(log(2/delta) / (2*numTests));
    m = ceil(eta * N);

    for k = 1:length(alpha_vals)

        alpha_max = alpha_vals(k);
        perc_events = zeros(numTests, 1);

        for t = 1:numTests

            % Satellites uniformes sur la sphere unite.
            % On travaille sur la sphere unite car d_max a deja ete converti en angle.
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

        % Probabilite exacte de lien entre deux points uniformes sur la sphere.
        p_link = (1 - cos(alpha_max)) / 2;

        %% Borne superieure 1 : satellites non isoles
        % Si Cmax >= m, alors il existe au moins m satellites non isoles.
        p_isolated = (1 - p_link)^(N-1);
        E_nonisolated = N * (1 - p_isolated);
        Bound_upper_nonisolated(k) = min(1, E_nonisolated / m);

        %% Borne superieure 2 : nombre d'aretes
        % Si Cmax >= m, alors il existe au moins m-1 aretes.
        if m <= 1
            Bound_upper_edges(k) = 1;
        else
            E_edges = (N * (N-1) / 2) * p_link;
            Bound_upper_edges(k) = min(1, E_edges / (m - 1));
        end

        %% Borne superieure 3 : union bound sur arbres couvrants
        % Une composante de taille m contient un arbre couvrant sur m sommets.
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
        % Si le satellite 1 a au moins m-1 voisins, alors Cmax >= m.
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

    %% Approximation par seuils, sur l'axe alpha_max
    % Seuil de percolation : point où P ~= 1/2.
    % Seuil de connexité   : point où P ~= 1-epsilon.
    deg_perc = 4.512;
    deg_conn = (N-1) * (1 - N^(-1/(N-1)));

    p_perc = deg_perc / (N-1);
    p_conn = deg_conn / (N-1);

    alpha_perc = NaN;
    alpha_conn = NaN;

    if p_perc <= 1
        alpha_perc = acos(1 - 2*p_perc);
    end

    if p_conn <= 1
        alpha_conn = acos(1 - 2*p_conn);
    end

    epsilon_threshold = 0.01;

    if isfinite(alpha_perc) && isfinite(alpha_conn) && abs(alpha_conn-alpha_perc) > eps
        k_sig = log((1-epsilon_threshold)/epsilon_threshold) / abs(alpha_conn-alpha_perc);
        s_sig = sign(alpha_conn-alpha_perc);
        res.P_threshold_approx = 1 ./ (1 + exp(-s_sig*k_sig*(alpha_vals-alpha_perc)));
    else
        res.P_threshold_approx = NaN(size(alpha_vals));
    end

    res.threshold_percolation = alpha_perc;
    res.threshold_connectivity = alpha_conn;
end
