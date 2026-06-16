function res = percolation_N_sweep(N_vals, alpha_max, numTests, eta, delta)
    % Calcule les courbes quand alpha_max est fixe et N varie.

    P_perc_hat = zeros(size(N_vals));
    Hoeffding_low = zeros(size(N_vals));
    Hoeffding_up  = zeros(size(N_vals));

    Bound_upper_nonisolated = zeros(size(N_vals));
    Bound_upper_edges = zeros(size(N_vals));
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

        % Borne supérieure 1 : satellites non isolés
        p_isolated = (1 - p_link)^(N-1);
        E_nonisolated = N * (1 - p_isolated);
        Bound_upper_nonisolated(k) = min(1, E_nonisolated / m);

        % Borne supérieure 2 : nombre d'arêtes
        E_edges = (N * (N-1) / 2) * p_link;
        Bound_upper_edges(k) = min(1, E_edges / (m - 1));

        Bound_upper_math(k) = min(Bound_upper_nonisolated(k), Bound_upper_edges(k));

        % Borne inférieure conservative par étoile
        Bound_lower_star(k) = 1 - binocdf(m-2, N-1, p_link);
    end

    res.P_perc_hat = P_perc_hat;
    res.Hoeffding_low = Hoeffding_low;
    res.Hoeffding_up = Hoeffding_up;
    res.Bound_upper_math = Bound_upper_math;
    res.Bound_lower_star = Bound_lower_star;
end
