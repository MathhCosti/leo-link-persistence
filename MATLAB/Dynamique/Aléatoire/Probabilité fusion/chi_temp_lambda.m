clear; clc; close all;

%% Paramètres
Re = 6371;            % rayon Terre [km]
h = 550;              % altitude [km]
R = Re + h;           % rayon orbital [km]

d_max = 1500;                 % distance max [km]
alpha_max = 20*pi/180;        % angle max [rad]

% Densité : satellites par 10^6 km^2
lambda_scaled_values = linspace(0.01, 2, 20);

n_iter = 3;         % nombre de réalisations Monte-Carlo
rng(1);

%% Stockage Monte-Carlo
Betti0_all = zeros(n_iter, length(lambda_scaled_values));
Betti1_graph_all = zeros(n_iter, length(lambda_scaled_values));
Betti1_complex_all = zeros(n_iter, length(lambda_scaled_values));

% Nombre d'arêtes empiriques
E_all = zeros(n_iter, length(lambda_scaled_values));

% Fraction corrective empirique : (N - beta0) / E
Frac_bridge_all = zeros(n_iter, length(lambda_scaled_values));

N_values = zeros(size(lambda_scaled_values));

%% Stockage des approximations théoriques
Beta0_theory_sparse = zeros(size(lambda_scaled_values));
Beta0_theory_isolated = zeros(size(lambda_scaled_values));
Beta0_theory_dimers = zeros(size(lambda_scaled_values));
Beta0_theory_trimers = zeros(size(lambda_scaled_values));

Beta1_graph_theory_sparse = zeros(size(lambda_scaled_values));
Beta1_graph_theory_isolated = zeros(size(lambda_scaled_values));
Beta1_graph_theory_dimers = zeros(size(lambda_scaled_values));
Beta1_graph_theory_trimers = zeros(size(lambda_scaled_values));

% Nombre d'arêtes théorique
E_theory_values = zeros(size(lambda_scaled_values));

% Fractions correctives théoriques
Frac_bridge_theory_sparse = zeros(size(lambda_scaled_values));
Frac_bridge_theory_isolated = zeros(size(lambda_scaled_values));
Frac_bridge_theory_dimers = zeros(size(lambda_scaled_values));
Frac_bridge_theory_trimers = zeros(size(lambda_scaled_values));

%% Paramètre effectif de connexion
alpha_from_dmax = 2*asin(min(d_max/(2*R), 1));
alpha_eff = min(alpha_max, alpha_from_dmax);
p_link = (1 - cos(alpha_eff))/2;

% Constantes géométriques utilisées aussi dans betti_dmax.m
c2_union = 1.4135;
c3_conn = 1.827;
c3_union = 1.80;

%% Théorie et valeurs de N : ne dépendent pas de la réalisation aléatoire
for k = 1:length(lambda_scaled_values)

    lambda_scaled = lambda_scaled_values(k);

    % Conversion en satellites / km^2
    lambda = lambda_scaled / 1e6;

    % Nombre moyen de satellites sur la sphère
    N = round(lambda * 4*pi*R^2);
    N = max(N, 1);

    N_values(k) = N;

    % Espérance du nombre d'arêtes
    E_theory = N*(N-1)/2 * p_link;
    E_theory_values(k) = E_theory;

    % ========================================================
    % Théorie géométrique utilisée dans betti_dmax.m
    % ========================================================

    % Composantes isolées (taille 1)
    C1_theory = N * (1 - p_link)^(N-1);

    % Dimères géométriques :
    % deux sommets reliés et aucun sommet extérieur dans
    % l'union de leurs voisinages.
    if N >= 2
        C2_theory = nchoosek(N,2) * p_link ...
            * max(1 - c2_union*p_link, 0)^(N-2);
    else
        C2_theory = 0;
    end

    % Trimères géométriques :
    % facteur de connexité interne c3_conn et facteur d'aire
    % d'union extérieure c3_union.
    if N >= 3
        C3_theory = nchoosek(N,3) * c3_conn * p_link^2 ...
            * max(1 - c3_union*p_link, 0)^(N-3);
    else
        C3_theory = 0;
    end

    % Modèle "connectés" / peu dense : beta0 approx N - E
    beta0_sparse = max(N - E_theory, 1);

    % Même convention que dans betti_dmax.m :
    % deux composantes résiduelles + petites composantes.
    beta0_isolated = min(max(2 + C1_theory, 1), N);
    beta0_dimers = min(max(2 + C1_theory + C2_theory, 1), N);
    beta0_trimers = min(max(2 + C1_theory + C2_theory + C3_theory, 1), N);

    Beta0_theory_sparse(k) = beta0_sparse;
    Beta0_theory_isolated(k) = beta0_isolated;
    Beta0_theory_dimers(k) = beta0_dimers;
    Beta0_theory_trimers(k) = beta0_trimers;

    Beta1_graph_theory_sparse(k) = max(E_theory - N + beta0_sparse, 0);
    Beta1_graph_theory_isolated(k) = max(E_theory - N + beta0_isolated, 0);
    Beta1_graph_theory_dimers(k) = max(E_theory - N + beta0_dimers, 0);
    Beta1_graph_theory_trimers(k) = max(E_theory - N + beta0_trimers, 0);

    % Fraction corrective théorique : (N - beta0) / E
    if E_theory > 0
        Frac_bridge_theory_sparse(k) = max(0, min(1, ...
            (N - beta0_sparse) / E_theory));
        Frac_bridge_theory_isolated(k) = max(0, min(1, ...
            (N - beta0_isolated) / E_theory));
        Frac_bridge_theory_dimers(k) = max(0, min(1, ...
            (N - beta0_dimers) / E_theory));
        Frac_bridge_theory_trimers(k) = max(0, min(1, ...
            (N - beta0_trimers) / E_theory));
    else
        Frac_bridge_theory_sparse(k) = NaN;
        Frac_bridge_theory_isolated(k) = NaN;
        Frac_bridge_theory_dimers(k) = NaN;
        Frac_bridge_theory_trimers(k) = NaN;
    end
end

%% Boucle Monte-Carlo
for it = 1:n_iter

    %% Boucle sur lambda
    for k = 1:length(lambda_scaled_values)

        N = N_values(k);

        %% Génération uniforme des satellites sur la sphère
        u = rand(N,1);
        v = rand(N,1);

        theta = 2*pi*u;
        phi = acos(2*v - 1);

        x = R * sin(phi).*cos(theta);
        y = R * sin(phi).*sin(theta);
        z = R * cos(phi);

        P = [x y z];

        %% Matrices de distance et d'angle
        if N >= 2
            D = squareform(pdist(P));
        else
            D = 0;
        end

        U = P ./ vecnorm(P,2,2);
        CosAlpha = U * U.';
        CosAlpha = max(min(CosAlpha,1),-1);
        Alpha = acos(CosAlpha);

        %% Graphe de liens
        A = (D <= d_max) & (Alpha <= alpha_max);
        A(1:N+1:end) = false;
        A = A | A.';

        %% Betti simulés + nombre d'arêtes
        [b0, b1_graph, b1_complex, E] = compute_betti_0_1(A);

        Betti0_all(it,k) = b0;
        Betti1_graph_all(it,k) = b1_graph;
        Betti1_complex_all(it,k) = b1_complex;
        E_all(it,k) = E;

        % Fraction corrective empirique : (N - beta0) / E
        % Si E = 0, on met NaN car il n'y a aucun lien à corriger.
        if E > 0
            Frac_bridge_all(it,k) = (N - b0) / E;
        else
            Frac_bridge_all(it,k) = NaN;
        end
    end

    fprintf('Itération %d / %d terminée\n', it, n_iter);
end

%% Moyennes et écarts-types Monte-Carlo
Betti0 = mean(Betti0_all, 1);
Betti1_graph = mean(Betti1_graph_all, 1);
Betti1_complex = mean(Betti1_complex_all, 1);

Betti0_std = std(Betti0_all, 0, 1);
Betti1_graph_std = std(Betti1_graph_all, 0, 1);
Betti1_complex_std = std(Betti1_complex_all, 0, 1);

E_mean = mean(E_all, 1);
E_std = std(E_all, 0, 1);

% Moyenne du ratio réalisation par réalisation
Frac_bridge = mean(Frac_bridge_all, 1, 'omitnan');
Frac_bridge_std = std(Frac_bridge_all, 0, 1, 'omitnan');

% Ratio des moyennes : (N - E[beta0]) / E[E]
Frac_bridge_ratio_means = (N_values - Betti0) ./ E_mean;
Frac_bridge_ratio_means(E_mean == 0) = NaN;
Frac_bridge_ratio_means = max(0, min(1, Frac_bridge_ratio_means));

%% Figure 1 : beta0 simulation vs théorie
figure;
errorbar(lambda_scaled_values, Betti0, Betti0_std, 'LineWidth', 1.5); hold on;
plot(lambda_scaled_values, Beta0_theory_sparse, '--', 'LineWidth', 2);
plot(lambda_scaled_values, Beta0_theory_isolated, ':', 'LineWidth', 2);
plot(lambda_scaled_values, Beta0_theory_dimers, '-.', 'LineWidth', 2);
plot(lambda_scaled_values, Beta0_theory_trimers, '-', 'LineWidth', 2);
grid on;
xlabel('\lambda [satellites / 10^6 km^2]');
ylabel('\beta_0');
legend('Simulation moyenne \pm écart-type', ...
       'Théorie connectés', ...
       'Théorie isolés', ...
       'Théorie dimères géométriques', ...
       'Théorie trimères géométriques', ...
       'Location', 'best');
title(sprintf('\\beta_0 moyen en fonction de \\lambda — %d itérations', n_iter));

%% Figure 2 : beta1 graphe seul simulation vs théorie
figure;
errorbar(lambda_scaled_values, Betti1_graph, Betti1_graph_std, 'LineWidth', 1.5); hold on;
plot(lambda_scaled_values, Beta1_graph_theory_sparse, '--', 'LineWidth', 2);
plot(lambda_scaled_values, Beta1_graph_theory_isolated, ':', 'LineWidth', 2);
plot(lambda_scaled_values, Beta1_graph_theory_dimers, '-.', 'LineWidth', 2);
plot(lambda_scaled_values, Beta1_graph_theory_trimers, '-', 'LineWidth', 2);
grid on;
xlabel('\lambda [satellites / 10^6 km^2]');
ylabel('\beta_1^{graphe}');
legend('Simulation moyenne \pm écart-type', ...
       'Théorie connectés', ...
       'Théorie isolés', ...
       'Théorie dimères géométriques', ...
       'Théorie trimères géométriques', ...
       'Location', 'best');
title(sprintf('\\beta_1 du graphe moyen en fonction de \\lambda — %d itérations', n_iter));

%% Figure 3 : nombre moyen d'arêtes
figure;
errorbar(lambda_scaled_values, E_mean, E_std, 'LineWidth', 1.5); hold on;
plot(lambda_scaled_values, E_theory_values, '--', 'LineWidth', 2);
grid on;
xlabel('\lambda [satellites / 10^6 km^2]');
ylabel('|E|');
legend('Simulation moyenne \pm écart-type', ...
       'Théorie', ...
       'Location', 'best');
title(sprintf('Nombre moyen de liens |E| en fonction de \\lambda — %d itérations', n_iter));

%% Figure 4 : fraction corrective (N - beta0) / |E|
figure;
errorbar(lambda_scaled_values, Frac_bridge, Frac_bridge_std, 'LineWidth', 1.5); hold on;
plot(lambda_scaled_values, Frac_bridge_ratio_means, '-.', 'LineWidth', 2);
plot(lambda_scaled_values, Frac_bridge_theory_sparse, '--', 'LineWidth', 2);
plot(lambda_scaled_values, Frac_bridge_theory_isolated, ':', 'LineWidth', 2);
plot(lambda_scaled_values, Frac_bridge_theory_dimers, '-.', 'LineWidth', 2);
plot(lambda_scaled_values, Frac_bridge_theory_trimers, '-', 'LineWidth', 2);
grid on;
xlabel('\lambda [satellites / 10^6 km^2]');
ylabel('(N - \beta_0)/|E|');
ylim([0 1.05]);
legend('Simulation : moyenne de (N-\beta_0)/|E|', ...
       'Simulation : (N-\langle\beta_0\rangle)/\langle|E|\rangle', ...
       'Théorie connectés', ...
       'Théorie isolés', ...
       'Théorie dimères géométriques', ...
       'Théorie trimères géométriques', ...
       'Location', 'best');
title(sprintf('Fraction corrective moyenne en fonction de \\lambda — %d itérations', n_iter));

%% Affichage console
T = table(lambda_scaled_values(:), N_values(:), Betti0(:), E_mean(:), ...
          Beta0_theory_isolated(:), Beta0_theory_dimers(:), ...
          Beta0_theory_trimers(:), ...
          Frac_bridge(:), Frac_bridge_ratio_means(:), ...
          Frac_bridge_theory_sparse(:), Frac_bridge_theory_isolated(:), ...
          Frac_bridge_theory_dimers(:), Frac_bridge_theory_trimers(:), ...
          'VariableNames', {'lambda_scaled','N','beta0_emp','E_emp', ...
                            'beta0_theory_isoles', ...
                            'beta0_theory_dimeres', ...
                            'beta0_theory_trimeres', ...
                            'frac_emp_mean_ratio', ...
                            'frac_emp_ratio_means', ...
                            'frac_theory_connectes', ...
                            'frac_theory_isoles', ...
                            'frac_theory_dimeres', ...
                            'frac_theory_trimeres'});
disp(T);

%% Sauvegarde
save('chi_temp_lambda_geom_results.mat', ...
    'lambda_scaled_values', 'N_values', ...
    'Betti0', 'Betti0_std', ...
    'Betti1_graph', 'Betti1_graph_std', ...
    'E_mean', 'E_std', ...
    'Beta0_theory_sparse', 'Beta0_theory_isolated', ...
    'Beta0_theory_dimers', 'Beta0_theory_trimers', ...
    'Beta1_graph_theory_sparse', ...
    'Beta1_graph_theory_isolated', ...
    'Beta1_graph_theory_dimers', ...
    'Beta1_graph_theory_trimers', ...
    'Frac_bridge', 'Frac_bridge_std', ...
    'Frac_bridge_ratio_means', ...
    'Frac_bridge_theory_sparse', ...
    'Frac_bridge_theory_isolated', ...
    'Frac_bridge_theory_dimers', ...
    'Frac_bridge_theory_trimers', ...
    'p_link', 'c2_union', 'c3_conn', 'c3_union');

fprintf('Résultats sauvegardés dans chi_temp_lambda_geom_results.mat\n');


%% ============================================================
%% Fonction Betti
%% ============================================================

function [beta0, beta1_graph, beta1_complex, E] = compute_betti_0_1(A)

    N = size(A,1);

    G = graph(A);
    comp = conncomp(G);
    beta0 = max(comp);

    [I,J] = find(triu(A,1));
    edges = [I J];
    E = size(edges,1);

    beta1_graph = E - N + beta0;

    triangles = [];

    for i = 1:N
        neigh = find(A(i,:) & (1:N > i));

        for a = 1:length(neigh)
            j = neigh(a);

            for b = a+1:length(neigh)
                k = neigh(b);

                if A(j,k)
                    triangles = [triangles; i j k]; %#ok<AGROW>
                end
            end
        end
    end

    T = size(triangles,1);

    if T == 0
        beta1_complex = beta1_graph;
        return;
    end

    edge_map = containers.Map;

    for e = 1:E
        edge_map(edge_key(edges(e,1), edges(e,2))) = e;
    end

    B2 = false(E,T);

    for t = 1:T
        tri = triangles(t,:);

        e1 = edge_map(edge_key(tri(1), tri(2)));
        e2 = edge_map(edge_key(tri(1), tri(3)));
        e3 = edge_map(edge_key(tri(2), tri(3)));

        B2(e1,t) = true;
        B2(e2,t) = true;
        B2(e3,t) = true;
    end

    rankB2 = rank_mod2(B2);

    beta1_complex = beta1_graph - rankB2;
    beta1_complex = max(beta1_complex, 0);
end


function r = rank_mod2(M)

    M = logical(M);
    [m,n] = size(M);

    r = 0;
    row = 1;

    for col = 1:n

        if row > m
            break;
        end

        pivot = find(M(row:m,col), 1);

        if isempty(pivot)
            continue;
        end

        pivot = pivot + row - 1;

        temp = M(row,:);
        M(row,:) = M(pivot,:);
        M(pivot,:) = temp;

        for i = 1:m
            if i ~= row && M(i,col)
                M(i,:) = xor(M(i,:), M(row,:));
            end
        end

        r = r + 1;
        row = row + 1;
    end
end


function key = edge_key(i,j)
    a = min(i,j);
    b = max(i,j);
    key = sprintf('%d_%d', a, b);
end
