clear; clc; close all;

%% Paramètres
N = 200;              % nombre de satellites
Re = 6371;            % rayon Terre [km]
h = 550;              % altitude [km]
R = Re + h;           % rayon orbital [km]

alpha_max = 20*pi/180;                    % angle max fixé [rad]
dmax_values = linspace(300, 2500, 40);    % valeurs de d_max [km]

n_iter = 100;                             % nombre d'itérations Monte-Carlo
rng(1);                                   % reproductibilité

%% Stockage des résultats simulés pour toutes les itérations
Betti0_all = zeros(n_iter, length(dmax_values));
Betti1_graph_all = zeros(n_iter, length(dmax_values));
Betti1_complex_all = zeros(n_iter, length(dmax_values));

%% Seuil imposé par alpha_max exprimé en distance corde
% Au-delà de cette distance, c'est alpha_max qui limite les liens.
d_alpha_max = 2*R*sin(alpha_max/2); %#ok<NASGU>

%% Boucle Monte-Carlo
for it = 1:n_iter

    fprintf('Itération Monte-Carlo %d / %d\n', it, n_iter);

    %% Génération uniforme des satellites sur une sphère
    u = rand(N,1);
    v = rand(N,1);

    theta = 2*pi*u;
    phi = acos(2*v - 1);

    x = R * sin(phi).*cos(theta);
    y = R * sin(phi).*sin(theta);
    z = R * cos(phi);

    P = [x y z];

    %% Matrices de distance et d'angle
    D = squareform(pdist(P));

    U = P ./ vecnorm(P,2,2);
    CosAlpha = U * U.';
    CosAlpha = max(min(CosAlpha,1),-1);
    Alpha = acos(CosAlpha);

    %% Boucle sur d_max
    for k = 1:length(dmax_values)

        d_max = dmax_values(k);

        %% Graphe de liens simulé
        A = (D <= d_max) & (Alpha <= alpha_max);
        A(1:N+1:end) = false;
        A = A | A.';

        [b0, b1_graph, b1_complex] = compute_betti_0_1(A);

        Betti0_all(it,k) = b0;
        Betti1_graph_all(it,k) = b1_graph;
        Betti1_complex_all(it,k) = b1_complex;
    end
end

%% Moyenne Monte-Carlo
Betti0 = mean(Betti0_all, 1);
Betti1_graph = mean(Betti1_graph_all, 1);
Betti1_complex = mean(Betti1_complex_all, 1);

%% Écart-type Monte-Carlo, utile pour visualiser la dispersion
Betti0_std = std(Betti0_all, 0, 1);
Betti1_graph_std = std(Betti1_graph_all, 0, 1);
Betti1_complex_std = std(Betti1_complex_all, 0, 1);

%% Stockage des approximations théoriques
Beta0_theory_sparse = zeros(size(dmax_values));
Beta0_theory_isolated = zeros(size(dmax_values));
Beta1_graph_theory_sparse = zeros(size(dmax_values));
Beta1_graph_theory_isolated = zeros(size(dmax_values));

%% Calcul des approximations théoriques probabilistes
for k = 1:length(dmax_values)

    d_max = dmax_values(k);

    % La contrainte effective est l'angle le plus restrictif entre :
    %   - alpha_max fixé ;
    %   - l'angle équivalent à d_max.
    alpha_from_dmax = 2*asin(min(d_max/(2*R), 1));
    alpha_eff = min(alpha_max, alpha_from_dmax);

    % Probabilité qu'une paire de satellites soit reliée.
    p_link = (1 - cos(alpha_eff))/2;

    % Nombre moyen d'arêtes.
    E_theory = N*(N-1)/2 * p_link;

    % Nombre moyen de satellites isolés.
    I_theory = N * (1 - p_link)^(N-1);

    % Approximation sparse : chaque arête fusionne deux composantes.
    beta0_sparse = N - E_theory;
    beta0_sparse = max(beta0_sparse, 1);

    % Approximation proche de la connectivité : composante géante + isolés.
    beta0_isolated = 1 + I_theory;

    Beta0_theory_sparse(k) = beta0_sparse;
    Beta0_theory_isolated(k) = beta0_isolated;

    % beta1 graphe = E - N + beta0.
    Beta1_graph_theory_sparse(k) = max(E_theory - N + beta0_sparse, 0);
    Beta1_graph_theory_isolated(k) = max(E_theory - N + beta0_isolated, 0);
end

%% Figure 1 : beta0 simulation moyenne vs théorie
figure;
errorbar(dmax_values, Betti0, Betti0_std, 'LineWidth', 1.5); hold on;
plot(dmax_values, Beta0_theory_sparse, '--', 'LineWidth', 2);
plot(dmax_values, Beta0_theory_isolated, ':', 'LineWidth', 2);
grid on;
xlabel('d_{max} [km]');
ylabel('beta 0');
legend('Simulation moyenne \pm écart-type', 'Théorie isolé', ...
    'Théorie connecté', 'Location', 'best');
title(sprintf('Moyenne Monte-Carlo de beta 0 en fonction de d_{max} (%d itérations)', n_iter));

%% Figure 2 : beta1 graphe seul simulation moyenne vs théorie
figure;
errorbar(dmax_values, Betti1_graph, Betti1_graph_std, 'LineWidth', 1.5); hold on;
plot(dmax_values, Beta1_graph_theory_sparse, '--', 'LineWidth', 2);
plot(dmax_values, Beta1_graph_theory_isolated, ':', 'LineWidth', 2);
grid on;
xlabel('d_{max} [km]');
ylabel('beta 1^{graphe}');
legend('Simulation moyenne \pm écart-type', 'Théorie isolés', 'Théorie connectés', ...
    'Location', 'best');
title(sprintf('Moyenne Monte-Carlo de beta 1 du graphe en fonction de d_{max} (%d itérations)', n_iter));

%% ============================================================
%% Fonction Betti
%% ============================================================

function [beta0, beta1_graph, beta1_complex] = compute_betti_0_1(A)

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
