clear; clc; close all;

%% Paramètres
N = 200;              % nombre de satellites
Re = 6371;            % rayon Terre [km]
h = 550;              % altitude [km]
R = Re + h;           % rayon orbital [km]

d_max = 4000;                         % distance max fixée [km]
alpha_values = linspace(1, 40, 40)*pi/180;   % valeurs de alpha_max [rad]

n_iter = 5;         % nombre de réalisations Monte-Carlo
rng(1);

%% Stockage des résultats simulés Monte-Carlo
Betti0_all = zeros(n_iter, length(alpha_values));
Betti1_graph_all = zeros(n_iter, length(alpha_values));
Betti1_complex_all = zeros(n_iter, length(alpha_values));

%% Stockage des approximations théoriques
E_theory = zeros(size(alpha_values));
I_theory = zeros(size(alpha_values));
Beta0_theory_sparse = zeros(size(alpha_values));
Beta0_theory_isolated = zeros(size(alpha_values));
Beta1_graph_theory_sparse = zeros(size(alpha_values));
Beta1_graph_theory_isolated = zeros(size(alpha_values));
Beta1_complex_theory_ER = zeros(size(alpha_values));

%% Angle imposé par la contrainte de distance
% d = 2R sin(alpha/2), donc d <= d_max équivaut à :
alpha_dmax = 2*asin(min(d_max/(2*R),1));

%% Théorie : ne dépend pas de la réalisation aléatoire
for k = 1:length(alpha_values)

    alpha_max = alpha_values(k);

    % Le graphe dépend en fait de l'angle effectif :
    % alpha_eff = min(alpha_max, angle associé à d_max)
    alpha_eff = min(alpha_max, alpha_dmax);

    % Probabilité que deux satellites soient liés
    p_link = (1 - cos(alpha_eff))/2;

    % Nombre moyen d'arêtes
    E_th = nchoosek(N,2) * p_link;

    % Nombre moyen de satellites isolés
    I_th = N * (1 - p_link)^(N-1);

    % Approximation sparse : beta0 ≈ N - E
    beta0_sparse = N - E_th;
    beta0_sparse = max(beta0_sparse, 1);
    beta0_sparse = min(beta0_sparse, N);

    % Approximation par satellites isolés, utile près du seuil de connectivité
    beta0_isolated = 1 + I_th;
    beta0_isolated = max(beta0_isolated, 1);
    beta0_isolated = min(beta0_isolated, N);

    % Betti-1 du graphe seul : beta1 = E - N + beta0
    beta1_graph_sparse = E_th - N + beta0_sparse;
    beta1_graph_sparse = max(beta1_graph_sparse, 0);

    beta1_graph_isolated = E_th - N + beta0_isolated;
    beta1_graph_isolated = max(beta1_graph_isolated, 0);

    % Approximation Erdős-Rényi grossière pour le complexe de clique :
    % on retire le nombre moyen de triangles remplis
    T_theory_ER = nchoosek(N,3) * p_link^3;
    beta1_complex_ER = beta1_graph_isolated - T_theory_ER;
    beta1_complex_ER = max(beta1_complex_ER, 0);

    % Stockage
    E_theory(k) = E_th;
    I_theory(k) = I_th;
    Beta0_theory_sparse(k) = beta0_sparse;
    Beta0_theory_isolated(k) = beta0_isolated;
    Beta1_graph_theory_sparse(k) = beta1_graph_sparse;
    Beta1_graph_theory_isolated(k) = beta1_graph_isolated;
    Beta1_complex_theory_ER(k) = beta1_complex_ER;
end

%% Boucle Monte-Carlo
for it = 1:n_iter

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

    %% Boucle sur alpha_max
    for k = 1:length(alpha_values)

        alpha_max = alpha_values(k);

        % Graphe de liens
        A = (D <= d_max) & (Alpha <= alpha_max);
        A(1:N+1:end) = false;
        A = A | A.';

        % Betti simulés
        [b0, b1_graph, b1_complex] = compute_betti_0_1(A);

        Betti0_all(it,k) = b0;
        Betti1_graph_all(it,k) = b1_graph;
        Betti1_complex_all(it,k) = b1_complex;
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

%% Axe en degrés
alpha_deg = alpha_values*180/pi;
alpha_dmax_deg = alpha_dmax*180/pi;

%% Affichage beta0
figure;
errorbar(alpha_deg, Betti0, Betti0_std, 'LineWidth', 1.5); hold on;
plot(alpha_deg, Beta0_theory_sparse, '--', 'LineWidth', 2);
plot(alpha_deg, Beta0_theory_isolated, ':', 'LineWidth', 2);
grid on;
xlabel('\alpha_{max} [deg]');
ylabel('\beta_0');
legend('Simulation moyenne \pm écart-type', ...
       'Théorie isolés', ...
       'Théorie connectés', ...
       'Location', 'best');
title(sprintf('\\beta_0 moyen en fonction de \\alpha_{max} — %d itérations', n_iter));

%% Affichage beta1 graphe seul
figure;
errorbar(alpha_deg, Betti1_graph, Betti1_graph_std, 'LineWidth', 1.5); hold on;
plot(alpha_deg, Beta1_graph_theory_sparse, '--', 'LineWidth', 2);
plot(alpha_deg, Beta1_graph_theory_isolated, ':', 'LineWidth', 2);
grid on;
xlabel('\alpha_{max} [deg]');
ylabel('\beta_1^{graphe}');
legend('Simulation moyenne \pm écart-type', ...
       'Théorie isolés', ...
       'Théorie connectés', ...
       'Location', 'best');
title(sprintf('\\beta_1 du graphe moyen en fonction de \\alpha_{max} — %d itérations', n_iter));

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
