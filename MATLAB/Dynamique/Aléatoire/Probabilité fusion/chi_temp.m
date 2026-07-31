clear; clc; close all;

%% ============================================================
%  ETUDE MONTE-CARLO POUR UNE VALEUR UNIQUE DE lambda
%
%  Le code :
%    1) fixe une seule densité lambda ;
%    2) génère plusieurs réalisations indépendantes du graphe ;
%    3) calcule les grandeurs empiriques à chaque réalisation ;
%    4) conserve leur moyenne et leur écart-type ;
%    5) compare ces moyennes aux approximations théoriques.
%% ============================================================

%% Paramètres géométriques
Re = 6371;            % rayon terrestre [km]
h = 550;              % altitude [km]
R = Re + h;           % rayon orbital [km]

d_max = 1500;                 % distance maximale de liaison [km]
alpha_max = 20*pi/180;        % angle maximal [rad]

%% Densité unique
% lambda_scaled est exprimé en satellites / 10^6 km^2.
% Exemple : lambda_scaled = 0.4 correspond à lambda = 4e-7 sat/km^2.
lambda_scaled = 0.4;
lambda = lambda_scaled / 1e6;     % satellites / km^2

%% Nombre de réalisations Monte-Carlo
n_iter = 100;
rng(1);

%% Nombre de satellites associé à la densité
N = round(lambda * 4*pi*R^2);
N = max(N, 1);

%% Paramètre effectif de connexion
alpha_from_dmax = 2*asin(min(d_max/(2*R), 1));
alpha_eff = min(alpha_max, alpha_from_dmax);
p_link = (1 - cos(alpha_eff))/2;

%% Constantes géométriques
c2_union = 1.4135;
c3_conn  = 1.827;
c3_union = 1.80;

%% ============================================================
%  1. Valeurs théoriques
%% ============================================================

% Espérance du nombre d'arêtes
E_theory = N*(N-1)/2 * p_link;

% Composantes isolées
C1_theory = N * (1 - p_link)^(N-1);

% Dimères
if N >= 2
    C2_theory = nchoosek(N,2) * p_link ...
        * max(1 - c2_union*p_link, 0)^(N-2);
else
    C2_theory = 0;
end

% Trimères
if N >= 3
    C3_theory = nchoosek(N,3) * c3_conn * p_link^2 ...
        * max(1 - c3_union*p_link, 0)^(N-3);
else
    C3_theory = 0;
end

% Approximations de beta0
Beta0_theory_sparse = max(N - E_theory, 1);
Beta0_theory_isolated = min(max(2 + C1_theory, 1), N);
Beta0_theory_dimers = min(max(2 + C1_theory + C2_theory, 1), N);
Beta0_theory_trimers = min(max(2 + C1_theory + C2_theory + C3_theory, 1), N);

% Approximations de beta1 du graphe
Beta1_graph_theory_sparse = ...
    max(E_theory - N + Beta0_theory_sparse, 0);
Beta1_graph_theory_isolated = ...
    max(E_theory - N + Beta0_theory_isolated, 0);
Beta1_graph_theory_dimers = ...
    max(E_theory - N + Beta0_theory_dimers, 0);
Beta1_graph_theory_trimers = ...
    max(E_theory - N + Beta0_theory_trimers, 0);

% Fractions correctives théoriques : (N-beta0)/E
if E_theory > 0
    Frac_bridge_theory_sparse = max(0, min(1, ...
        (N - Beta0_theory_sparse) / E_theory));
    Frac_bridge_theory_isolated = max(0, min(1, ...
        (N - Beta0_theory_isolated) / E_theory));
    Frac_bridge_theory_dimers = max(0, min(1, ...
        (N - Beta0_theory_dimers) / E_theory));
    Frac_bridge_theory_trimers = max(0, min(1, ...
        (N - Beta0_theory_trimers) / E_theory));
else
    Frac_bridge_theory_sparse = NaN;
    Frac_bridge_theory_isolated = NaN;
    Frac_bridge_theory_dimers = NaN;
    Frac_bridge_theory_trimers = NaN;
end

%% ============================================================
%  2. Stockage des réalisations empiriques
%% ============================================================

Betti0_all = zeros(n_iter,1);
Betti1_graph_all = zeros(n_iter,1);
Betti1_complex_all = zeros(n_iter,1);
E_all = zeros(n_iter,1);
Frac_bridge_all = nan(n_iter,1);

%% ============================================================
%  3. Boucle Monte-Carlo
%% ============================================================

for it = 1:n_iter

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

    %% Betti et nombre d'arêtes
    [b0, b1_graph, b1_complex, E] = compute_betti_0_1(A);

    Betti0_all(it) = b0;
    Betti1_graph_all(it) = b1_graph;
    Betti1_complex_all(it) = b1_complex;
    E_all(it) = E;

    if E > 0
        Frac_bridge_all(it) = (N - b0) / E;
    end

    if mod(it, max(1, floor(n_iter/10))) == 0 || it == n_iter
        fprintf('Itération %d / %d terminée\n', it, n_iter);
    end
end

%% ============================================================
%  4. Moyennes empiriques
%% ============================================================

Betti0_emp = mean(Betti0_all);
Betti0_emp_std = std(Betti0_all);

Betti1_graph_emp = mean(Betti1_graph_all);
Betti1_graph_emp_std = std(Betti1_graph_all);

Betti1_complex_emp = mean(Betti1_complex_all);
Betti1_complex_emp_std = std(Betti1_complex_all);

E_emp = mean(E_all);
E_emp_std = std(E_all);

% Moyenne du rapport calculé réalisation par réalisation
Frac_bridge_emp_mean_ratio = mean(Frac_bridge_all, 'omitnan');
Frac_bridge_emp_std = std(Frac_bridge_all, 0, 'omitnan');

% Rapport des valeurs moyennes :
% (N - E[beta0]) / E[E]
if E_emp > 0
    Frac_bridge_emp_ratio_means = ...
        (N - Betti0_emp) / E_emp;
    Frac_bridge_emp_ratio_means = ...
        max(0, min(1, Frac_bridge_emp_ratio_means));
else
    Frac_bridge_emp_ratio_means = NaN;
end

%% ============================================================
%  5. Affichage console
%% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' MONTE-CARLO POUR UNE VALEUR UNIQUE DE lambda\n');
fprintf('============================================================\n');
fprintf('lambda_scaled                         : %.6f sat/10^6 km^2\n', ...
    lambda_scaled);
fprintf('lambda                                : %.8e sat/km^2\n', lambda);
fprintf('N                                     : %d\n', N);
fprintf('Nombre d''itérations                   : %d\n', n_iter);
fprintf('p_link                                : %.8e\n', p_link);
fprintf('------------------------------------------------------------\n');
fprintf('beta0 empirique moyen                 : %.8f\n', Betti0_emp);
fprintf('écart-type beta0                      : %.8f\n', Betti0_emp_std);
fprintf('E empirique moyen                     : %.8f\n', E_emp);
fprintf('écart-type E                          : %.8f\n', E_emp_std);
fprintf('beta1 graphe empirique moyen          : %.8f\n', Betti1_graph_emp);
fprintf('beta1 complexe empirique moyen        : %.8f\n', Betti1_complex_emp);
fprintf('------------------------------------------------------------\n');
fprintf('Fraction emp. moyenne des rapports    : %.8f\n', ...
    Frac_bridge_emp_mean_ratio);
fprintf('Fraction emp. rapport des moyennes    : %.8f\n', ...
    Frac_bridge_emp_ratio_means);
fprintf('------------------------------------------------------------\n');
fprintf('beta0 théorie connectés               : %.8f\n', ...
    Beta0_theory_sparse);
fprintf('beta0 théorie isolés                  : %.8f\n', ...
    Beta0_theory_isolated);
fprintf('beta0 théorie dimères                 : %.8f\n', ...
    Beta0_theory_dimers);
fprintf('beta0 théorie trimères                : %.8f\n', ...
    Beta0_theory_trimers);
fprintf('E théorique                           : %.8f\n', E_theory);
fprintf('============================================================\n\n');

%% Tableau récapitulatif
Modele = ["Empirique"; "Connectés"; "Isolés"; "Dimères"; "Trimères"];

Beta0_resume = [
    Betti0_emp
    Beta0_theory_sparse
    Beta0_theory_isolated
    Beta0_theory_dimers
    Beta0_theory_trimers
];

Fraction_resume = [
    Frac_bridge_emp_ratio_means
    Frac_bridge_theory_sparse
    Frac_bridge_theory_isolated
    Frac_bridge_theory_dimers
    Frac_bridge_theory_trimers
];

T = table(Modele, Beta0_resume, Fraction_resume);
disp(T);

%% ============================================================
%  6. Figures
%% ============================================================

% Distribution empirique de beta0
figure;
histogram(Betti0_all);
grid on;
xlabel('\beta_0');
ylabel('Nombre de réalisations');
title(sprintf(['Distribution empirique de \\beta_0 ', ...
    '(\\lambda = %.3f, %d itérations)'], ...
    lambda_scaled, n_iter));

xline(Betti0_emp, '--', ...
    sprintf('Moyenne = %.3f', Betti0_emp), ...
    'LabelHorizontalAlignment', 'left');

% Comparaison de beta0
figure;
beta0_values = [
    Betti0_emp
    Beta0_theory_sparse
    Beta0_theory_isolated
    Beta0_theory_dimers
    Beta0_theory_trimers
];

bar(beta0_values);
hold on;
errorbar(1, Betti0_emp, Betti0_emp_std, ...
    'LineStyle', 'none', 'LineWidth', 1.5);
grid on;
xticks(1:5);
xticklabels({'Empirique','Connectés','Isolés','Dimères','Trimères'});
ylabel('\beta_0');
title(sprintf('\\beta_0 pour \\lambda = %.3f', lambda_scaled));

% Comparaison des fractions correctives
figure;
frac_values = [
    Frac_bridge_emp_ratio_means
    Frac_bridge_theory_sparse
    Frac_bridge_theory_isolated
    Frac_bridge_theory_dimers
    Frac_bridge_theory_trimers
];

bar(frac_values);
hold on;
errorbar(1, Frac_bridge_emp_ratio_means, Frac_bridge_emp_std, ...
    'LineStyle', 'none', 'LineWidth', 1.5);
grid on;
xticks(1:5);
xticklabels({'Empirique','Connectés','Isolés','Dimères','Trimères'});
ylabel('(N-\beta_0)/|E|');
ylim([0 1.05]);
title(sprintf('Fraction corrective pour \\lambda = %.3f', lambda_scaled));

%% ============================================================
%  7. Sauvegarde
%% ============================================================

save('chi_temp_results.mat', ...
    'lambda_scaled', 'lambda', 'N', 'n_iter', ...
    'R', 'd_max', 'alpha_max', 'alpha_eff', 'p_link', ...
    'Betti0_all', 'Betti1_graph_all', ...
    'Betti1_complex_all', 'E_all', 'Frac_bridge_all', ...
    'Betti0_emp', 'Betti0_emp_std', ...
    'Betti1_graph_emp', 'Betti1_graph_emp_std', ...
    'Betti1_complex_emp', 'Betti1_complex_emp_std', ...
    'E_emp', 'E_emp_std', ...
    'Frac_bridge_emp_mean_ratio', ...
    'Frac_bridge_emp_ratio_means', ...
    'Frac_bridge_emp_std', ...
    'E_theory', ...
    'C1_theory', 'C2_theory', 'C3_theory', ...
    'Beta0_theory_sparse', ...
    'Beta0_theory_isolated', ...
    'Beta0_theory_dimers', ...
    'Beta0_theory_trimers', ...
    'Beta1_graph_theory_sparse', ...
    'Beta1_graph_theory_isolated', ...
    'Beta1_graph_theory_dimers', ...
    'Beta1_graph_theory_trimers', ...
    'Frac_bridge_theory_sparse', ...
    'Frac_bridge_theory_isolated', ...
    'Frac_bridge_theory_dimers', ...
    'Frac_bridge_theory_trimers', ...
    'c2_union', 'c3_conn', 'c3_union');

fprintf('Résultats sauvegardés dans chi_temp_lambda_single_results.mat\n');

%% ============================================================
%  Fonctions
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
