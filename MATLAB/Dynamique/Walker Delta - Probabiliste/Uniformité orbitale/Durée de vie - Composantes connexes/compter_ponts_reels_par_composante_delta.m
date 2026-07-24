%% compter_ponts_reels_par_composante_delta.m
% Compte les vrais ponts de chaque composante connexe à chaque instant
% dans le fichier leo_zigzag_analysis_results_delta.mat.
%
% Un pont est une arête dont la suppression augmente le nombre de
% composantes connexes.
%
% Le script calcule :
% - le nombre de ponts de chaque composante ;
% - le nombre moyen de ponts par composante ;
% - la fraction réelle de liens qui sont des ponts ;
% - les moyennes conditionnelles aux composantes non triviales ;
% - l'évolution temporelle de ces grandeurs.

clearvars;
clc;
close all;

%% ============================ FICHIER ================================

script_dir = fileparts(mfilename('fullpath'));
data_file = fullfile(script_dir, 'leo_zigzag_analysis_results_delta.mat');

if ~isfile(data_file)
    error('Fichier introuvable : %s', data_file);
end

S = load(data_file, 'Adjacency', 'time_values', 'N', 'dmax', 'dt');

required_fields = {'Adjacency','time_values','N'};
for k = 1:numel(required_fields)
    if ~isfield(S, required_fields{k})
        error('Variable manquante dans le .mat : %s', required_fields{k});
    end
end

Adjacency = S.Adjacency;
time_values = S.time_values(:);
N = double(S.N);

Nt = numel(time_values);

if numel(Adjacency) ~= Nt
    error('Adjacency et time_values ont des tailles incompatibles.');
end

%% ===================== VARIABLES DE STOCKAGE ========================

n_components_t = zeros(Nt,1);
n_nontrivial_components_t = zeros(Nt,1);
n_components_with_bridge_t = zeros(Nt,1);

n_edges_t = zeros(Nt,1);
n_bridges_t = zeros(Nt,1);

mean_bridges_per_component_t = NaN(Nt,1);
mean_bridges_per_nontrivial_component_t = NaN(Nt,1);
mean_bridges_given_positive_t = NaN(Nt,1);

bridge_fraction_t = NaN(Nt,1);
mean_edges_per_component_t = NaN(Nt,1);

% Une ligne par composante observée :
% [temps_index, composante_locale, taille, nb_liens, nb_ponts]
component_table_raw = [];

%% ========================= BOUCLE TEMPORELLE ========================

for t = 1:Nt

    A = logical(Adjacency{t});

    if ~isequal(size(A), [N N])
        error('Adjacency{%d} doit être une matrice N x N.', t);
    end

    % Sécurité : matrice simple, symétrique, sans boucle.
    A = A | A.';
    A(1:N+1:end) = false;

    G = graph(A);
    labels = conncomp(G);
    component_ids = unique(labels);

    n_components_t(t) = numel(component_ids);
    n_edges_t(t) = numedges(G);

    bridges_per_component = zeros(numel(component_ids),1);
    edges_per_component = zeros(numel(component_ids),1);
    sizes_per_component = zeros(numel(component_ids),1);

    for c = 1:numel(component_ids)

        vertices = find(labels == component_ids(c));
        Gc = subgraph(G, vertices);

        n_vertices_c = numnodes(Gc);
        n_edges_c = numedges(Gc);

        % MATLAB récent : findedge après biconncomp(...,'OutputForm','cell')
        % n'est pas nécessaire. On utilise une fonction Tarjan locale,
        % compatible avec les versions ne disposant pas d'une fonction
        % bridges dédiée.
        bridge_pairs_local = find_bridges_tarjan(adjacency(Gc));

        n_bridges_c = size(bridge_pairs_local,1);

        sizes_per_component(c) = n_vertices_c;
        edges_per_component(c) = n_edges_c;
        bridges_per_component(c) = n_bridges_c;

        component_table_raw = [component_table_raw; ...
            t, c, n_vertices_c, n_edges_c, n_bridges_c]; %#ok<AGROW>
    end

    n_bridges_t(t) = sum(bridges_per_component);

    nontrivial_mask = sizes_per_component >= 2;
    positive_bridge_mask = bridges_per_component > 0;

    n_nontrivial_components_t(t) = sum(nontrivial_mask);
    n_components_with_bridge_t(t) = sum(positive_bridge_mask);

    mean_bridges_per_component_t(t) = ...
        mean(bridges_per_component, 'omitnan');

    if any(nontrivial_mask)
        mean_bridges_per_nontrivial_component_t(t) = ...
            mean(bridges_per_component(nontrivial_mask), 'omitnan');
    end

    if any(positive_bridge_mask)
        mean_bridges_given_positive_t(t) = ...
            mean(bridges_per_component(positive_bridge_mask), 'omitnan');
    end

    if n_edges_t(t) > 0
        bridge_fraction_t(t) = n_bridges_t(t) / n_edges_t(t);
    end

    if n_components_t(t) > 0
        mean_edges_per_component_t(t) = ...
            n_edges_t(t) / n_components_t(t);
    end
end

%% =========================== TABLE ================================

component_table = array2table(component_table_raw, ...
    'VariableNames', { ...
    'time_index', ...
    'component_index', ...
    'component_size', ...
    'n_edges', ...
    'n_bridges'});

component_table.time = time_values(component_table.time_index);

%% ====================== STATISTIQUES GLOBALES =======================

total_components = sum(n_components_t);
total_nontrivial_components = sum(n_nontrivial_components_t);
total_components_with_bridge = sum(n_components_with_bridge_t);

total_edges = sum(n_edges_t);
total_bridges = sum(n_bridges_t);

% Moyenne globale pondérée par le nombre de composantes observées.
mean_bridges_per_component = ...
    total_bridges / max(total_components,1);

nontrivial_rows = component_table.component_size >= 2;
positive_rows = component_table.n_bridges > 0;

mean_bridges_per_nontrivial_component = ...
    mean(component_table.n_bridges(nontrivial_rows), 'omitnan');

mean_bridges_given_positive = ...
    mean(component_table.n_bridges(positive_rows), 'omitnan');

bridge_fraction = total_bridges / max(total_edges,1);

mean_edges_per_component = total_edges / max(total_components,1);

fraction_components_with_bridge = ...
    total_components_with_bridge / max(total_components,1);

fraction_nontrivial_components_with_bridge = ...
    sum(nontrivial_rows & positive_rows) / max(sum(nontrivial_rows),1);

%% ============================ AFFICHAGE =============================

fprintf('\n');
fprintf('====================================================================\n');
fprintf(' Ponts réels par composante - Walker Delta temporel\n');
fprintf('====================================================================\n');
fprintf('N                                      : %d\n', N);
fprintf('Nombre d''instants                     : %d\n', Nt);
fprintf('Nombre total de composantes observées  : %d\n', total_components);
fprintf('Composantes non triviales observées    : %d\n', ...
    total_nontrivial_components);
fprintf('Nombre total de liens observés         : %d\n', total_edges);
fprintf('Nombre total de ponts observés         : %d\n', total_bridges);
fprintf('--------------------------------------------------------------------\n');
fprintf('Ponts moyens / composante              : %.6f\n', ...
    mean_bridges_per_component);
fprintf('Ponts moyens / composante non triviale : %.6f\n', ...
    mean_bridges_per_nontrivial_component);
fprintf('Ponts moyens sachant B>0               : %.6f\n', ...
    mean_bridges_given_positive);
fprintf('Liens moyens / composante              : %.6f\n', ...
    mean_edges_per_component);
fprintf('Fraction réelle de liens ponts         : %.6f\n', ...
    bridge_fraction);
fprintf('Fraction des composantes avec un pont  : %.6f\n', ...
    fraction_components_with_bridge);
fprintf('Fraction non triviales avec un pont    : %.6f\n', ...
    fraction_nontrivial_components_with_bridge);
fprintf('====================================================================\n\n');

%% ============================= GRAPHES ==============================

figure;
plot(time_values, mean_bridges_per_component_t, 'LineWidth', 1.5);
grid on;
xlabel('Temps (s)');
ylabel('Nombre moyen de ponts');
title('Nombre moyen de ponts par composante');

figure;
plot(time_values, bridge_fraction_t, 'LineWidth', 1.5);
grid on;
xlabel('Temps (s)');
ylabel('Fraction de liens qui sont des ponts');
title('Fraction réelle de liens critiques');

figure;
histogram(component_table.n_bridges, ...
    'BinMethod', 'integers', ...
    'Normalization', 'probability');
grid on;
xlabel('Nombre de ponts dans la composante');
ylabel('Probabilité empirique');
title('Distribution du nombre de ponts par composante');

figure;
scatter(component_table.component_size, ...
    component_table.n_bridges, 12, 'filled');
grid on;
xlabel('Taille de la composante');
ylabel('Nombre de ponts');
title('Ponts en fonction de la taille des composantes');

%% ============================= SAUVEGARDE ===========================

results = struct();

results.N = N;
results.time_values = time_values;

if isfield(S,'dmax')
    results.dmax = S.dmax;
end

if isfield(S,'dt')
    results.dt = S.dt;
end

results.component_table = component_table;

results.n_components_t = n_components_t;
results.n_nontrivial_components_t = n_nontrivial_components_t;
results.n_components_with_bridge_t = n_components_with_bridge_t;
results.n_edges_t = n_edges_t;
results.n_bridges_t = n_bridges_t;

results.mean_bridges_per_component_t = ...
    mean_bridges_per_component_t;
results.mean_bridges_per_nontrivial_component_t = ...
    mean_bridges_per_nontrivial_component_t;
results.mean_bridges_given_positive_t = ...
    mean_bridges_given_positive_t;
results.bridge_fraction_t = bridge_fraction_t;
results.mean_edges_per_component_t = mean_edges_per_component_t;

results.total_components = total_components;
results.total_nontrivial_components = total_nontrivial_components;
results.total_edges = total_edges;
results.total_bridges = total_bridges;

results.mean_bridges_per_component = ...
    mean_bridges_per_component;
results.mean_bridges_per_nontrivial_component = ...
    mean_bridges_per_nontrivial_component;
results.mean_bridges_given_positive = ...
    mean_bridges_given_positive;
results.mean_edges_per_component = ...
    mean_edges_per_component;
results.bridge_fraction = bridge_fraction;
results.fraction_components_with_bridge = ...
    fraction_components_with_bridge;
results.fraction_nontrivial_components_with_bridge = ...
    fraction_nontrivial_components_with_bridge;

save('ponts_reels_par_composante_delta_results.mat', ...
    'results', '-v7.3');

writetable(component_table, ...
    'ponts_reels_par_composante_delta.csv');

fprintf('Résultats sauvegardés dans :\n');
fprintf('  ponts_reels_par_composante_delta_results.mat\n');
fprintf('  ponts_reels_par_composante_delta.csv\n');

%% ========================== FONCTION LOCALE =========================

function bridge_pairs = find_bridges_tarjan(A)
% Trouve tous les ponts d'un graphe simple non orienté par l'algorithme
% de Tarjan.
%
% Entrée :
%   A : matrice d'adjacence logique ou numérique.
%
% Sortie :
%   bridge_pairs : matrice B x 2 contenant les indices locaux des ponts.

    A = logical(A);
    n = size(A,1);

    visited = false(n,1);
    discovery = zeros(n,1);
    low = zeros(n,1);
    parent = zeros(n,1);

    timer = 0;
    bridge_pairs = zeros(0,2);

    for root = 1:n
        if ~visited(root)
            dfs(root);
        end
    end

    function dfs(u)
        visited(u) = true;
        timer = timer + 1;
        discovery(u) = timer;
        low(u) = timer;

        neighbors = find(A(u,:));

        for idx = 1:numel(neighbors)
            v = neighbors(idx);

            if ~visited(v)
                parent(v) = u;
                dfs(v);

                low(u) = min(low(u), low(v));

                if low(v) > discovery(u)
                    bridge_pairs(end+1,:) = sort([u v]); %#ok<AGROW>
                end

            elseif v ~= parent(u)
                low(u) = min(low(u), discovery(v));
            end
        end
    end
end
