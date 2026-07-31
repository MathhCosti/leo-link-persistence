%% p_bridge_emp.m
% Calcule empiriquement :
%
%   p_bridge_emp = P(pont | rupture)
%
% c'est-à-dire la proportion des liens qui disparaissent entre t et t+dt
% qui étaient des ponts dans le graphe juste avant leur disparition.
%
% Le script utilise directement le fichier analysis_temp_results.mat.

clear; clc; close all;

%% 1. Chargement des résultats de analysis_temp
this_dir = fileparts(mfilename('fullpath'));

candidate_files = {
    fullfile(this_dir, '..', 'analysis_temp_results.mat')
    fullfile(this_dir, 'analysis_temp_results(2).mat')
};

results_file = '';
for i_file = 1:numel(candidate_files)
    if isfile(candidate_files{i_file})
        results_file = candidate_files{i_file};
        break;
    end
end

if isempty(results_file)
    error(['Impossible de trouver analysis_temp_results.mat ou ', ...
           'analysis_temp_results(2).mat dans le dossier du script.']);
end

S = load(results_file);

if ~isfield(S, 'Adjacency') || ~iscell(S.Adjacency)
    error('Le fichier doit contenir la cellule Adjacency.');
end

Adjacency = S.Adjacency;
Nt = numel(Adjacency);

if Nt < 2
    error('Adjacency doit contenir au moins deux instants temporels.');
end

%% 2. Initialisation
n_transitions = Nt - 1;

n_removed_edges_t   = zeros(n_transitions, 1);
n_removed_bridges_t = zeros(n_transitions, 1);
p_bridge_bord_t     = nan(n_transitions, 1);

%% 3. Parcours des transitions t -> t + dt
for k = 1:n_transitions

    % Graphe juste avant la transition
    A_t = normalize_adjacency(Adjacency{k});

    % Graphe juste après la transition
    A_next = normalize_adjacency(Adjacency{k+1});

    if ~isequal(size(A_t), size(A_next))
        error('Les matrices Adjacency{%d} et Adjacency{%d} ont des tailles différentes.', ...
              k, k+1);
    end

    % Liens présents à t mais absents à t+dt
    removed_mask = A_t & ~A_next;

    % On ne compte chaque arête non orientée qu'une fois
    [removed_i, removed_j] = find(triu(removed_mask, 1));
    n_removed_edges_t(k) = numel(removed_i);

    % Détermination des ponts du graphe à l'instant t
    bridge_mask_t = bridge_mask_tarjan(A_t);

    % Parmi les liens supprimés, combien étaient des ponts à t ?
    if n_removed_edges_t(k) > 0
        linear_removed = sub2ind(size(A_t), removed_i, removed_j);
        n_removed_bridges_t(k) = nnz(bridge_mask_t(linear_removed));

        p_bridge_bord_t(k) = ...
            n_removed_bridges_t(k) / n_removed_edges_t(k);
    end
end

%% 4. Agrégation
total_removed_edges   = sum(n_removed_edges_t);
total_removed_bridges = sum(n_removed_bridges_t);

if total_removed_edges > 0
    p_bridge_emp = total_removed_bridges / total_removed_edges;
else
    p_bridge_emp = NaN;
    warning('Aucun lien rompu n''a été observé.');
end

% Moyenne temporelle non pondérée, donnée seulement pour comparaison
valid_t = ~isnan(p_bridge_bord_t);
if any(valid_t)
    p_bridge_bord_mean_t = mean(p_bridge_bord_t(valid_t));
else
    p_bridge_bord_mean_t = NaN;
end

%% 5. Calcul cohérent de chi_bridge_emp = P(pont | lien)
%
% Les ponts et les arêtes doivent être comptés sur les mêmes graphes G_t.
% On utilise ici les graphes avant chaque transition, k = 1,...,Nt-1,
% exactement comme pour le calcul de p_bridge_emp.

n_edges_t   = zeros(n_transitions, 1);
n_bridges_t = zeros(n_transitions, 1);
chi_bridge_t = nan(n_transitions, 1);

for k = 1:n_transitions

    A_t = normalize_adjacency(Adjacency{k});

    n_edges_t(k) = nnz(triu(A_t, 1));

    bridge_mask_t = bridge_mask_tarjan(A_t);
    n_bridges_t(k) = nnz(triu(bridge_mask_t, 1));

    if n_edges_t(k) > 0
        chi_bridge_t(k) = ...
            n_bridges_t(k) / n_edges_t(k);
    end
end

total_edges_observed   = sum(n_edges_t);
total_bridges_observed = sum(n_bridges_t);

if total_edges_observed > 0
    chi_bridge_emp = ...
        total_bridges_observed / total_edges_observed;
else
    chi_bridge_emp = NaN;
    warning('Aucun lien n''a été observé pour calculer chi_bridge_emp.');
end

valid_chi_t = ~isnan(chi_bridge_t);
if any(valid_chi_t)
    chi_bridge_mean_t = mean(chi_bridge_t(valid_chi_t));
else
    chi_bridge_mean_t = NaN;
end

%% 6. Affichage
fprintf('\n');
fprintf('=================================================================\n');
fprintf(' CALCUL EMPIRIQUE DE p_bridge_bord = P(pont | rupture)\n');
fprintf('=================================================================\n');
fprintf('Fichier chargé                              : %s\n', results_file);
fprintf('Nombre de transitions                      : %d\n', n_transitions);
fprintf('Nombre total de liens rompus               : %d\n', total_removed_edges);
fprintf('Nombre de liens rompus qui étaient des ponts: %d\n', total_removed_bridges);
fprintf('-----------------------------------------------------------------\n');
fprintf('p_bridge_emp = P(pont | rupture)       : %.6f\n', ...
        p_bridge_emp);
fprintf('Moyenne temporelle des rapports             : %.6f\n', ...
        p_bridge_bord_mean_t);

if ~isnan(chi_bridge_emp)
    fprintf('chi_bridge_emp = P(pont | lien)             : %.6f\n', ...
            chi_bridge_emp);
    fprintf('Moyenne temporelle de chi_bridge(t)         : %.6f\n', ...
            chi_bridge_mean_t);
    fprintf('Nombre total de liens observés              : %d\n', ...
            total_edges_observed);
    fprintf('Nombre total de ponts observés              : %d\n', ...
            total_bridges_observed);
    fprintf('Rapport p_bridge_emp / chi_bridge_emp       : %.6f\n', ...
            p_bridge_emp / chi_bridge_emp);
end

fprintf('=================================================================\n\n');

%% 7. Axe temporel
if isfield(S, 'time_values') && numel(S.time_values) == Nt
    x = S.time_values(1:end-1);
    x_label = 'Temps avant transition (s)';
else
    x = (1:n_transitions).';
    x_label = 'Indice de transition';
end

%% 8. Tracés
figure;
plot(x, p_bridge_bord_t, 'LineWidth', 1.4);
grid on;
xlabel(x_label);
ylabel('P(pont | rupture)');
title('Estimation temporelle de p_{bridge,bord}^{emp}');

if ~isnan(p_bridge_emp)
    yline(p_bridge_emp, '--', ...
        sprintf('Agrégée = %.4f', p_bridge_emp), ...
        'LabelHorizontalAlignment', 'left');
end

figure;
plot(x, chi_bridge_t, 'LineWidth', 1.4);
grid on;
xlabel(x_label);
ylabel('P(pont | lien)');
title('Estimation temporelle de \chi_{bridge}^{emp}');

if ~isnan(chi_bridge_emp)
    yline(chi_bridge_emp, '--', ...
        sprintf('Agrégée = %.4f', chi_bridge_emp), ...
        'LabelHorizontalAlignment', 'left');
end

figure;
plot(x, n_removed_edges_t, 'LineWidth', 1.2);
hold on;
plot(x, n_removed_bridges_t, 'LineWidth', 1.2);
grid on;
xlabel(x_label);
ylabel('Nombre de liens');
title('Liens rompus et ponts rompus par transition');
legend('Liens rompus', 'Ponts rompus', 'Location', 'best');

%% 9. Sauvegarde
save(fullfile(this_dir, 'p_bridge_emp_results.mat'), ...
    'p_bridge_emp', 'p_bridge_bord_mean_t', ...
    'chi_bridge_emp', 'chi_bridge_mean_t', ...
    'total_removed_edges', 'total_removed_bridges', ...
    'total_edges_observed', 'total_bridges_observed', ...
    'n_removed_edges_t', 'n_removed_bridges_t', ...
    'n_edges_t', 'n_bridges_t', ...
    'p_bridge_bord_t', 'chi_bridge_t', 'x');

fprintf('Résultats sauvegardés dans p_bridge_emp_results.mat\n');

%% ========================================================================
% Fonctions locales
%% ========================================================================

function A = normalize_adjacency(A_in)
    % Convertit une matrice d'adjacence en matrice logique, symétrique,
    % sans boucle et creuse.

    A = spones(A_in);
    A = logical(A | A.');
    A(1:size(A,1)+1:end) = false;
    A = sparse(A);
end

function bridge_mask = bridge_mask_tarjan(A)
    % Renvoie une matrice logique symétrique dont les entrées vraies
    % correspondent aux ponts du graphe non orienté.

    n = size(A, 1);

    visited   = false(n, 1);
    discovery = zeros(n, 1);
    low       = zeros(n, 1);
    parent    = zeros(n, 1);
    timer     = 0;

    bridge_i = zeros(max(nnz(triu(A,1)),1), 1);
    bridge_j = zeros(max(nnz(triu(A,1)),1), 1);
    n_bridges = 0;

    for start_node = 1:n
        if ~visited(start_node)
            dfs_bridge(start_node);
        end
    end

    bridge_i = bridge_i(1:n_bridges);
    bridge_j = bridge_j(1:n_bridges);

    bridge_mask = sparse( ...
        [bridge_i; bridge_j], ...
        [bridge_j; bridge_i], ...
        true, n, n);

    function dfs_bridge(u)
        visited(u) = true;
        timer = timer + 1;
        discovery(u) = timer;
        low(u) = timer;

        neighbors = find(A(u, :));

        for idx = 1:numel(neighbors)
            v = neighbors(idx);

            if ~visited(v)
                parent(v) = u;
                dfs_bridge(v);

                low(u) = min(low(u), low(v));

                if low(v) > discovery(u)
                    n_bridges = n_bridges + 1;
                    bridge_i(n_bridges) = u;
                    bridge_j(n_bridges) = v;
                end

            elseif v ~= parent(u)
                low(u) = min(low(u), discovery(v));
            end
        end
    end
end
