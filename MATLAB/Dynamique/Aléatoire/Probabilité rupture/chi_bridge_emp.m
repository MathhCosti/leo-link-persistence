%% facteur_liens_deleteres_depuis_analysis_temp.m
% Exécute analysis_temp, puis calcule automatiquement les grandeurs
% nécessaires au facteur correctif de p_break.
%
% Le facteur principal affiché est :
%   F_liens = nombre moyen de liens par composante
%
% La proportion de liens délétères est estimée par la proportion de ponts :
%   p_link_A = nombre total de ponts / nombre total de liens
%
% Le nombre moyen de liens délétères par composante vaut alors :
%   F_deletere = F_liens * p_link_A
%              = nombre total de ponts / nombre total de composantes

clear; clc; close all;

%% 1. Exécution du code de simulation
% Placer ce fichier dans le même dossier que analysis_temp(9).m,
% ou modifier le chemin ci-dessous.
analysis_file = fullfile(fileparts(mfilename('fullpath')), '..', 'analysis_temp.m');

if ~isfile(analysis_file)
    error('Fichier introuvable : %s', analysis_file);
end

% analysis_temp commence par clear : il efface donc analysis_file,
% mais son chemin a déjà été transmis à run.
run(analysis_file);

%% 2. Vérification des sorties nécessaires
if ~exist('Adjacency', 'var') || ~iscell(Adjacency)
    error('analysis_temp doit produire la cellule Adjacency.');
end

if ~exist('N', 'var') || ~isscalar(N) || N < 1
    error('analysis_temp doit produire le nombre de satellites N.');
end

Nt_factor = numel(Adjacency);

%% 3. Parcours de tous les graphes temporels
all_component_sizes = [];

n_components_t = zeros(Nt_factor,1);
n_edges_t      = zeros(Nt_factor,1);
n_bridges_t    = zeros(Nt_factor,1);
p_link_t       = zeros(Nt_factor,1);

for k_factor = 1:Nt_factor

    A_factor = spones(Adjacency{k_factor});
    A_factor = A_factor | A_factor.';
    A_factor = A_factor - diag(diag(A_factor));

    G_factor = graph(A_factor);

    % Composantes et tailles
    labels_factor = conncomp(G_factor);
    sizes_factor = accumarray(labels_factor(:), 1);

    n_components_t(k_factor) = numel(sizes_factor);
    all_component_sizes = [all_component_sizes; sizes_factor]; %#ok<AGROW>

    % Nombre de liens
    n_edges_t(k_factor) = numedges(G_factor);

    % Les ponts sont précisément les liens dont la suppression
    % augmente le nombre de composantes connexes.
    % MATLAB ne fournit pas bridges(G) dans toutes les versions.
    % On utilise dfsearch pour récupérer les arêtes de type "back" puis
    % un algorithme DFS de Tarjan implémenté dans la fonction locale.
    n_bridges_t(k_factor) = count_bridges_tarjan(A_factor);

    % Probabilité géométrique empirique qu'une paire soit liée
    if N >= 2
        p_link_t(k_factor) = n_edges_t(k_factor) / nchoosek(N,2);
    end
end

%% 4. Moments de la taille d'une composante typique
E_NA  = mean(all_component_sizes);
E_NA2 = mean(all_component_sizes.^2);

%% 5. Grandeurs agrégées cohérentes au niveau composante
N_comp_total   = sum(n_components_t);
N_edges_total  = sum(n_edges_t);
N_bridges_total = sum(n_bridges_t);

% Nombre moyen RÉEL de liens par composante
N_liens_moyen_emp = N_edges_total / N_comp_total;

% Proportion de liens qui sont des ponts (liens délétères)
if N_edges_total > 0
    p_link_A = N_bridges_total / N_edges_total;
else
    p_link_A = 0;
end

% Nombre moyen de liens délétères par composante
N_liens_deleteres = N_bridges_total / N_comp_total;

% Probabilité moyenne qu'une paire quelconque soit liée
p_link = mean(p_link_t);

%% 6. Estimation par la formule proposée
% Attention : cette expression traite les paires d'une composante comme
% si elles avaient la probabilité globale p_link d'être connectées.
% Elle est donc affichée à titre de comparaison.
N_paires_moyen = (E_NA2 - E_NA) / 2;
N_liens_moyen_formule = p_link * N_paires_moyen;
N_deleteres_formule = N_liens_moyen_formule * p_link_A;

%% 7. Affichage
fprintf('\n');
fprintf('=================================================================\n');
fprintf(' FACTEUR p_break CALCULE DEPUIS analysis_temp\n');
fprintf('=================================================================\n');
fprintf('Nombre de graphes temporels                : %d\n', Nt_factor);
fprintf('Nombre total de composantes observées      : %d\n', N_comp_total);
fprintf('Nombre total de liens observés             : %d\n', N_edges_total);
fprintf('Nombre total de ponts observés             : %d\n', N_bridges_total);
fprintf('-----------------------------------------------------------------\n');
fprintf('E[N_A]                                     : %.6f\n', E_NA);
fprintf('E[N_A^2]                                   : %.6f\n', E_NA2);
fprintf('p_link moyen                               : %.6e\n', p_link);
fprintf('p_link_A = proportion de ponts             : %.6f\n', p_link_A);
fprintf('-----------------------------------------------------------------\n');
fprintf('FACTEUR RECOMMANDE = liens/composante      : %.6f\n', N_liens_moyen_emp);
fprintf('Liens délétères moyens/composante          : %.6f\n', N_liens_deleteres);
fprintf('-----------------------------------------------------------------\n');
fprintf('Facteur obtenu par la formule analytique   : %.6f\n', ...
        N_liens_moyen_formule);
fprintf('Liens délétères par la formule             : %.6f\n', ...
        N_deleteres_formule);
fprintf('=================================================================\n\n');

%% 8. Tracés temporels
if exist('time_values', 'var') && numel(time_values) == Nt_factor
    x_factor = time_values(:);
    x_label_factor = 'Temps (s)';
else
    x_factor = (1:Nt_factor).';
    x_label_factor = 'Indice temporel';
end

figure;
plot(x_factor, n_edges_t ./ max(n_components_t,1), 'LineWidth', 1.5);
grid on;
xlabel(x_label_factor);
ylabel('Liens moyens par composante');
title('Facteur correctif temporel de p_{break}');

yline(N_liens_moyen_emp, '--', ...
    sprintf('Moyenne = %.3f', N_liens_moyen_emp), ...
    'LabelHorizontalAlignment', 'left');

figure;
plot(x_factor, n_bridges_t ./ max(n_components_t,1), 'LineWidth', 1.5);
grid on;
xlabel(x_label_factor);
ylabel('Ponts moyens par composante');
title('Nombre moyen de liens délétères');

yline(N_liens_deleteres, '--', ...
    sprintf('Moyenne = %.3f', N_liens_deleteres), ...
    'LabelHorizontalAlignment', 'left');

%% 9. Sauvegarde des résultats
save('chi_bridge_emp_results.mat', ...
    'E_NA', 'E_NA2', 'p_link', 'p_link_A', ...
    'N_liens_moyen_emp', 'N_liens_deleteres', ...
    'N_liens_moyen_formule', 'N_deleteres_formule', ...
    'n_components_t', 'n_edges_t', 'n_bridges_t', 'p_link_t');

fprintf('Résultats sauvegardés dans facteur_liens_deleteres_resultats.mat\n');

%% =================================================================
% Fonction locale : comptage des ponts par l'algorithme de Tarjan
%% =================================================================
function nBridges = count_bridges_tarjan(A)

    n = size(A,1);
    visited = false(n,1);
    discovery = zeros(n,1);
    low = zeros(n,1);
    parent = zeros(n,1);
    timer = 0;
    nBridges = 0;

    for startNode = 1:n
        if ~visited(startNode)
            dfs_bridge(startNode);
        end
    end

    function dfs_bridge(u)
        visited(u) = true;
        timer = timer + 1;
        discovery(u) = timer;
        low(u) = timer;

        neighbors = find(A(u,:));

        for idx = 1:numel(neighbors)
            v = neighbors(idx);

            if ~visited(v)
                parent(v) = u;
                dfs_bridge(v);
                low(u) = min(low(u), low(v));

                if low(v) > discovery(u)
                    nBridges = nBridges + 1;
                end

            elseif v ~= parent(u)
                low(u) = min(low(u), discovery(v));
            end
        end
    end
end
