clear; clc; close all;

%% ============================================================
%  ÉTUDE TOPOLOGIQUE TEMPORELLE D'UN RÉSEAU LEO
%  Mouvement avec vecteurs aléatoires tangentiels
%  Sans animation
%
%  Sorties :
%  - beta0(t) : nombre de composantes connexes
%  - beta1_graph(t) : nombre de cycles du graphe non rempli
%  - beta0 et beta1 sur la suite zigzag par unions
%
%  Zigzag construit :
%  G1 -> G1 union G2 <- G2 -> G2 union G3 <- G3 ...
%% ============================================================

%% Paramètres physiques
R_earth = 6371;      % km
h = 550;             % km
R = R_earth + h;     % rayon orbital

mu = 398600;              % km^3/s^2
omega = sqrt(mu / R^3);   % vitesse angulaire rad/s

%% Paramètres du processus de Poisson
lambda = 4e-7;       % satellites / km^2
surface_sphere = 4*pi*R^2;

N = poissrnd(lambda * surface_sphere);

fprintf('Nombre de satellites générés : N = %d\n', N);

%% Génération uniforme des positions initiales sur la sphère
u = rand(N,1);
phi = 2*pi*rand(N,1);
theta = acos(1 - 2*u);

x = R * sin(theta) .* cos(phi);
y = R * sin(theta) .* sin(phi);
z = R * cos(theta);

positions0 = [x y z];

%% ============================================================
%  Directions aléatoires tangentes à la sphère
%% ============================================================

% Vecteurs radiaux unitaires
r0 = positions0 / R;              % N x 3

% Vecteurs aléatoires 3D quelconques
a = randn(N,3);

% Projection dans le plan tangent à la sphère :
% v = a - (a.r0) r0
v = a - sum(a .* r0, 2) .* r0;

% Normalisation : chaque v est unitaire et tangent à la sphère
v = v ./ vecnorm(v, 2, 2);

% Optionnel : sens aléatoire +1 / -1
sens = sign(rand(N,1) - 0.5);
sens(sens == 0) = 1;
v = sens .* v;

%% Paramètres des liens et du temps
dmax = 1500;     % km
dt = 3;         % pas temporel en secondes
Tmax = 1500;    % durée totale de simulation

time_values = 0:dt:Tmax;
Nt = length(time_values);

%% Stockage
Positions = cell(Nt,1);
Adjacency = cell(Nt,1);

num_edges = zeros(Nt,1);
beta0 = zeros(Nt,1);
beta1_graph = zeros(Nt,1);
largest_component = zeros(Nt,1);

% Distribution empirique des tailles de composantes :
% component_size_counts(k,s) = nombre de composantes de taille s au temps k
component_size_counts = zeros(Nt,N);

%% ============================================================
%  1. CONSTRUCTION DES GRAPHES TEMPORELS G(t)
%% ============================================================

for k = 1:Nt

    t = time_values(k);

    %% Mouvement sur grand cercle avec direction tangentielle aléatoire
    % Formule géodésique sur la sphère :
    % r(t) = R [ r0 cos(omega t) + v sin(omega t) ]
    positions_t = R * (r0 * cos(omega*t) + v * sin(omega*t));

    %% Graphe de lien
    D = squareform(pdist(positions_t));
    A = (D <= dmax) & (D > 0);
    A = sparse(A);

    %% Stockage
    Positions{k} = positions_t;
    Adjacency{k} = A;

    %% Mesures topologiques sur le graphe
    G = graph(A);

    comp = conncomp(G);
    beta0(k) = max(comp);

    comp_sizes = accumarray(comp', 1);
    largest_component(k) = max(comp_sizes);

    % Histogramme exact des tailles de composantes a cet instant
    % Ex. si les tailles sont [1 1 2 5], on stocke
    % C_1=2, C_2=1, C_5=1.
    counts_by_size = accumarray(comp_sizes,1,[N 1]);
    component_size_counts(k,:) = counts_by_size.';

    E = nnz(triu(A,1));
    num_edges(k) = E;

    % Nombre cyclomatique du graphe :
    % beta1 = E - V + C
    beta1_graph(k) = E - N + beta0(k);
end

%% ============================================================
%  1.b DISTRIBUTION MOYENNE DES TAILLES DE COMPOSANTES
%% ============================================================

% Nombre moyen de composantes de taille s :
mean_component_count_by_size = mean(component_size_counts,1);

% Nombre moyen de sommets appartenant a des composantes de taille s :
s_values_component = 1:N;
mean_nodes_by_component_size = ...
    s_values_component .* mean_component_count_by_size;

% Verification exacte :
% pour chaque temps, somme_s C_s(t) = beta0(t)
beta0_from_component_sizes = sum(component_size_counts,2);
max_component_count_error = ...
    max(abs(beta0_from_component_sizes-beta0));

% Et somme_s s*C_s(t) = N
N_from_component_sizes = ...
    component_size_counts * (1:N).';
max_component_mass_error = ...
    max(abs(N_from_component_sizes-N));

fprintf('\n');
fprintf('====================================================================\n');
fprintf(' Distribution des tailles de composantes\n');
fprintf('====================================================================\n');
fprintf('beta0 moyen                         : %.8f\n', mean(beta0));
fprintf('Somme_s E[C_s]                     : %.8f\n', ...
    sum(mean_component_count_by_size));
fprintf('Erreur max sum_s C_s(t) - beta0(t) : %.3e\n', ...
    max_component_count_error);
fprintf('Erreur max sum_s s C_s(t) - N      : %.3e\n', ...
    max_component_mass_error);
fprintf('====================================================================\n');

% On ne trace que jusqu'a la derniere taille effectivement observee.
last_nonzero_size = find(mean_component_count_by_size>0,1,'last');

if ~isempty(last_nonzero_size)
    figure;
    stem(1:last_nonzero_size, ...
        mean_component_count_by_size(1:last_nonzero_size), ...
        'filled');
    grid on;
    xlabel('Taille s de la composante');
    ylabel('Nombre moyen de composantes E[C_s]');
    title('Distribution empirique des tailles de composantes');

    figure;
    stem(1:last_nonzero_size, ...
        mean_nodes_by_component_size(1:last_nonzero_size), ...
        'filled');
    grid on;
    xlabel('Taille s de la composante');
    ylabel('Nombre moyen de satellites s E[C_s]');
    title('Masse de sommets par taille de composante');
end

%% ============================================================
%  2. GRAPHES TEMPORELS CLASSIQUES
%% ============================================================

figure;
plot(time_values, beta0, 'LineWidth', 2);
grid on;
xlabel('Temps (s)');
ylabel('\beta_0');
title('\beta_0(t) : nombre de composantes connexes');

figure;
plot(time_values, beta1_graph, 'LineWidth', 2);
grid on;
xlabel('Temps (s)');
ylabel('\beta_1 graphe');
title('\beta_1(t) du graphe non rempli');

figure;
plot(time_values, largest_component / N, 'LineWidth', 2);
grid on;
xlabel('Temps (s)');
ylabel('|C_{max}| / N');
title('Fraction de satellites dans la plus grande composante');

figure;
plot(time_values, num_edges, 'LineWidth', 2);
grid on;
xlabel('Temps (s)');
ylabel('Nombre de liens');
title('Nombre de liens inter-satellites');

%% ============================================================
%  3. CONSTRUCTION DU ZIGZAG PAR UNIONS
%
%  G1 -> G1 U G2 <- G2 -> G2 U G3 <- G3 ...
%% ============================================================

Nz = 2*Nt - 1;

ZigzagAdjacency = cell(Nz,1);
ZigzagLabels = zeros(Nz,1);

idx = 1;

for k = 1:Nt

    % Graphe réel G_k
    ZigzagAdjacency{idx} = Adjacency{k};
    ZigzagLabels(idx) = k;
    idx = idx + 1;

    % Graphe union G_k U G_{k+1}
    if k < Nt
        ZigzagAdjacency{idx} = Adjacency{k} | Adjacency{k+1};
        ZigzagLabels(idx) = k + 0.5;
        idx = idx + 1;
    end
end

%% ============================================================
%  4. BETTI SUR LA SUITE ZIGZAG
%% ============================================================

beta0_zigzag = zeros(Nz,1);
beta1_zigzag_graph = zeros(Nz,1);
num_edges_zigzag = zeros(Nz,1);
largest_component_zigzag = zeros(Nz,1);

for k = 1:Nz

    A = ZigzagAdjacency{k};
    G = graph(A);

    comp = conncomp(G);
    beta0_zigzag(k) = max(comp);

    comp_sizes = accumarray(comp', 1);
    largest_component_zigzag(k) = max(comp_sizes);

    E = nnz(triu(A,1));
    num_edges_zigzag(k) = E;

    beta1_zigzag_graph(k) = E - N + beta0_zigzag(k);
end

%% ============================================================
%  4.b PONTS EXPOSES DANS LES GRAPHES UNION
%
%  Pour chaque transition
%
%      G_k -> U_k <- G_{k+1},
%      U_k = G_k union G_{k+1},
%
%  on calcule :
%
%      B_k       = nombre de vrais ponts de U_k,
%      beta0(U_k)= nombre de composantes exposees,
%
%  puis
%
%      Bbar_expose
%        = sum_k B_k / sum_k beta0(U_k).
%
%  On mesure aussi
%
%      q_break|pont
%        = nombre de ponts de U_k absents de G_{k+1}
%          / nombre total de ponts de U_k.
%
%  Enfin :
%
%      p_break_ponts
%        = 1-(1-q_break|pont)^Bbar_expose.
%% ============================================================

n_bridges_union = zeros(Nt-1,1);
n_removed_bridges_union = zeros(Nt-1,1);
beta0_union = zeros(Nt-1,1);

mean_bridges_per_union_component_t = NaN(Nt-1,1);
q_break_given_bridge_t = NaN(Nt-1,1);

for k = 1:Nt-1

    A_union = logical(Adjacency{k} | Adjacency{k+1});
    A_next  = logical(Adjacency{k+1});

    G_union = graph(A_union);
    comp_union = conncomp(G_union);
    beta0_union(k) = max(comp_union);

    bridge_pairs_union = find_bridges_tarjan(A_union);
    n_bridges_union(k) = size(bridge_pairs_union,1);

    if beta0_union(k) > 0
        mean_bridges_per_union_component_t(k) = ...
            n_bridges_union(k) / beta0_union(k);
    end

    if n_bridges_union(k) > 0
        bridge_idx = sub2ind([N N], ...
            bridge_pairs_union(:,1), ...
            bridge_pairs_union(:,2));

        removed_mask = ~A_next(bridge_idx);
        n_removed_bridges_union(k) = sum(removed_mask);

        q_break_given_bridge_t(k) = ...
            n_removed_bridges_union(k) / n_bridges_union(k);
    end
end

total_union_components = sum(beta0_union);
total_union_bridges = sum(n_bridges_union);
total_removed_union_bridges = sum(n_removed_bridges_union);

mean_bridges_per_exposed_component = ...
    total_union_bridges / max(total_union_components,1);

q_break_given_bridge_global = ...
    total_removed_union_bridges / max(total_union_bridges,1);

p_break_from_bridges = ...
    1 - (1-q_break_given_bridge_global)^ ...
    mean_bridges_per_exposed_component;

p_break_from_bridges_linear = ...
    mean_bridges_per_exposed_component * ...
    q_break_given_bridge_global;

fprintf('\n');
fprintf('====================================================================\n');
fprintf(' Ponts exposes dans les graphes union\n');
fprintf('====================================================================\n');
fprintf('Composantes union exposees au total    : %d\n', ...
    total_union_components);
fprintf('Ponts union exposes au total           : %d\n', ...
    total_union_bridges);
fprintf('Ponts union retires au total           : %d\n', ...
    total_removed_union_bridges);
fprintf('Ponts moyens / composante exposee      : %.8f\n', ...
    mean_bridges_per_exposed_component);
fprintf('q_break empirique conditionnel au pont : %.8f\n', ...
    q_break_given_bridge_global);
fprintf('p_break empirique reconstruit probabiliste       : %.8f\n', ...
    p_break_from_bridges);
fprintf('p_break empirique reconstruit lineaire           : %.8f\n', ...
    p_break_from_bridges_linear);
fprintf('====================================================================\n');

figure;
plot(time_values(1:end-1), ...
    mean_bridges_per_union_component_t, ...
    'LineWidth',1.4);
hold on;
yline(mean_bridges_per_exposed_component,'--', ...
    sprintf('moyenne globale = %.4f', ...
    mean_bridges_per_exposed_component), ...
    'LineWidth',1.5);
grid on;
xlabel('Temps (s)');
ylabel('Ponts / composante union');
title('Nombre moyen de ponts par composante exposee');
hold off;

figure;
plot(time_values(1:end-1), ...
    q_break_given_bridge_t, ...
    'LineWidth',1.4);
hold on;
yline(q_break_given_bridge_global,'--', ...
    sprintf('moyenne globale = %.4f', ...
    q_break_given_bridge_global), ...
    'LineWidth',1.5);
grid on;
xlabel('Temps (s)');
ylabel('P(retrait | pont de U_k)');
title('Probabilite conditionnelle de rupture d''un pont');
hold off;


%% ============================================================
%  4.c DECOMPOSITION EXACTE RUPTURES / FUSIONS
%
%  Pour chaque transition G_k -> G_{k+1}, on sépare :
%
%    A_removed = arêtes présentes dans G_k et absentes de G_{k+1}
%    A_added   = arêtes absentes de G_k et présentes dans G_{k+1}
%
%  On construit ensuite le graphe intermédiaire
%
%    G_k_minus = G_k \ A_removed
%
%  qui contient uniquement l'effet des suppressions.
%
%  La fragmentation brute due aux ruptures vaut :
%
%    delta_beta0_break(k)
%      = beta0(G_k_minus) - beta0(G_k) >= 0
%
%  Puis l'effet des créations de liens vaut :
%
%    delta_beta0_merge(k)
%      = beta0(G_{k+1}) - beta0(G_k_minus) <= 0
%
%  On vérifie ainsi exactement :
%
%    beta0(k+1)-beta0(k)
%      = delta_beta0_break(k) + delta_beta0_merge(k)
%
%  On compare enfin delta_beta0_break au nombre de ponts de G_k
%  effectivement supprimés.
%% ============================================================

delta_beta0_net = diff(beta0);

beta0_after_removals = zeros(Nt-1,1);
delta_beta0_break = zeros(Nt-1,1);
delta_beta0_merge = zeros(Nt-1,1);

n_removed_edges = zeros(Nt-1,1);
n_added_edges = zeros(Nt-1,1);

n_bridges_current = zeros(Nt-1,1);
n_removed_bridges_current = zeros(Nt-1,1);

for k = 1:Nt-1

    A_current = logical(Adjacency{k});
    A_next = logical(Adjacency{k+1});

    % Arêtes supprimées et ajoutées pendant la transition
    A_removed = A_current & ~A_next;
    A_added = ~A_current & A_next;

    n_removed_edges(k) = nnz(triu(A_removed,1));
    n_added_edges(k) = nnz(triu(A_added,1));

    % Graphe intermédiaire ne contenant que l'effet des suppressions
    A_after_removals = A_current & ~A_removed;

    G_after_removals = graph(sparse(A_after_removals));
    comp_after_removals = conncomp(G_after_removals);
    beta0_after_removals(k) = max(comp_after_removals);

    % Effet brut des seules ruptures
    delta_beta0_break(k) = ...
        beta0_after_removals(k) - beta0(k);

    % Effet des créations de liens appliquées ensuite
    delta_beta0_merge(k) = ...
        beta0(k+1) - beta0_after_removals(k);

    % Ponts présents dans G_k
    bridge_pairs_current = find_bridges_tarjan(A_current);
    n_bridges_current(k) = size(bridge_pairs_current,1);

    if n_bridges_current(k) > 0
        bridge_idx_current = sub2ind([N N], ...
            bridge_pairs_current(:,1), ...
            bridge_pairs_current(:,2));

        removed_current_mask = A_removed(bridge_idx_current);

        n_removed_bridges_current(k) = ...
            sum(removed_current_mask);
    end
end

% Vérification numérique de la décomposition
decomposition_error = ...
    delta_beta0_net - ...
    (delta_beta0_break + delta_beta0_merge);

max_decomposition_error = max(abs(decomposition_error));

% Agrégats globaux
total_delta_beta0_break = sum(delta_beta0_break);
total_delta_beta0_merge = sum(delta_beta0_merge);
total_removed_bridges_current = sum(n_removed_bridges_current);

% Rapport entre fragmentation brute et ponts initiaux supprimés
gamma_frag_break = ...
    total_delta_beta0_break / ...
    max(total_removed_bridges_current,1);

% Rapport temporel
gamma_frag_break_t = NaN(Nt-1,1);
mask_removed_bridge = n_removed_bridges_current > 0;

gamma_frag_break_t(mask_removed_bridge) = ...
    delta_beta0_break(mask_removed_bridge) ./ ...
    n_removed_bridges_current(mask_removed_bridge);

fprintf('\n');
fprintf('====================================================================\n');
fprintf(' Décomposition exacte ruptures / fusions\n');
fprintf('====================================================================\n');
fprintf('Fragmentation brute totale due aux ruptures : %d\n', ...
    total_delta_beta0_break);
fprintf('Effet total des créations de liens          : %d\n', ...
    total_delta_beta0_merge);
fprintf('Variation nette totale de beta0             : %d\n', ...
    sum(delta_beta0_net));
fprintf('Ponts de G_k effectivement supprimés        : %d\n', ...
    total_removed_bridges_current);
fprintf('gamma_frag_break                            : %.8f\n', ...
    gamma_frag_break);
fprintf('Erreur max de décomposition                 : %.3e\n', ...
    max_decomposition_error);
fprintf('====================================================================\n');

%% Graphe 1 : décomposition de la variation de beta0
figure;
stairs(time_values(1:end-1), ...
    delta_beta0_net, ...
    'LineWidth',1.4);
hold on;
stairs(time_values(1:end-1), ...
    delta_beta0_break, ...
    'LineWidth',1.4);
stairs(time_values(1:end-1), ...
    delta_beta0_merge, ...
    'LineWidth',1.4);
grid on;
xlabel('Temps (s)');
ylabel('Variation de \beta_0');
title('Décomposition de \Delta\beta_0 : ruptures et créations de liens');
legend('\Delta\beta_0 net', ...
       '\Delta\beta_0 dû aux ruptures', ...
       '\Delta\beta_0 dû aux créations', ...
       'Location','best');
hold off;

%% Graphe 2 : fragmentation brute et ponts supprimés
figure;
stairs(time_values(1:end-1), ...
    delta_beta0_break, ...
    'LineWidth',1.5);
hold on;
stairs(time_values(1:end-1), ...
    n_removed_bridges_current, ...
    'LineWidth',1.5);
grid on;
xlabel('Temps (s)');
ylabel('Nombre par transition');
title('Fragmentation brute et ponts de G_k supprimés');
legend('\Delta\beta_0 dû aux seules ruptures', ...
       'Ponts de G_k supprimés', ...
       'Location','best');
hold off;

%% Graphe 3 : rapport fragmentation brute / ponts supprimés
figure;
plot(time_values(1:end-1), ...
    gamma_frag_break_t, ...
    'LineWidth',1.4);
hold on;
yline(gamma_frag_break,'--', ...
    sprintf('\\gamma_{frag,rupt} global = %.4f', ...
    gamma_frag_break), ...
    'LineWidth',1.5);
yline(1,':','Valeur de référence = 1', ...
    'LineWidth',1.2);
grid on;
xlabel('Temps (s)');
ylabel('\gamma_{frag,rupt}(t)');
title('Fragmentation brute rapportée aux ponts initiaux supprimés');
hold off;

%% ============================================================
%  5. SAUVEGARDE DES DONNÉES
%% ============================================================

save('analysis_temp_results.mat', ...
    'N', 'R', 'h', 'lambda', 'dmax', 'dt', 'Tmax', ...
    'time_values', ...
    'positions0', 'r0', 'v', ...
    'Positions', 'Adjacency', ...
    'beta0', 'beta1_graph', 'largest_component', 'num_edges', ...
    'component_size_counts', ...
    'mean_component_count_by_size', ...
    'mean_nodes_by_component_size', ...
    's_values_component', ...
    'beta0_from_component_sizes', ...
    'N_from_component_sizes', ...
    'max_component_count_error', ...
    'max_component_mass_error', ...
    'ZigzagAdjacency', 'ZigzagLabels', ...
    'beta0_zigzag', 'beta1_zigzag_graph', ...
    'largest_component_zigzag', 'num_edges_zigzag', ...
    'n_bridges_union', 'n_removed_bridges_union', ...
    'beta0_union', ...
    'mean_bridges_per_union_component_t', ...
    'q_break_given_bridge_t', ...
    'mean_bridges_per_exposed_component', ...
    'q_break_given_bridge_global', ...
    'p_break_from_bridges', ...
    'p_break_from_bridges_linear', ...
    'delta_beta0_net', ...
    'beta0_after_removals', ...
    'delta_beta0_break', ...
    'delta_beta0_merge', ...
    'n_removed_edges', ...
    'n_added_edges', ...
    'n_bridges_current', ...
    'n_removed_bridges_current', ...
    'total_delta_beta0_break', ...
    'total_delta_beta0_merge', ...
    'total_removed_bridges_current', ...
    'gamma_frag_break_t', ...
    'gamma_frag_break', ...
    'decomposition_error', ...
    'max_decomposition_error');

fprintf('\nAnalyse terminée.\n');
fprintf('Résultats sauvegardés dans analysis_temp_results.mat\n');


%% ============================================================
%  FONCTION LOCALE : PONTS EXACTS PAR TARJAN
%% ============================================================

function bridge_pairs = find_bridges_tarjan(A)
% Retourne les ponts exacts d'un graphe simple non oriente.

    A = logical(A | A.');
    n = size(A,1);
    A(1:n+1:end) = false;

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

        for idx_n = 1:numel(neighbors)
            w = neighbors(idx_n);

            if ~visited(w)
                parent(w) = u;
                dfs(w);

                low(u) = min(low(u),low(w));

                if low(w) > discovery(u)
                    bridge_pairs(end+1,:) = sort([u w]); %#ok<AGROW>
                end

            elseif w ~= parent(u)
                low(u) = min(low(u),discovery(w));
            end
        end
    end
end
