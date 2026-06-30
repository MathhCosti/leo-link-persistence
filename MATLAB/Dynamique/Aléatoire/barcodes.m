clc; close all;

%% ============================================================
%  BARCODE ZIGZAG H0 POUR RÉSEAU LEO
%  Mouvement avec vecteurs aléatoires tangentiels
%
%  Entrée :
%  - leo_zigzag_analysis_random_vectors_results.mat
%
%  Sortie :
%  - intervalles de persistance H0
%  - figure barcode
%  - histogramme des durées de vie
%
%  Suite utilisée :
%  G1 -> G1 U G2 <- G2 -> G2 U G3 <- G3 ...
%% ============================================================

load('leo_zigzag_analysis_random_vectors_results.mat', ...
    'ZigzagAdjacency', 'ZigzagLabels', 'time_values', ...
    'N', 'R', 'lambda', 'dmax', 'dt');

Nz = length(ZigzagAdjacency);

fprintf('Nombre d''objets dans le zigzag : %d\n', Nz);

%% Conversion des labels zigzag en temps physiques
% Exemple :
% label 1   -> t1
% label 1.5 -> milieu entre t1 et t2
% label 2   -> t2

ZigzagTime = zeros(Nz,1);

for k = 1:Nz
    lab = ZigzagLabels(k);

    if abs(lab - round(lab)) < 1e-12
        idx = round(lab);
        ZigzagTime(k) = time_values(idx);
    else
        idx = floor(lab);
        ZigzagTime(k) = 0.5 * (time_values(idx) + time_values(idx+1));
    end
end

%% ============================================================
%  1. Calcul des espaces H0
%
%  H0(G) est de dimension = nombre de composantes connexes.
%  Une base de H0 est donnée par les composantes connexes.
%% ============================================================

component_labels = cell(Nz,1);
h0_dims = zeros(Nz,1);

for k = 1:Nz
    A = ZigzagAdjacency{k};
    G = graph(A);

    comp = conncomp(G);
    comp = comp(:);

    component_labels{k} = comp;
    h0_dims(k) = max(comp);
end

fprintf('Dimensions H0 min/max : %d / %d\n', min(h0_dims), max(h0_dims));

%% ============================================================
%  2. Construction des matrices du module zigzag H0
%
%  Si k impair :
%     G_i -> G_i U G_{i+1}
%
%  Si k pair :
%     G_i U G_{i+1} <- G_{i+1}
%
%  Dans les deux cas, la matrice envoie les composantes du graphe
%  inclus vers les composantes du graphe union.
%% ============================================================

maps = cell(Nz-1,1);

for k = 1:Nz-1

    if mod(k,2) == 1
        % Flèche vers la droite : V_k -> V_{k+1}
        maps{k}.type = 'f';

        maps{k}.mat = build_H0_map( ...
            component_labels{k}, ...
            component_labels{k+1}, ...
            h0_dims(k), ...
            h0_dims(k+1));

    else
        % Flèche vers la gauche : V_k <- V_{k+1}
        % La matrice représente l'application V_{k+1} -> V_k
        maps{k}.type = 'g';

        maps{k}.mat = build_H0_map( ...
            component_labels{k+1}, ...
            component_labels{k}, ...
            h0_dims(k+1), ...
            h0_dims(k));
    end
end

%% ============================================================
%  3. Calcul du barcode zigzag H0
%% ============================================================

intervals = zigzag_barcode_from_module_mod2(h0_dims, maps);

fprintf('Nombre total de barres H0 : %d\n', size(intervals,1));

birth_index = intervals(:,1);
death_index = intervals(:,2);

birth_time = ZigzagTime(birth_index);
death_time = ZigzagTime(death_index);

lifetimes = death_time - birth_time;

%% ============================================================
%  Suppression de la composante H0 persistante globale
%  Elle correspond normalement à la barre née au début et morte à la fin
%% ============================================================

persistent_idx = (birth_index == 1) & (death_index == Nz);

fprintf('Nombre de barres persistantes globales supprimées : %d\n', sum(persistent_idx));

% On garde uniquement les autres barres
keep_idx = ~persistent_idx;

intervals   = intervals(keep_idx, :);
birth_index = birth_index(keep_idx);
death_index = death_index(keep_idx);

birth_time  = birth_time(keep_idx);
death_time  = death_time(keep_idx);
lifetimes   = lifetimes(keep_idx);

%% Sauvegarde
save('leo_H0_zigzag_random_vectors_barcodes.mat', ...
    'intervals', 'birth_index', 'death_index', ...
    'birth_time', 'death_time', 'lifetimes', ...
    'ZigzagTime', 'ZigzagLabels', 'h0_dims');

fprintf('Barcodes sauvegardés dans leo_H0_zigzag_random_vectors_barcodes.mat\n');

%% ============================================================
%  TEST DE QUEUE EXPONENTIELLE DES DUREES DE BARRES H0
%% ============================================================

positive_lifetimes = lifetimes(lifetimes > 0);

Lvals = unique(sort(positive_lifetimes));
survival = zeros(size(Lvals));

for i = 1:length(Lvals)
    survival(i) = mean(positive_lifetimes >= Lvals(i));
end

% Tracé de la fonction de survie empirique
figure;
semilogy(Lvals, survival, 'o-', 'LineWidth', 1.5);
grid on;
hold on;

% On recalcule omega pour éviter de dépendre d'une variable non sauvée.
mu = 398600;                    % km^3/s^2
omega = sqrt(mu / R^3);         % rad/s

% Dans le modèle à vecteurs tangentiels aléatoires, tous les satellites
% ont la même norme de vitesse v_orb, mais des directions indépendantes.
% Pour deux directions indépendantes dans un plan tangent local,
% une approximation classique de la vitesse relative moyenne est :
% E[|v1 - v2|] ≈ 4/pi * v_orb.
% L'ancienne valeur 4/3 est proche mais venait plutôt d'une approximation grossière.
v_orb = R * omega;              % km/s
v_rel = (4/pi) * v_orb;         % approximation moyenne avec directions aléatoires

% Aire balayée pendant un pas de temps
A_sweep = 2 * dmax * v_rel * dt;

% Probabilité théorique de fusion pendant un pas
p_merge = 1 - exp(-lambda * A_sweep);

% Probabilité théorique de rupture d'un lien individuel pendant un pas.
% Approximation : seuls les liens proches de la frontière dmax peuvent casser.
% q_break ≈ (2/pi) * (v_rel * dt / dmax)
q_break_raw = (2/pi) * (v_rel * dt / dmax);

% On borne la probabilité pour éviter des valeurs non physiques si dt est trop grand.
% Si q_break_raw devient proche de 1, l'approximation linéaire n'est plus fiable.
q_break = min(max(q_break_raw, 0), 1 - eps);

% Probabilité totale de disparition/modification pendant un pas :
% p_death = 1 - P(pas de fusion et pas de rupture)
p_death = 1 - (1 - p_merge) * (1 - q_break);

% Temps caractéristique théorique
tau_th = -dt / log(1 - p_death);

% Courbe de survie théorique
survival_th = exp(-Lvals / tau_th);

semilogy(Lvals, survival_th, '--', 'LineWidth', 2);

xlabel('Durée des barres (s)');
ylabel('Probabilité de survie');
title('Survie des barres H0 — vecteurs aléatoires');

legend('Données simulées', 'Modèle exponentiel approx. (merge + break)', 'Location', 'best');

hold off;

fprintf('\n--- Analyse des durées de barres H0 ---\n');
fprintf('Temps caractéristique théorique tau_th : %.2f s\n', tau_th);
fprintf('Vitesse orbitale v_orb : %.3f km/s\n', v_orb);
fprintf('Vitesse relative moyenne approx. v_rel : %.3f km/s\n', v_rel);
fprintf('Probabilité de fusion p_merge : %.6f\n', p_merge);
fprintf('Probabilité de rupture q_break : %.6f\n', q_break);
fprintf('Probabilité totale p_death : %.6f\n', p_death);
if q_break_raw >= 1
    fprintf('Attention : q_break_raw = %.3f >= 1, approximation linéaire invalide pour ce dt.\n', q_break_raw);
end

%% ============================================================
%  4. Affichage du barcode
%% ============================================================

% Tri par durée décroissante
[~, order] = sort(lifetimes, 'descend');

% Pour éviter une figure illisible si beaucoup de barres
maxBarsToPlot = 150;
order = order(1:min(maxBarsToPlot, length(order)));

figure;
hold on;
grid on;

for ii = 1:length(order)
    id = order(ii);

    x0 = birth_time(id);
    x1 = death_time(id);
    y = ii;

    if abs(x1 - x0) < 1e-12
        plot(x0, y, 'ko', 'MarkerSize', 4);
    else
        plot([x0 x1], [y y], 'k-', 'LineWidth', 1.2);
    end
end

xlabel('Temps (s)');
ylabel('Barres H_0 triées par durée décroissante');
title(sprintf('Barcode zigzag H_0 — vecteurs aléatoires — %d plus longues barres', length(order)));

hold off;

%% ============================================================
%  5. Histogramme des durées
%% ============================================================

figure;
histogram(lifetimes, 30);
grid on;
xlabel('Durée de vie des composantes (s)');
ylabel('Nombre de barres');
title('Distribution des durées de vie des composantes H_0 — vecteurs aléatoires');

%% ============================================================
%  6. Quelques statistiques
%% ============================================================

fprintf('\nStatistiques sur les barres H0 :\n');
fprintf('Durée moyenne : %.2f s\n', mean(lifetimes));
fprintf('Durée médiane : %.2f s\n', median(lifetimes));
fprintf('Durée max : %.2f s\n', max(lifetimes));
fprintf('Nombre de barres de durée nulle : %d\n', sum(lifetimes == 0));

long_threshold = 0.2 * max(ZigzagTime);
fprintf('Nombre de barres longues (> %.1f s) : %d\n', ...
    long_threshold, sum(lifetimes > long_threshold));

%% ============================================================
%  FONCTIONS LOCALES
%% ============================================================

function M = build_H0_map(labels_source, labels_target, dim_source, dim_target)
    % Construit la matrice induite en H0 par une inclusion de graphes.
    %
    % Chaque composante source est envoyée vers la composante cible
    % qui la contient.
    %
    % M est de taille dim_target x dim_source.

    M = zeros(dim_target, dim_source);

    for c = 1:dim_source
        vertices = find(labels_source == c);

        target_comps = unique(labels_target(vertices));

        if length(target_comps) ~= 1
            error(['Inclusion invalide pour H0 : une composante source ', ...
                   'est envoyée dans plusieurs composantes cibles.']);
        end

        target_c = target_comps(1);
        M(target_c, c) = 1;
    end

    M = mod(M,2);
end

function intervals = zigzag_barcode_from_module_mod2(dims, maps)
    % Calcule le barcode zigzag d'un module de type quelconque
    % sur F2, à partir des dimensions et des matrices.
    %
    % Sortie :
    % intervals : matrice nb_intervalles x 2
    %             chaque ligne est [birth_index, death_index]

    n = length(dims);

    % Right-filtration initiale : R = (0, V1)
    R = cell(2,1);
    R{1} = zeros(dims(1),0);
    R{2} = eye(dims(1));

    % Liste des temps de naissance associés aux quotients
    b = 1;

    % Dimensions des sous-quotients
    r = filtration_quotient_dims(R);

    intervals = [];

    for k = 1:n-1

        current_type = maps{k}.type;

        if current_type == 'f'
            % V_k -> V_{k+1}
            M = maps{k}.mat;

            Rnext = cell(length(R)+1,1);

            for i = 1:length(R)
                Rnext{i} = gf2_col_basis(M * R{i});
            end

            Rnext{end} = eye(dims(k+1));

            bnext = [b, k+1];

            rnext = filtration_quotient_dims(Rnext);

            % Features qui meurent à k :
            % c_i^k = r_i^k - r_i^{k+1}
            for i = 1:length(r)
                c = r(i) - rnext(i);

                if c > 0
                    intervals = [intervals; repmat([b(i), k], c, 1)];
                end
            end

        elseif current_type == 'g'
            % V_k <- V_{k+1}
            % La matrice N représente V_{k+1} -> V_k
            N = maps{k}.mat;

            Rnext = cell(length(R)+1,1);
            Rnext{1} = zeros(dims(k+1),0);

            for i = 1:length(R)
                Rnext{i+1} = gf2_preimage(N, R{i});
            end

            bnext = [k+1, b];

            rnext = filtration_quotient_dims(Rnext);

            % Features qui meurent à k :
            % c_i^k = r_i^k - r_{i+1}^{k+1}
            for i = 1:length(r)
                c = r(i) - rnext(i+1);

                if c > 0
                    intervals = [intervals; repmat([b(i), k], c, 1)];
                end
            end

        else
            error('Type de flèche inconnu.');
        end

        R = Rnext;
        b = bnext;
        r = rnext;
    end

    % À la fin, toutes les features encore vivantes meurent à n
    for i = 1:length(r)
        c = r(i);

        if c > 0
            intervals = [intervals; repmat([b(i), n], c, 1)];
        end
    end
end

function dims = filtration_quotient_dims(R)
    % Calcule dim(R_i / R_{i-1}) pour une filtration
    % R = {R0, R1, ..., Rn}

    m = length(R) - 1;
    dims = zeros(1,m);

    for i = 1:m
        dims(i) = gf2_rank(R{i+1}) - gf2_rank(R{i});
    end
end

function P = gf2_preimage(A, S)
    % Calcule la préimage A^{-1}(S) sur F2.
    %
    % A : matrice m x n
    % S : base d'un sous-espace de F2^m, taille m x s
    %
    % On cherche x tel que A x appartient à Span(S).
    % Donc il existe y tel que A x + S y = 0.
    %
    % On calcule le noyau de [A S], puis on projette sur les coordonnées x.

    A = mod(full(A),2);
    S = mod(full(S),2);

    n = size(A,2);

    Big = [A S];
    Z = gf2_null(Big);

    X = Z(1:n,:);

    P = gf2_col_basis(X);
end

function B = gf2_col_basis(A)
    % Extrait une base de colonnes indépendantes de A sur F2.

    A = mod(full(A),2);

    if isempty(A)
        B = zeros(size(A,1),0);
        return;
    end

    [~, pivots] = gf2_rref(A);

    if isempty(pivots)
        B = zeros(size(A,1),0);
    else
        B = mod(A(:,pivots),2);
    end
end

function r = gf2_rank(A)
    % Rang sur F2.

    A = mod(full(A),2);

    if isempty(A)
        r = 0;
        return;
    end

    [~, pivots] = gf2_rref(A);
    r = length(pivots);
end

function Z = gf2_null(A)
    % Base du noyau de A sur F2.
    %
    % A est m x n.
    % Z est n x d, les colonnes forment une base de Ker(A).

    A = mod(full(A),2);
    [R, pivots] = gf2_rref(A);

    n = size(A,2);
    free_cols = setdiff(1:n, pivots);

    if isempty(free_cols)
        Z = zeros(n,0);
        return;
    end

    Z = zeros(n, length(free_cols));

    for j = 1:length(free_cols)
        f = free_cols(j);

        z = zeros(n,1);
        z(f) = 1;

        for p = 1:length(pivots)
            col = pivots(p);
            z(col) = R(p,f);
        end

        Z(:,j) = mod(z,2);
    end
end

function [R, pivots] = gf2_rref(A)
    % Forme échelonnée réduite sur F2.
    %
    % Retourne aussi les colonnes pivot.

    A = mod(full(A),2);
    [m,n] = size(A);

    R = A;
    pivots = [];

    row = 1;

    for col = 1:n
        if row > m
            break;
        end

        pivot_rel = find(R(row:m,col), 1);

        if isempty(pivot_rel)
            continue;
        end

        pivot = pivot_rel + row - 1;

        % échange de lignes
        if pivot ~= row
            tmp = R(row,:);
            R(row,:) = R(pivot,:);
            R(pivot,:) = tmp;
        end

        % élimination sur toutes les autres lignes
        for rr = 1:m
            if rr ~= row && R(rr,col) == 1
                R(rr,:) = mod(R(rr,:) + R(row,:), 2);
            end
        end

        pivots(end+1) = col;
        row = row + 1;
    end
end