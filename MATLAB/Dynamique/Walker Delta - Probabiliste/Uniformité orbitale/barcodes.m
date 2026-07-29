clc; close all;

%% ============================================================
%  BARCODE ZIGZAG H0 POUR RÉSEAU LEO - WALKER-DELTA
%
%  Entrée :
%  - leo_zigzag_analysis_results_delta.mat
%
%  Sortie :
%  - intervalles de persistance H0
%  - figure barcode
%  - histogramme des durées de vie
%
%  Suite utilisée :
%  G1 -> G1 U G2 <- G2 -> G2 U G3 <- G3 ...
%% ============================================================

load('analysis_temp_results.mat', ...
    'ZigzagAdjacency', 'ZigzagLabels', 'time_values', 'N', ...
    'R', 'lambda', 'dmax', 'dt', 'omega', 'inc_deg', 'P');

% Securite si un ancien fichier de resultats ne contient pas omega.
if ~exist('omega', 'var')
    mu = 398600;
    omega = sqrt(mu / R^3);
end
if ~exist('inc_deg', 'var')
    inc_deg = NaN;
end
if ~exist('P', 'var')
    P = NaN;
end

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

%% Chargement des probabilites theoriques p_merge et p_break
%
% Les deux codes theoriques doivent avoir ete executes auparavant :
%   - pmerge_th.m
%   - pbreak_th.m
%
% Ils produisent respectivement :
%   pmerge_theorique_walker_delta_results.mat
%   pbreak_theorique_walker_delta_results.mat

script_dir = fileparts(mfilename('fullpath'));

pmerge_results_file = fullfile(script_dir, 'Probabilité fusion', ...
    'pmerge_th_results.mat');

pbreak_results_file = fullfile(script_dir, 'Probabilité rupture', ...
    'pbreak_th_results.mat');

if ~isfile(pmerge_results_file)
    error(['Fichier introuvable : %s\n' ...
        'Executer d''abord pmerge_th.m.'],pmerge_results_file);
end

if ~isfile(pbreak_results_file)
    error(['Fichier introuvable : %s\n' ...
        'Executer d''abord pbreak_th.m.'],pbreak_results_file);
end

Smerge = load(pmerge_results_file);
Sbreak = load(pbreak_results_file);

%% Recuperation de p_merge theorique
% On retient la moyenne orbitale obtenue par integration sur u.
if isfield(Smerge,'p_merge_mean')
    p_merge_th = Smerge.p_merge_mean;
elseif isfield(Smerge,'p_merge_delta')
    p_merge_th = Smerge.p_merge_delta;
else
    error(['Aucune variable p_merge_mean ou p_merge_delta ' ...
        'dans %s.'],pmerge_results_file);
end

%% Recuperation de p_break theorique
if isfield(Sbreak,'p_break_delta')
    p_break_th = Sbreak.p_break_delta;
else
    error('La variable p_break_delta est absente de %s.', ...
        pbreak_results_file);
end

%% Verification de coherence du pas temporel
if isfield(Smerge,'Delta_t') && abs(Smerge.Delta_t-dt) > 1e-12
    warning(['Le pas temporel du calcul p_merge (%.4f s) ' ...
        'differe de celui des barcodes (%.4f s).'], ...
        Smerge.Delta_t,dt);
end

if isfield(Sbreak,'Delta_t') && abs(Sbreak.Delta_t-dt) > 1e-12
    warning(['Le pas temporel du calcul p_break (%.4f s) ' ...
        'differe de celui des barcodes (%.4f s).'], ...
        Sbreak.Delta_t,dt);
end

%% Probabilite theorique totale de disparition
% Sous l'hypothese d'independance conditionnelle :
%   p_disp = 1-(1-p_merge)(1-p_break).
p_disp_th = 1-(1-p_merge_th)*(1-p_break_th);
p_disp_th = p_disp_th ./ 2;
p_disp_th = min(max(p_disp_th,0),1);

% Temps caracteristique theorique associe a une loi geometrique
% par pas de temps, puis exponentielle en temps continu.
if p_disp_th <= 0
    tau_th = Inf;
    survival_th = ones(size(Lvals));
elseif p_disp_th >= 1
    tau_th = 0;
    survival_th = double(Lvals <= 0);
else
    tau_th = -dt/log(1-p_disp_th);
    survival_th = exp(-Lvals/tau_th);
end

semilogy(Lvals, survival_th, '--', 'LineWidth', 2);

xlabel('Durée des barres (s)');
ylabel('Probabilité de survie');
title(sprintf('Survie des barres H0 - Walker-Delta, i = %.1f deg', inc_deg));

legend('Données simulées', 'Modèle exponentiel', 'Location', 'best');

hold off;

fprintf('\n--- Analyse theorique des durees de barres H0 ---\n');
fprintf('p_merge theorique                     : %.10f\n',p_merge_th);
fprintf('p_break theorique                     : %.10f\n',p_break_th);
fprintf('p_disp theorique sous independance    : %.10f\n',p_disp_th);
fprintf('Temps caracteristique theorique       : %.2f s\n',tau_th);

%% Sauvegarde
save('barcodes_results.mat', ...
    'intervals', 'birth_index', 'death_index', ...
    'birth_time', 'death_time', 'lifetimes', ...
    'ZigzagTime', 'ZigzagLabels', 'h0_dims',...
    'p_merge_th','p_break_th','p_disp_th','tau_th','survival_th');

fprintf('Barcodes sauvegardes dans barcodes_results.mat\n');

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
title(sprintf('Barcode zigzag H_0 Walker-Delta — %d plus longues barres', length(order)));

hold off;

%% ============================================================
%  5. Histogramme des durées
%% ============================================================

figure;
histogram(lifetimes, 30);
grid on;
xlabel('Durée de vie des composantes (s)');
ylabel('Nombre de barres');
title('Distribution des durees de vie des composantes H_0 - Walker-Delta');

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