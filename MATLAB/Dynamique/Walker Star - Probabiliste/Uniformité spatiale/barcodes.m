clc; close all;

%% ============================================================
%  BARCODE ZIGZAG H0 POUR RÉSEAU LEO
%
%  Entrée :
%  - leo_zigzag_analysis_results.mat
%
%  Sortie :
%  - intervalles de persistance H0
%  - figure barcode
%  - histogramme des durées de vie
%
%  Suite utilisée :
%  G1 -> G1 U G2 <- G2 -> G2 U G3 <- G3 ...
%% ============================================================

%% Parametres physiques
R_earth = 6371;      % km
h = 550;             % km
R = R_earth + h;     % rayon orbital

mu = 398600;              % km^3/s^2
omega = sqrt(mu / R^3);   % vitesse angulaire orbitale rad/s

T_half = pi / omega;

load('leo_zigzag_analysis_results.mat', ...
    'ZigzagAdjacency', 'ZigzagLabels', 'time_values', 'N');

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

%% Sauvegarde
save('leo_H0_zigzag_barcodes.mat', ...
    'intervals', 'birth_index', 'death_index', ...
    'birth_time', 'death_time', 'lifetimes', ...
    'ZigzagTime', 'ZigzagLabels', 'h0_dims');

fprintf('Barcodes sauvegardés dans leo_H0_zigzag_barcodes.mat\n');

%% ============================================================
%  TEST DE QUEUE EXPONENTIELLE DES DUREES DE BARRES H0
%  avec p_disp théorique temporel
%% ============================================================

positive_lifetimes_all = lifetimes(lifetimes > 0);

%% Chargement du p_disp théorique temporel
% Le fichier doit contenir au minimum :
%   - p_th_v : p_disp^th(t)
%   - tv ou t_uniform : instants associés
script_dir = fileparts(mfilename('fullpath'));
pdisp_file = fullfile(script_dir, 'Probabilité disparition', 'comparaison_pdisp_th_emp_spectres_results.mat');

% Repli si le fichier est dans le dossier courant
if ~isfile(pdisp_file)
    pdisp_file_alt = fullfile(script_dir, 'comparaison_pdisp_th_emp_spectres_results.mat');
    if isfile(pdisp_file_alt)
        pdisp_file = pdisp_file_alt;
    end
end

if ~isfile(pdisp_file)
    error(['Fichier introuvable : ', pdisp_file, newline, ...
           'Lance d''abord le script qui calcule p_{disp}^{th}(t), ', ...
           'ou place comparaison_pdisp_th_emp_spectres_results.mat dans le dossier courant.']);
end

Sdisp = load(pdisp_file);

if isfield(Sdisp, 'p_th_v')
    pdisp_th_t = Sdisp.p_emp_v(:);
else
    error('Le fichier %s ne contient pas la variable p_th_v.', pdisp_file);
end

if isfield(Sdisp, 'tv')
    t_pdisp = Sdisp.tv(:);
elseif isfield(Sdisp, 't_uniform')
    t_pdisp = Sdisp.t_uniform(:);
else
    % Repli : on utilise le pas de temps des graphes zigzag
    dt_tmp = median(diff(time_values));
    t_pdisp = (0:length(pdisp_th_t)-1)' * dt_tmp;
    warning('Aucun vecteur temps trouvé pour p_th_v : temps reconstruit avec time_values.');
end

% Nettoyage et mise au bon format
valid = isfinite(t_pdisp) & isfinite(pdisp_th_t);
t_pdisp = t_pdisp(valid);
pdisp_th_t = pdisp_th_t(valid);

[t_pdisp, order_p] = sort(t_pdisp);
pdisp_th_t = pdisp_th_t(order_p);

% On force la probabilité dans [0,1[
pdisp_th_t = min(max(pdisp_th_t, 0), 1 - eps);

% Pas de temps associé au p_disp temporel
dt_pdisp = median(diff(t_pdisp));

%% Restriction explicite à une demi-révolution
% On ne compare la survie que pour des durées L inférieures à une
% demi-révolution. Les tranches quasi stationnaires sont également
% construites uniquement sur cet intervalle [0, T_half].
%
% Par défaut, on utilise la période associée à la fréquence dominante de
% p_disp^th(t), si elle est disponible dans le fichier de résultats. Cette
% période correspond à la périodicité demi-orbitale observée dans les courbes.

% Sécurité : on ne dépasse pas la fenêtre effectivement couverte par p_disp(t)
T_pdisp_available = max(t_pdisp) - min(t_pdisp);
T_half = min(T_half, T_pdisp_available);

% Points de survie uniquement jusqu'à la demi-révolution.
% Important : la probabilité de survie reste calculée par rapport à toutes
% les barres positives, mais on n'affiche/évalue que les durées L <= T_half.
Lvals = unique(sort(positive_lifetimes_all(positive_lifetimes_all <= T_half)));

if isempty(Lvals)
    error('Aucune durée de barre positive inférieure à T_half = %.3f s.', T_half);
end

survival = zeros(size(Lvals));
for i = 1:length(Lvals)
    survival(i) = mean(positive_lifetimes_all >= Lvals(i));
end

% Tracé de la fonction de survie empirique restreinte à la demi-révolution
figure;
semilogy(Lvals, survival, 'o-', 'LineWidth', 1.5);
grid on;
hold on;

%% Survie théorique par intervalles quasi stationnaires
% Au lieu d'appliquer directement p_disp^th(t) à chaque pas de temps,
% on découpe l'intervalle de durée étudié en tranches. Sur chaque tranche,
% p_disp^th(t) est remplacé par sa moyenne temporelle. Cela donne une
% survie exponentielle par morceaux, c'est-à-dire des droites par morceaux
% en échelle semi-logarithmique.

nb_intervalles = 12;      % à ajuster si besoin
Lmax_model = T_half;
interval_edges = linspace(0, Lmax_model, nb_intervalles + 1);

% Temps relatif du modèle p_disp(t), utilisé pour moyenner sur les tranches
t_rel_pdisp = t_pdisp - min(t_pdisp);

pdisp_piece_mean = zeros(nb_intervalles,1);

for k = 1:nb_intervalles
    a = interval_edges(k);
    b = interval_edges(k+1);

    if k < nb_intervalles
        mask = (t_rel_pdisp >= a) & (t_rel_pdisp < b);
    else
        mask = (t_rel_pdisp >= a) & (t_rel_pdisp <= b);
    end

    if any(mask)
        pdisp_piece_mean(k) = mean(pdisp_th_t(mask));
    else
        % Si aucun point p_disp ne tombe exactement dans la tranche,
        % on interpole au centre de la tranche.
        mid = 0.5 * (a + b);
        pdisp_piece_mean(k) = interp1(t_rel_pdisp, pdisp_th_t, mid, 'linear', 'extrap');
    end
end

% On force les moyennes dans [0,1[
pdisp_piece_mean = min(max(pdisp_piece_mean, 0), 1 - eps);

survival_th_temp = zeros(size(Lvals));

for i = 1:length(Lvals)
    L = Lvals(i);
    logS = 0;

    for k = 1:nb_intervalles
        a = interval_edges(k);
        b = interval_edges(k+1);

        % Portion de la durée L qui tombe dans la tranche k
        overlap = max(0, min(L, b) - a);

        if overlap > 0
            logS = logS + (overlap / dt_pdisp) * log(1 - pdisp_piece_mean(k));
        end
    end

    survival_th_temp(i) = exp(logS);
end

% Pour visualiser les droites par morceaux, on évalue aussi le modèle aux
% bords des intervalles.
survival_piece_edges = ones(size(interval_edges));
logS_edge = 0;
for k = 1:nb_intervalles
    duration_k = interval_edges(k+1) - interval_edges(k);
    logS_edge = logS_edge + (duration_k / dt_pdisp) * log(1 - pdisp_piece_mean(k));
    survival_piece_edges(k+1) = exp(logS_edge);
end

semilogy(Lvals, survival_th_temp, '--', 'LineWidth', 2);
semilogy(interval_edges, survival_piece_edges, 's--', 'LineWidth', 1.2, 'MarkerSize', 4);

xlabel('Durée des barres (s)');
ylabel('Probabilité de survie');
title('Survie des barres H_0 — approximation par tranches');

legend('Données simulées', 'Modèle par tranches', 'Bords des tranches', 'Location', 'best');

hold off;

% Temps caractéristique équivalent à partir de la moyenne temporelle
pdisp_th_mean = mean(pdisp_th_t);
tau_th_equiv = -dt_pdisp / log(1 - pdisp_th_mean);

fprintf('\n--- Analyse des durées de barres H0 ---\n');
fprintf('Durée moyenne positive empirique (toutes barres) : %.2f s\n', mean(positive_lifetimes_all));
fprintf('Demi-révolution utilisée T_half : %.2f s\n', T_half);
fprintf('Nombre de points de survie conservés (L <= T_half) : %d / %d\n', length(Lvals), length(unique(positive_lifetimes_all)));
fprintf('p_disp^th moyen temporel utilise : %.6f\n', pdisp_th_mean);
fprintf('Nombre de tranches quasi stationnaires : %d\n', nb_intervalles);
fprintf('p_disp moyen par tranche : '); fprintf('%.4f ', pdisp_piece_mean); fprintf('\n');
fprintf('tau équivalent depuis p_disp^th moyen : %.2f s\n', tau_th_equiv);

% Sauvegarde complémentaire des quantités liées au modèle temporel
save('leo_H0_zigzag_barcodes.mat', ...
    'pdisp_th_t', 't_pdisp', 'survival_th_temp', ...
    'pdisp_th_mean', 'tau_th_equiv', ...
    'nb_intervalles', 'interval_edges', 'pdisp_piece_mean', 'survival_piece_edges', ...
    'T_half', 'positive_lifetimes_all', ...
    '-append');

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
title(sprintf('Barcode zigzag H_0 — %d plus longues barres', length(order)));

hold off;

%% ============================================================
%  5. Histogramme des durées
%% ============================================================

figure;
histogram(lifetimes, 30);
grid on;
xlabel('Durée de vie des composantes (s)');
ylabel('Nombre de barres');
title('Distribution des durées de vie des composantes H_0');

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