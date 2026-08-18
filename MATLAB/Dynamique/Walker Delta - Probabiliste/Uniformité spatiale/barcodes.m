clc; close all;

%% ============================================================
%  BARCODE ZIGZAG H0 POUR RÉSEAU LEO - WALKER-DELTA SPATIAL
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
%  SURVIE THEORIQUE A PARTIR DE p_merge ET p_break
%
%  Le modele p_disp n'est plus utilise.
%
%  On charge les series temporelles theoriques :
%
%      p_merge_global_th(t)
%      p_break_global_th(t)
%
%  puis on calcule leurs moyennes temporelles :
%
%      p_merge_mean_th = mean(p_merge_global_th)
%      p_break_mean_th = mean(p_break_global_th)
%
%  La probabilite constante utilisee pour la courbe de survie est
%  definie, conformement au choix retenu ici, comme la moyenne des
%  deux valeurs :
%
%      p_change_mean_th =
%          (p_merge_mean_th + p_break_mean_th)/2.
%
%  La survie theorique est ensuite :
%
%      tau_th = -dt_model/log(1-p_change_mean_th)
%      S_th(L) = exp(-L/tau_th)
%% ============================================================

positive_lifetimes = lifetimes(lifetimes > 0);

Lvals = unique(sort(positive_lifetimes));
survival = zeros(size(Lvals));

for i = 1:length(Lvals)
    survival(i) = mean(positive_lifetimes >= Lvals(i));
end

%% Localisation des deux fichiers theoriques t,phi
script_dir = fileparts(mfilename('fullpath'));

pmerge_file = fullfile( ...
    script_dir, ...
    'Probabilité fusion', ...
    'pmerge_t_phi_th_results.mat');

pbreak_file = fullfile( ...
    script_dir, ...
    'Probabilité rupture', ...
    'pbreak_t_phi_th_results.mat');

if ~isfile(pmerge_file)
    error(['Fichier pmerge_t_phi_th_results.mat introuvable : %s. ', ...
        'Lancez d''abord pmerge_t_phi_th.m.'],pmerge_file);
end

if ~isfile(pbreak_file)
    error(['Fichier pbreak_t_phi_th_results.mat introuvable : %s. ', ...
        'Lancez d''abord pbreak_t_phi_th.m.'],pbreak_file);
end

Smerge = load(pmerge_file);
Sbreak = load(pbreak_file);

% Les scripts t,phi sauvegardent directement les probabilites
% globales ponderees par le nombre local de composantes.
if ~isfield(Smerge,'p_merge_global_th')
    error(['La variable p_merge_global_th est absente de %s. ', ...
        'Regenerer le fichier avec pmerge_t_phi_th.m.'],pmerge_file);
end

if ~isfield(Sbreak,'p_break_global_th')
    error(['La variable p_break_global_th est absente de %s. ', ...
        'Regenerer le fichier avec pbreak_t_phi_th.m.'],pbreak_file);
end

p_merge_th_t = double(Smerge.p_merge_global_th(:));
p_break_th_t = double(Sbreak.p_break_global_th(:));

% On conserve uniquement les valeurs finies pour les moyennes.
p_merge_th_t = p_merge_th_t(isfinite(p_merge_th_t));
p_break_th_t = p_break_th_t(isfinite(p_break_th_t));

if isempty(p_merge_th_t)
    error('La serie p_merge_global_th est vide ou invalide.');
end

if isempty(p_break_th_t)
    error('La serie p_break_global_th est vide ou invalide.');
end

%% Moyennes temporelles
p_merge_mean_th = mean(p_merge_th_t,'omitnan');
p_break_mean_th = mean(p_break_th_t,'omitnan');

% Choix demande : moyenne entre les moyennes de p_merge et p_break.
p_change_mean_th = 0.5*(p_merge_mean_th+p_break_mean_th);
p_change_mean_th = min(max(p_change_mean_th,0),1-eps);

%% Pas temporel du modele
dt_candidates = [];

% Les deux scripts t,phi sauvegardent DeltaT.
if isfield(Smerge,'DeltaT')
    dt_candidates(end+1) = double(Smerge.DeltaT); %#ok<SAGROW>
elseif isfield(Smerge,'time_values') && numel(Smerge.time_values) >= 2
    dt_candidates(end+1) = ...
        median(diff(double(Smerge.time_values(:)))); %#ok<SAGROW>
end

if isfield(Sbreak,'DeltaT')
    dt_candidates(end+1) = double(Sbreak.DeltaT); %#ok<SAGROW>
elseif isfield(Sbreak,'time_values') && numel(Sbreak.time_values) >= 2
    dt_candidates(end+1) = ...
        median(diff(double(Sbreak.time_values(:)))); %#ok<SAGROW>
end

if isempty(dt_candidates)
    dt_model = dt;
else
    dt_model = mean(dt_candidates,'omitnan');

    if numel(dt_candidates) == 2 && ...
            abs(dt_candidates(1)-dt_candidates(2)) > ...
            1e-9*max(1,abs(dt_model))
        warning(['Les pas temporels de p_merge et p_break diffèrent ', ...
            '(%.6g s et %.6g s). Leur moyenne %.6g s est utilisée.'], ...
            dt_candidates(1),dt_candidates(2),dt_model);
    end
end

%% Temps caracteristique et survie theorique
if p_change_mean_th > 0
    tau_th = -dt_model/log(1-p_change_mean_th);
    survival_th = exp(-Lvals/tau_th);
else
    tau_th = Inf;
    survival_th = ones(size(Lvals));
end

%% Trace empirique / theorie
figure;
semilogy(Lvals,survival,'o-', ...
    'LineWidth',1.5, ...
    'DisplayName','Données simulées');

hold on;
grid on;

semilogy(Lvals,survival_th,'--', ...
    'LineWidth',2, ...
    'DisplayName',sprintf( ...
        ['Modèle exponentiel, ', ...
         '\bar{p}=(\bar{p}_{merge}+\bar{p}_{break})/2=%.4g'], ...
        p_change_mean_th));

xlabel('Durée des barres (s)');
ylabel('Probabilité de survie');
title(sprintf( ...
    'Survie des barres H0 - Delta spatial, i = %.1f deg', ...
    inc_deg));

legend('Location','best');
hold off;

fprintf('\n--- Survie theorique avec p_merge et p_break moyens ---\n');
fprintf('Fichier p_merge charge             : %s\n',pmerge_file);
fprintf('Fichier p_break charge             : %s\n',pbreak_file);
fprintf('Nombre de valeurs de p_merge_t     : %d\n',numel(p_merge_th_t));
fprintf('Nombre de valeurs de p_break_t     : %d\n',numel(p_break_th_t));
fprintf('Pas temporel retenu                : %.6f s\n',dt_model);
fprintf('p_merge theorique moyen            : %.8f\n',p_merge_mean_th);
fprintf('p_break theorique moyen            : %.8f\n',p_break_mean_th);
fprintf('Moyenne entre les deux             : %.8f\n',p_change_mean_th);
fprintf('Temps caracteristique theorique    : %.6f s\n',tau_th);

%% Sauvegarde des resultats de survie
save('barcodes_results.mat', ...
    'Lvals','survival','survival_th', ...
    'p_merge_th_t','p_break_th_t', ...
    'p_merge_mean_th','p_break_mean_th','p_change_mean_th', ...
    'dt_model','tau_th','pmerge_file','pbreak_file', ...
    'intervals','birth_index','death_index', ...
    'birth_time','death_time','lifetimes', ...
    'ZigzagTime','ZigzagLabels','h0_dims');

fprintf('Resultats de survie sauvegardes dans barcodes_results.mat\n');

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
title(sprintf('Barcode zigzag H_0 Delta spatial — %d plus longues barres', length(order)));

hold off;

%% ============================================================
%  5. Histogramme des durées
%% ============================================================

figure;
histogram(lifetimes, 30);
grid on;
xlabel('Durée de vie des composantes (s)');
ylabel('Nombre de barres');
title('Distribution des durees de vie des composantes H_0 - Delta spatial');

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

function file = first_existing_file(candidates)
    file = '';

    for kk = 1:numel(candidates)
        if isfile(candidates{kk})
            file = candidates{kk};
            return;
        end
    end
end

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