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

load('leo_zigzag_analysis_results_delta.mat', ...
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

%% Sauvegarde
save('leo_H0_zigzag_barcodes_delta_spatial.mat', ...
    'intervals', 'birth_index', 'death_index', ...
    'birth_time', 'death_time', 'lifetimes', ...
    'ZigzagTime', 'ZigzagLabels', 'h0_dims');

fprintf('Barcodes sauvegardes dans leo_H0_zigzag_barcodes_delta_spatial.mat\n');

%% ============================================================
%  MODELE THEORIQUE DE SURVIE A p_disp(t) VARIABLE
%% ============================================================

positive_mask = lifetimes > 0;
positive_lifetimes = lifetimes(positive_mask);

Lvals = unique(sort(positive_lifetimes));
survival_emp = zeros(size(Lvals));
for ii = 1:numel(Lvals)
    survival_emp(ii) = mean(positive_lifetimes >= Lvals(ii));
end

script_dir = fileparts(mfilename('fullpath'));
pdisp_file = fullfile(script_dir, 'Probabilité disparition', 'pdisp_modele_delta_spatial_results.mat');

if isempty(pdisp_file)
    error(['Fichier pdisp_modele_delta_spatial_results.mat introuvable. ', ...
        'Lance d''abord pdisp_modele_delta_spatial.m.']);
end

Sp = load(pdisp_file);
required_pdisp = {'p_disp_t','p_break_t'};
for qq = 1:numel(required_pdisp)
    if ~isfield(Sp, required_pdisp{qq})
        error('Variable %s absente de %s.', required_pdisp{qq}, pdisp_file);
    end
end

p_disp_th_t = double(Sp.p_disp_t(:));
p_break_th_t = double(Sp.p_break_t(:));

if isfield(Sp, 't_transition')
    t_pdisp_th = double(Sp.t_transition(:));
elseif isfield(Sp, 'time_values')
    t_tmp = double(Sp.time_values(:));
    t_pdisp_th = t_tmp(1:numel(p_disp_th_t));
else
    error('Aucune grille temporelle trouvée dans %s.', pdisp_file);
end

if numel(t_pdisp_th) ~= numel(p_disp_th_t)
    error('Tailles incompatibles entre t_transition et p_disp_t.');
end
if numel(p_break_th_t) ~= numel(p_disp_th_t)
    error('p_break_t et p_disp_t doivent avoir la même taille.');
end
if numel(t_pdisp_th) < 2
    error('La grille temporelle théorique doit contenir au moins deux points.');
end

dt_pdisp = median(diff(t_pdisp_th));

if isfield(Sp, 'beta0_geom_t')
    beta0_th_t = double(Sp.beta0_geom_t(:));
    if numel(beta0_th_t) ~= numel(p_disp_th_t)
        beta0_th_t = interp1( ...
            linspace(t_pdisp_th(1), t_pdisp_th(end), numel(beta0_th_t)), ...
            beta0_th_t, t_pdisp_th, 'linear', 'extrap');
    end
    beta0_source = 'beta0 théorique temporel';
elseif isfield(Sp, 'beta0_geom_merge')
    beta0_geom = double(Sp.beta0_geom_merge);
    beta0_th_t = beta0_geom * ones(size(p_disp_th_t));
    beta0_source = 'beta0 théorique constant du modèle de fusion';
else
    warning(['Aucun beta0 théorique trouvé. ', ...
        'Utilisation de beta0 = moyenne temporelle empirique du zigzag.']);
    beta0_th_t = mean(h0_dims) * ones(size(p_disp_th_t));
    beta0_source = 'beta0 moyen empirique de repli';
end

B_birth_th = max(beta0_th_t - 1, 0) .* max(p_break_th_t, 0);
if sum(B_birth_th, 'omitnan') <= 0
    warning(['Les poids analytiques de naissance sont nuls. ', ...
        'Utilisation de poids uniformes.']);
    birth_weights = ones(size(B_birth_th)) / numel(B_birth_th);
else
    birth_weights = B_birth_th / sum(B_birth_th, 'omitnan');
end

p_disp_th_t = max(0, min(1-1e-12, 2*p_disp_th_t));
hazard_th_t = -log(1-p_disp_th_t) / dt_pdisp;
period_pdisp = numel(hazard_th_t) * dt_pdisp;

survival_th_variable = zeros(size(Lvals));
for ii = 1:numel(Lvals)
    ell = Lvals(ii);
    S_cond = zeros(size(t_pdisp_th));
    for jj = 1:numel(t_pdisp_th)
        integrated_hazard = periodic_piecewise_hazard_integral( ...
            hazard_th_t, dt_pdisp, jj, ell);
        S_cond(jj) = exp(-integrated_hazard);
    end
    survival_th_variable(ii) = sum(birth_weights .* S_cond, 'omitnan');
end

hazard_mean = sum(birth_weights .* hazard_th_t, 'omitnan');
survival_th_constant = exp(-hazard_mean * Lvals);
if hazard_mean > 0
    tau_mean = 1 / hazard_mean;
else
    tau_mean = Inf;
end

figure;
semilogy(Lvals, survival_emp, 'o-', 'LineWidth', 1.5, ...
    'DisplayName', 'Données simulées');
hold on; grid on;
semilogy(Lvals, survival_th_variable, '--', 'LineWidth', 2.2, ...
    'DisplayName', 'Modèle p_{disp}(t) variable');
semilogy(Lvals, survival_th_constant, ':', 'LineWidth', 1.8, ...
    'DisplayName', 'Exponentielle au risque moyen');
xlabel('Durée des barres (s)');
ylabel('Probabilité de survie');
title(sprintf('Survie des barres H0 - Delta spatial, i = %.1f deg', inc_deg));
legend('Location', 'best');
hold off;

figure;
yyaxis left;
plot(t_pdisp_th, B_birth_th, 'LineWidth', 1.5);
ylabel('B_j^{th}');
yyaxis right;
plot(t_pdisp_th, birth_weights, '--', 'LineWidth', 1.4);
ylabel('Poids de naissance w_j');
grid on;
xlabel('Temps / phase de naissance (s)');
title('Poids analytiques des phases de naissance');

figure;
yyaxis left;
plot(t_pdisp_th, p_disp_th_t, 'LineWidth', 1.5);
ylabel('p_{disp}^{th}(t)');
yyaxis right;
plot(t_pdisp_th, hazard_th_t, '--', 'LineWidth', 1.4);
ylabel('h(t) (s^{-1})');
grid on;
xlabel('Temps (s)');
title('Probabilité de disparition et taux instantané équivalent');

fprintf(' --- Modèle de survie à p_disp(t) variable --- ');
fprintf('Fichier p_disp chargé : %s', pdisp_file);
fprintf('Source de beta0 : %s', beta0_source);
fprintf('Pas temporel du modèle p_disp : %.2f s', dt_pdisp);
fprintf('Période utilisée pour le prolongement : %.2f s', period_pdisp);
fprintf('Somme des naissances théoriques B_j : %.6f', sum(B_birth_th, 'omitnan'));
fprintf('Taux moyen pondéré : %.6e s^-1', hazard_mean);
fprintf('Temps caractéristique au risque moyen : %.2f s', tau_mean);

save('survie_H0_delta_pdisp_variable_results.mat', ...
    'Lvals', 'survival_emp', 'survival_th_variable', 'survival_th_constant', ...
    't_pdisp_th', 'p_disp_th_t', 'p_break_th_t', 'beta0_th_t', ...
    'B_birth_th', 'birth_weights', 'hazard_th_t', 'hazard_mean', ...
    'tau_mean', 'dt_pdisp', 'period_pdisp', 'pdisp_file', 'beta0_source');

fprintf('Résultats sauvegardés dans survie_H0_delta_pdisp_variable_results.mat');

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

function integral_h = periodic_piecewise_hazard_integral(hazard_values, dt_h, start_index, duration)
    if duration <= 0
        integral_h = 0;
        return;
    end
    M = numel(hazard_values);
    remaining = duration;
    idx_h = start_index;
    integral_h = 0;
    period_duration = M * dt_h;
    if remaining >= period_duration
        n_periods = floor(remaining / period_duration);
        integral_h = integral_h + n_periods * dt_h * sum(hazard_values);
        remaining = remaining - n_periods * period_duration;
    end
    while remaining > 1e-12
        segment_duration = min(dt_h, remaining);
        integral_h = integral_h + hazard_values(idx_h) * segment_duration;
        remaining = remaining - segment_duration;
        idx_h = idx_h + 1;
        if idx_h > M
            idx_h = 1;
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