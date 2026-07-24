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
%  MODELE DE SURVIE A p_disp CONSTANT PAR MORCEAUX
%
%  Le signal p_disp^th(t) est d'abord ramené à une demi-période orbitale,
%  puis découpé en M segments. Dans chaque segment m :
%
%      p_disp^th(t) ~= p_m.
%
%  Pour une barre née à la phase t_j et survivant pendant une durée ell,
%  la survie vaut :
%
%      S(ell | t_j)
%      = produit_m (1-p_m)^{n_m(t_j,ell)},
%
%  où n_m(t_j,ell) est le nombre de pas temporels passés dans le segment m.
%
%  La survie globale est obtenue en moyennant sur les phases de naissance.
%% ============================================================

positive_lifetimes = lifetimes(lifetimes > 0);

Lvals = unique(sort(positive_lifetimes));
survival_emp = zeros(size(Lvals));

for ii = 1:numel(Lvals)
    survival_emp(ii) = mean(positive_lifetimes >= Lvals(ii));
end

%% ============================================================
%  CHARGEMENT DE p_disp THEORIQUE
%% ============================================================

script_dir = fileparts(mfilename('fullpath'));
pdisp_file = fullfile(script_dir, 'Probabilité disparition', 'pdisp_modele_delta_spatial_results.mat');


if ~isfile(pdisp_file)
    error('Fichier introuvable : %s', pdisp_file);
end

Sp = load(pdisp_file);

if ~isfield(Sp,'p_disp_t')
    error('La variable p_disp_t est absente de %s.', pdisp_file);
end

p_disp_th_full = double(Sp.p_disp_t(:));

if isfield(Sp,'t_transition')
    t_pdisp_full = double(Sp.t_transition(:));
elseif isfield(Sp,'time_values')
    tt = double(Sp.time_values(:));
    t_pdisp_full = tt(1:numel(p_disp_th_full));
else
    error('Aucune grille temporelle trouvée dans %s.', pdisp_file);
end

if numel(t_pdisp_full) ~= numel(p_disp_th_full)
    error('Tailles incompatibles entre t_transition et p_disp_t.');
end

dt_pdisp = median(diff(t_pdisp_full));

%% ============================================================
%  EXTRACTION D'UNE DEMI-PERIODE ORBITALE
%% ============================================================

T_orb = 2*pi/omega;
T_period = T_orb/2;

M_period = max(1, round(T_period/dt_pdisp));

if M_period > numel(p_disp_th_full)
    error(['Le signal p_disp théorique est plus court qu''une ', ...
        'demi-période orbitale.']);
end

p_period = p_disp_th_full(1:M_period);
t_period = (0:M_period-1)'*dt_pdisp;

p_period = max(0,min(1-1e-12,p_period));

%% ============================================================
%  APPROXIMATION CONSTANTE PAR MORCEAUX
%% ============================================================

% Nombre de segments sur une demi-période.
% 4 correspond naturellement aux phases :
% minimum, croissance, maximum, décroissance.
n_segments = 24;

segment_edges_idx = round(linspace(1,M_period+1,n_segments+1));
segment_edges_idx(1) = 1;
segment_edges_idx(end) = M_period+1;

p_segment = zeros(n_segments,1);
segment_start_time = zeros(n_segments,1);
segment_end_time = zeros(n_segments,1);
segment_duration = zeros(n_segments,1);

segment_index_of_step = zeros(M_period,1);

for m = 1:n_segments
    idx1 = segment_edges_idx(m);
    idx2 = segment_edges_idx(m+1)-1;

    idx2 = max(idx1,min(idx2,M_period));

    idx_seg = idx1:idx2;

    p_segment(m) = mean(p_period(idx_seg),'omitnan');

    segment_start_time(m) = (idx1-1)*dt_pdisp;
    segment_end_time(m) = idx2*dt_pdisp;
    segment_duration(m) = numel(idx_seg)*dt_pdisp;

    segment_index_of_step(idx_seg) = m;
end

% Signal constant par morceaux reconstruit sur la demi-période.
p_piecewise_period = p_segment(segment_index_of_step);

%% ============================================================
%  PHASES DE NAISSANCE
%% ============================================================

% Version analytique simple : les phases de naissance sont pondérées
% uniformément sur une demi-période.
%
% Pour tester les vraies phases de naissance, remplacer use_empirical_births
% par true.

use_empirical_births = false;

if use_empirical_births
    birth_phase = mod(birth_time,T_period);
    birth_step = floor(birth_phase/dt_pdisp)+1;
    birth_step = max(1,min(M_period,birth_step));

    birth_weights_step = accumarray( ...
        birth_step,1,[M_period,1],@sum,0);

    if sum(birth_weights_step)>0
        birth_weights_step = ...
            birth_weights_step/sum(birth_weights_step);
    else
        birth_weights_step = ones(M_period,1)/M_period;
    end

    birth_weight_type = 'phases de naissance empiriques';
else
    birth_weights_step = ones(M_period,1)/M_period;
    birth_weight_type = 'phases de naissance uniformes';
end

%% ============================================================
%  SURVIE THEORIQUE CONSTANTE PAR MORCEAUX
%% ============================================================

survival_piecewise = zeros(size(Lvals));

for ii = 1:numel(Lvals)
    ell = Lvals(ii);

    S_by_birth_phase = zeros(M_period,1);

    for jj = 1:M_period
        S_by_birth_phase(jj) = ...
            survival_piecewise_periodic( ...
                p_piecewise_period,dt_pdisp,jj,ell);
    end

    survival_piecewise(ii) = ...
        sum(birth_weights_step.*S_by_birth_phase,'omitnan');
end

%% Référence exponentielle utilisant la probabilité moyenne
p_mean = mean(p_period,'omitnan');

tau_mean = -dt_pdisp/log(max(1-2*p_mean,eps));
survival_constant = exp(-Lvals/tau_mean);

%% ============================================================
%  FIGURES
%% ============================================================

figure;
semilogy(Lvals,survival_emp,'o-', ...
    'LineWidth',1.5, ...
    'DisplayName','Données simulées');
hold on;
grid on;

semilogy(Lvals,survival_piecewise,'--', ...
    'LineWidth',2.2, ...
    'DisplayName',sprintf('p_{disp} constant par %d segments', ...
    n_segments));

semilogy(Lvals,survival_constant,':', ...
    'LineWidth',1.8, ...
    'DisplayName','Exponentielle au p moyen');

xlabel('Durée des barres (s)');
ylabel('Probabilité de survie');
title(sprintf(['Survie des barres H0 - Delta spatial, ', ...
    'i = %.1f deg'],inc_deg));
legend('Location','best');
hold off;

%% Comparaison du signal original et de l'approximation par morceaux
figure;
plot(t_period,p_period,'LineWidth',1.6, ...
    'DisplayName','p_{disp}^{th}(t)');
hold on;
stairs(t_period,p_piecewise_period,'--','LineWidth',2, ...
    'DisplayName','Approximation constante par morceaux');

for m = 2:n_segments
    xline(segment_start_time(m),':');
end

grid on;
xlabel('Temps dans une demi-période (s)');
ylabel('p_{disp}(t)');
title('Approximation de p_{disp}(t) par probabilités constantes');
legend('Location','best');
hold off;

%% Valeurs des probabilités par segment
segment_table = table( ...
    (1:n_segments)', ...
    segment_start_time, ...
    segment_end_time, ...
    segment_duration, ...
    p_segment, ...
    'VariableNames',{ ...
    'segment', ...
    't_start', ...
    't_end', ...
    'duration', ...
    'p_disp_segment'});

disp(segment_table);

fprintf('\n--- Modèle de survie constant par morceaux ---\n');
fprintf('Fichier p_disp chargé            : %s\n',pdisp_file);
fprintf('Demi-période orbitale            : %.2f s\n',T_period);
fprintf('Pas temporel p_disp              : %.2f s\n',dt_pdisp);
fprintf('Nombre de segments               : %d\n',n_segments);
fprintf('Pondération des naissances       : %s\n',birth_weight_type);
fprintf('p_disp moyen                     : %.6f\n',p_mean);
fprintf('Temps caractéristique moyen      : %.2f s\n',tau_mean);

%% Sauvegarde
save('survie_H0_delta_pdisp_constant_morceaux_results.mat', ...
    'Lvals','survival_emp','survival_piecewise', ...
    'survival_constant', ...
    'p_period','p_piecewise_period','p_segment', ...
    'segment_edges_idx','segment_start_time', ...
    'segment_end_time','segment_duration','segment_table', ...
    'n_segments','dt_pdisp','T_orb','T_period', ...
    'birth_weights_step','birth_weight_type', ...
    'p_mean','tau_mean','pdisp_file');

fprintf('\nRésultats sauvegardés dans ');
fprintf('survie_H0_delta_pdisp_constant_morceaux_results.mat\n');

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

function S = survival_piecewise_periodic( ...
    p_steps_period,dt_step,start_index,duration)
% Calcule la survie avec une probabilité constante sur chaque pas,
% prolongée périodiquement.
%
% Le signal p_steps_period contient déjà l'approximation constante
% par morceaux. La fonction multiplie les facteurs (1-p) rencontrés
% pendant la durée demandée.

    if duration <= 0
        S = 1;
        return;
    end

    M = numel(p_steps_period);

    % Nombre de pas entiers et fraction du dernier pas.
    n_full_steps = floor(duration/dt_step);
    remainder = duration-n_full_steps*dt_step;

    logS = 0;
    idx = start_index;

    % Accélération pour les périodes complètes.
    if n_full_steps >= M
        n_periods = floor(n_full_steps/M);

        logS = logS + ...
            n_periods*sum(log(max(1-p_steps_period,eps)));

        n_full_steps = n_full_steps-n_periods*M;
    end

    for kk = 1:n_full_steps
        logS = logS + log(max(1-p_steps_period(idx),eps));

        idx = idx+1;
        if idx>M
            idx = 1;
        end
    end

    % Fraction du dernier pas : hypothèse d'un taux constant sur le pas.
    if remainder > 1e-12
        fractional_power = remainder/dt_step;

        logS = logS + fractional_power* ...
            log(max(1-p_steps_period(idx),eps));
    end

    S = exp(logS);
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