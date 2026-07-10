clear; clc; close all;

%% ============================================================
%  p_disp(t) MOYENNE SUR PLUSIEURS ITERATIONS
%
%  Ce script calcule p_disp(t) a partir des barres H0 zigzag, exactement
%  comme pdisp_temp.m, mais repete le calcul sur plusieurs realisations.
%
%  Le script est autonome : si un fichier d'iteration n'existe pas,
%  il genere directement une nouvelle realisation dans ce meme fichier.
%
%  Grandeur mesuree :
%  p_disp_emp(t_k) = (# barres H0 vivantes a t_k qui meurent avant t_{k+1})
%                   / (# barres H0 vivantes a t_k)
%% ============================================================

%% Parametres utilisateur
n_iter = 3;

% Si true, le script genere une nouvelle realisation si le fichier
% iteratif correspondant n'existe pas encore.
auto_generate_if_missing = true;

iter_mat_pattern = 'leo_zigzag_analysis_results_iter_%03d.mat';

% Les fichiers de chaque iteration sont sauvegardes automatiquement.

%% Stockage
all_pdisp = [];
all_alive_counts = [];
all_death_counts = [];
all_pth = NaN(n_iter,1);
all_pmerge = NaN(n_iter,1);
all_pbreak = NaN(n_iter,1);
all_f = [];
all_P1 = [];
all_fdom = NaN(n_iter,1);
all_Tdom = NaN(n_iter,1);
all_N = NaN(n_iter,1);
all_nb_bars = NaN(n_iter,1);

time_ref = [];
f_ref = [];
R_ref = NaN;
dt_ref = NaN;
mu = 398600; % km^3/s^2

fprintf('\n=== Calcul de p_disp moyen sur %d iterations ===\n', n_iter);

for it = 1:n_iter

    fprintf('\n--- Iteration %d / %d ---\n', it, n_iter);

    iter_mat_file = sprintf(iter_mat_pattern, it);

    if isfile(iter_mat_file)
        mat_file = iter_mat_file;
        fprintf('Chargement du fichier existant : %s\n', mat_file);
    else
        if ~auto_generate_if_missing
            error('Fichier %s introuvable et auto_generate_if_missing = false.', iter_mat_file);
        end

        fprintf('Generation autonome d''une nouvelle realisation...\n');
        generate_zigzag_realization_mat(iter_mat_file);
        mat_file = iter_mat_file;
        fprintf('Realisation sauvegardee : %s\n', iter_mat_file);
    end

    result = compute_pdisp_from_mat(mat_file);

    all_N(it) = result.N;
    all_nb_bars(it) = result.nb_bars;
    all_pth(it) = result.p_disp_th;
    all_pmerge(it) = result.p_merge;
    all_pbreak(it) = result.p_break;
    all_fdom(it) = result.f_dom;
    all_Tdom(it) = result.T_dom;

    if isempty(time_ref)
        time_ref = result.t_plot(:);
        R_ref = result.R;
        dt_ref = result.dt;
        all_pdisp = NaN(length(time_ref), n_iter);
        all_alive_counts = zeros(length(time_ref), n_iter);
        all_death_counts = zeros(length(time_ref), n_iter);
    end

    % Si la grille temporelle est differente, interpolation sur la grille de reference.
    if length(result.t_plot) == length(time_ref) && max(abs(result.t_plot(:) - time_ref(:))) < 1e-9
        pdisp_it = result.p_disp_emp_t(:);
        alive_it = result.alive_counts(:);
        death_it = result.death_counts(:);
    else
        pdisp_it = interp1(result.t_plot(:), result.p_disp_emp_t(:), time_ref, 'linear', NaN);
        alive_it = interp1(result.t_plot(:), double(result.alive_counts(:)), time_ref, 'nearest', 0);
        death_it = interp1(result.t_plot(:), double(result.death_counts(:)), time_ref, 'nearest', 0);
    end

    all_pdisp(:,it) = pdisp_it;
    all_alive_counts(:,it) = alive_it;
    all_death_counts(:,it) = death_it;

    if isempty(f_ref)
        f_ref = result.f(:);
        all_f = f_ref;
        all_P1 = NaN(length(f_ref), n_iter);
    end

    if ~isempty(result.f)
        if length(result.f) == length(f_ref) && max(abs(result.f(:) - f_ref(:))) < 1e-12
            all_P1(:,it) = result.P1(:);
        else
            all_P1(:,it) = interp1(result.f(:), result.P1(:), f_ref, 'linear', NaN);
        end
    end
end

%% Moyennes temporelles et globales
p_disp_emp_mean_t = mean(all_pdisp, 2, 'omitnan');
p_disp_emp_std_t = std(all_pdisp, 0, 2, 'omitnan');

alive_total_t = sum(all_alive_counts, 2, 'omitnan');
death_total_t = sum(all_death_counts, 2, 'omitnan');
p_disp_emp_global_t = death_total_t ./ alive_total_t;
p_disp_emp_global_t(alive_total_t == 0) = NaN;

p_disp_emp_mean_iter = mean(all_pdisp, 1, 'omitnan');
p_disp_emp_mean_of_iter = mean(p_disp_emp_mean_iter, 'omitnan');

p_disp_emp_global = sum(death_total_t, 'omitnan') / sum(alive_total_t, 'omitnan');
p_disp_th_mean = mean(all_pth, 'omitnan');
p_merge_mean = mean(all_pmerge, 'omitnan');
p_break_mean = mean(all_pbreak, 'omitnan');

fprintf('\n=== Resultats moyens sur %d iterations ===\n', n_iter);
fprintf('N moyen                              : %.2f\n', mean(all_N, 'omitnan'));
fprintf('Nombre moyen de barres H0             : %.2f\n', mean(all_nb_bars, 'omitnan'));
fprintf('p_merge theorique moyen               : %.6f\n', p_merge_mean);
fprintf('p_break theorique moyen               : %.6f\n', p_break_mean);
fprintf('p_disp theorique moyen                : %.6f\n', p_disp_th_mean);
fprintf('p_disp empirique moyen des iterations : %.6f\n', p_disp_emp_mean_of_iter);
fprintf('p_disp empirique global               : %.6f\n', p_disp_emp_global);

%% Graphe temporel moyen de p_disp(t)
figure;
plot(time_ref, p_disp_emp_mean_t, 'b-o', 'LineWidth', 1.4, 'MarkerSize', 3); hold on;
plot(time_ref, p_disp_emp_mean_t + p_disp_emp_std_t, 'b:', 'LineWidth', 0.9);
plot(time_ref, p_disp_emp_mean_t - p_disp_emp_std_t, 'b:', 'LineWidth', 0.9);
yline(p_disp_th_mean, 'r--', ...
    sprintf('p_{disp}^{th} moyen = %.4f', p_disp_th_mean), ...
    'LineWidth', 1.8, 'LabelHorizontalAlignment', 'left');
yline(p_disp_emp_global, 'k:', ...
    sprintf('moyenne empirique globale = %.4f', p_disp_emp_global), ...
    'LineWidth', 2.0, 'LabelHorizontalAlignment', 'left');
grid on;
xlabel('Temps (s)');
ylabel('p_{disp}(t)');
title(sprintf('p_{disp}(t) moyen sur %d iterations', n_iter));
legend('Moyenne empirique', 'Moyenne + ecart-type', 'Moyenne - ecart-type', ...
       'Theorie constante', 'Moyenne empirique globale', 'Location', 'best');

%% Spectre moyen
if ~isempty(all_f) && ~isempty(all_P1)
    P1_mean = mean(all_P1, 2, 'omitnan');
    P1_std = std(all_P1, 0, 2, 'omitnan');

    % Frequence dominante du spectre moyen, hors composante nulle.
    if length(P1_mean) >= 2
        [amp_dom_mean, idx_rel] = max(P1_mean(2:end));
        idx_dom = idx_rel + 1;
        f_dom_mean = all_f(idx_dom);
        T_dom_mean = 1/f_dom_mean;
    else
        amp_dom_mean = NaN;
        f_dom_mean = NaN;
        T_dom_mean = NaN;
    end

    T_orb = 2*pi*sqrt(R_ref^3/mu);
    f_orb = 1/T_orb;
    f_half_orb = 2*f_orb;

    fprintf('\n--- Spectre moyen ---\n');
    fprintf('Frequence dominante du spectre moyen : %.6e Hz\n', f_dom_mean);
    fprintf('Periode dominante du spectre moyen   : %.2f s\n', T_dom_mean);
    fprintf('Amplitude dominante moyenne FFT      : %.6e\n', amp_dom_mean);
    fprintf('Frequence demi-revolution theorique  : %.6e Hz\n', f_half_orb);
    fprintf('Periode demi-revolution theorique    : %.2f s\n', T_orb/2);

    figure;
    plot(all_f, P1_mean, 'b-o', 'LineWidth', 1.3, 'MarkerSize', 3); hold on;
    plot(all_f, P1_mean + P1_std, 'b:', 'LineWidth', 0.8);
    plot(all_f, max(P1_mean - P1_std, 0), 'b:', 'LineWidth', 0.8);

    f1 = f_half_orb;
    f_max_plot = max(all_f);
    n_harm = floor(f_max_plot / f1);

    for n = 1:n_harm
        fh = n * f1;
        if n == 1
            xline(fh, 'r--', sprintf('H_%d', n), ...
                'LineWidth', 1.8, 'LabelHorizontalAlignment', 'left');
        else
            xline(fh, 'm:', sprintf('H_%d', n), ...
                'LineWidth', 1.2, 'LabelHorizontalAlignment', 'left');
        end
    end

    xline(f_dom_mean, 'k:', ...
        sprintf('dominante : %.3e Hz', f_dom_mean), ...
        'LineWidth', 1.8, 'LabelHorizontalAlignment', 'left');

    grid on;
    xlabel('Frequence (Hz)');
    ylabel('Amplitude moyenne');
    title(sprintf('Spectre moyen de p_{disp}(t) sur %d iterations', n_iter));
    legend('Spectre moyen', 'Moyenne + ecart-type', 'Moyenne - ecart-type', ...
           '1ere harmonique theorique', 'Harmoniques multiples', ...
           'Frequence dominante', 'Location', 'best');
else
    P1_mean = [];
    P1_std = [];
    f_dom_mean = NaN;
    T_dom_mean = NaN;
    T_orb = 2*pi*sqrt(R_ref^3/mu);
    f_orb = 1/T_orb;
    f_half_orb = 2*f_orb;
end

%% Sauvegarde
save('pdisp_moyenne_iterations_results.mat', ...
    'n_iter', 'time_ref', 'all_pdisp', 'p_disp_emp_mean_t', ...
    'p_disp_emp_std_t', 'p_disp_emp_global_t', ...
    'p_disp_emp_mean_of_iter', 'p_disp_emp_global', ...
    'p_disp_th_mean', 'p_merge_mean', 'p_break_mean', ...
    'all_pth', 'all_pmerge', 'all_pbreak', 'all_N', 'all_nb_bars', ...
    'all_f', 'all_P1', 'P1_mean', 'P1_std', ...
    'all_fdom', 'all_Tdom', 'f_dom_mean', 'T_dom_mean', ...
    'T_orb', 'f_orb', 'f_half_orb');

fprintf('\nResultats sauvegardes dans pdisp_moyenne_iterations_results.mat\n');


%% ============================================================
%  GENERATION AUTONOME D'UNE REALISATION LEO + ZIGZAG
%% ============================================================

function generate_zigzag_realization_mat(output_file)

    %% Parametres physiques
    R_earth = 6371;      % km
    h = 550;             % km
    R = R_earth + h;     % rayon orbital

    mu = 398600;              % km^3/s^2
    omega = sqrt(mu / R^3);   % rad/s

    %% Parametres du processus de Poisson
    lambda = 4e-7;       % satellites / km^2
    surface_sphere = 4*pi*R^2;
    N = poissrnd(lambda * surface_sphere);

    fprintf('Nombre de satellites generes : N = %d\n', N);

    %% Positions initiales uniformes sur la sphere
    u = rand(N,1);
    phi = 2*pi*rand(N,1);
    theta = acos(1 - 2*u);

    %% Sens de rotation defini par un plan separateur passant par les poles
    y0 = R * sin(theta) .* sin(phi);
    rotation_sign = ones(N, 1);
    rotation_sign(y0 < 0) = -1;

    %% Parametres des liens et du temps
    dmax = 1500;     % km
    dt = 60;         % s
    Tmax = 12000;    % s

    time_values = 0:dt:Tmax;
    Nt = length(time_values);

    Positions = cell(Nt,1);
    Adjacency = cell(Nt,1);

    %% Construction des graphes temporels
    for k = 1:Nt
        t = time_values(k);

        phi_t = phi;
        theta_t = theta + rotation_sign * omega * t;

        x_t = R * sin(theta_t) .* cos(phi_t);
        y_t = R * sin(theta_t) .* sin(phi_t);
        z_t = R * cos(theta_t);

        positions_t = [x_t y_t z_t];

        D = squareform(pdist(positions_t));
        A = (D <= dmax) & (D > 0);
        A = sparse(A);

        Positions{k} = positions_t;
        Adjacency{k} = A;
    end

    %% Construction du zigzag par unions
    Nz = 2*Nt - 1;
    ZigzagAdjacency = cell(Nz,1);
    ZigzagLabels = zeros(Nz,1);

    idx = 1;
    for k = 1:Nt
        ZigzagAdjacency{idx} = Adjacency{k};
        ZigzagLabels(idx) = k;
        idx = idx + 1;

        if k < Nt
            ZigzagAdjacency{idx} = Adjacency{k} | Adjacency{k+1};
            ZigzagLabels(idx) = k + 0.5;
            idx = idx + 1;
        end
    end

    save(output_file, ...
        'N', 'R', 'h', 'lambda', 'dmax', 'dt', 'Tmax', ...
        'time_values', 'Positions', 'Adjacency', ...
        'ZigzagAdjacency', 'ZigzagLabels');
end

%% ============================================================
%  FONCTION PRINCIPALE : calcul p_disp pour un fichier .mat
%% ============================================================

function result = compute_pdisp_from_mat(mat_file)

    load(mat_file, ...
        'ZigzagAdjacency', 'ZigzagLabels', 'time_values', ...
        'N', 'R', 'lambda', 'dmax', 'dt');

    Nz = length(ZigzagAdjacency);
    Nt = length(time_values);
    mu = 398600;              % km^3/s^2

    %% Theorie
    v_orb = sqrt(mu/R);       % km/s
    v_rel = (4/pi) * v_orb;   % km/s

    p_merge = 1 - exp(-2 * lambda * dmax * v_rel * dt);
    p_break = 2 * v_rel * dt / (pi * dmax);

    p_merge = min(max(p_merge, 0), 1);
    p_break = min(max(p_break, 0), 1);

    p_disp_th = 1 - (1 - p_merge) * (1 - p_break);

    %% Conversion labels -> temps physiques
    ZigzagTime = zeros(Nz,1);

    for k = 1:Nz
        lab = ZigzagLabels(k);

        if abs(lab - round(lab)) < 1e-12
            idx_time = round(lab);
            ZigzagTime(k) = time_values(idx_time);
        else
            idx_time = floor(lab);
            ZigzagTime(k) = 0.5 * (time_values(idx_time) + time_values(idx_time+1));
        end
    end

    %% Espaces H0 du zigzag
    component_labels_zigzag = cell(Nz,1);
    h0_dims = zeros(Nz,1);

    for k = 1:Nz
        A = ZigzagAdjacency{k};
        G = graph(A);

        comp = conncomp(G);
        comp = comp(:);

        component_labels_zigzag{k} = comp;
        h0_dims(k) = max(comp);
    end

    %% Applications H0
    maps = cell(Nz-1,1);

    for k = 1:Nz-1
        if mod(k,2) == 1
            maps{k}.type = 'f';
            maps{k}.mat = build_H0_map( ...
                component_labels_zigzag{k}, ...
                component_labels_zigzag{k+1}, ...
                h0_dims(k), ...
                h0_dims(k+1));
        else
            maps{k}.type = 'g';
            maps{k}.mat = build_H0_map( ...
                component_labels_zigzag{k+1}, ...
                component_labels_zigzag{k}, ...
                h0_dims(k+1), ...
                h0_dims(k));
        end
    end

    %% Barcode H0
    intervals_H0 = zigzag_barcode_from_module_mod2(h0_dims, maps);

    birth_index_H0 = intervals_H0(:,1);
    death_index_H0 = intervals_H0(:,2);

    birth_time_H0 = ZigzagTime(birth_index_H0);
    death_time_H0 = ZigzagTime(death_index_H0);
    lifetimes_H0 = death_time_H0 - birth_time_H0;

    %% p_disp empirique temporel
    p_disp_emp_t = NaN(Nt-1, 1);
    alive_counts = zeros(Nt-1, 1);
    death_counts = zeros(Nt-1, 1);

    for k = 1:Nt-1
        t0 = time_values(k);
        t1 = time_values(k+1);

        alive = (birth_time_H0 <= t0) & (death_time_H0 > t0);
        dying = alive & (death_time_H0 <= t1);

        alive_counts(k) = sum(alive);
        death_counts(k) = sum(dying);

        if alive_counts(k) > 0
            p_disp_emp_t(k) = death_counts(k) / alive_counts(k);
        end
    end

    %% FFT
    t_plot = time_values(1:end-1);
    t_freq = t_plot(:);
    x_freq = p_disp_emp_t(:);
    valid = isfinite(x_freq);

    if nnz(valid) >= 4
        x_interp = x_freq;
        if any(~valid)
            x_interp(~valid) = interp1(t_freq(valid), x_freq(valid), ...
                                       t_freq(~valid), 'linear', 'extrap');
        end

        x_centered = x_interp - mean(x_interp, 'omitnan');

        Fs = 1/dt;
        L = length(x_centered);

        Y = fft(x_centered);
        P2 = abs(Y/L);
        P1 = P2(1:floor(L/2)+1);

        if length(P1) > 2
            P1(2:end-1) = 2*P1(2:end-1);
        end

        f = Fs*(0:floor(L/2))/L;

        if length(P1) >= 2
            [~, idx_rel] = max(P1(2:end));
            idx_dom = idx_rel + 1;
            f_dom = f(idx_dom);
            T_dom = 1/f_dom;
        else
            f_dom = NaN;
            T_dom = NaN;
        end
    else
        f = [];
        P1 = [];
        f_dom = NaN;
        T_dom = NaN;
    end

    %% Sortie
    result.N = N;
    result.R = R;
    result.lambda = lambda;
    result.dmax = dmax;
    result.dt = dt;
    result.t_plot = t_plot(:);
    result.p_disp_emp_t = p_disp_emp_t(:);
    result.alive_counts = alive_counts(:);
    result.death_counts = death_counts(:);
    result.p_disp_th = p_disp_th;
    result.p_merge = p_merge;
    result.p_break = p_break;
    result.nb_bars = size(intervals_H0,1);
    result.birth_time_H0 = birth_time_H0;
    result.death_time_H0 = death_time_H0;
    result.lifetimes_H0 = lifetimes_H0;
    result.f = f(:);
    result.P1 = P1(:);
    result.f_dom = f_dom;
    result.T_dom = T_dom;
end


%% ============================================================
%  FONCTIONS LOCALES POUR LE BARCODE ZIGZAG H0
%% ============================================================

function M = build_H0_map(labels_source, labels_target, dim_source, dim_target)
    % Construit la matrice induite en H0 par une inclusion de graphes.
    % Chaque composante source est envoyee vers la composante cible
    % qui la contient. M est de taille dim_target x dim_source.

    M = zeros(dim_target, dim_source);

    for c = 1:dim_source
        vertices = find(labels_source == c);
        target_comps = unique(labels_target(vertices));

        if length(target_comps) ~= 1
            error(['Inclusion invalide pour H0 : une composante source ', ...
                   'est envoyee dans plusieurs composantes cibles.']);
        end

        target_c = target_comps(1);
        M(target_c, c) = 1;
    end

    M = mod(M,2);
end

function intervals = zigzag_barcode_from_module_mod2(dims, maps)
    % Calcule le barcode zigzag d'un module sur F2 a partir des dimensions
    % et des matrices d'applications.

    n = length(dims);

    R = cell(2,1);
    R{1} = zeros(dims(1),0);
    R{2} = eye(dims(1));

    b = 1;
    r = filtration_quotient_dims(R);
    intervals = [];

    for k = 1:n-1
        current_type = maps{k}.type;

        if current_type == 'f'
            M = maps{k}.mat;
            Rnext = cell(length(R)+1,1);

            for i = 1:length(R)
                Rnext{i} = gf2_col_basis(M * R{i});
            end

            Rnext{end} = eye(dims(k+1));
            bnext = [b, k+1];
            rnext = filtration_quotient_dims(Rnext);

            for i = 1:length(r)
                c = r(i) - rnext(i);
                if c > 0
                    intervals = [intervals; repmat([b(i), k], c, 1)]; %#ok<AGROW>
                end
            end

        elseif current_type == 'g'
            Nmat = maps{k}.mat;
            Rnext = cell(length(R)+1,1);
            Rnext{1} = zeros(dims(k+1),0);

            for i = 1:length(R)
                Rnext{i+1} = gf2_preimage(Nmat, R{i});
            end

            bnext = [k+1, b];
            rnext = filtration_quotient_dims(Rnext);

            for i = 1:length(r)
                c = r(i) - rnext(i+1);
                if c > 0
                    intervals = [intervals; repmat([b(i), k], c, 1)]; %#ok<AGROW>
                end
            end
        else
            error('Type de fleche inconnu.');
        end

        R = Rnext;
        b = bnext;
        r = rnext;
    end

    for i = 1:length(r)
        c = r(i);
        if c > 0
            intervals = [intervals; repmat([b(i), n], c, 1)]; %#ok<AGROW>
        end
    end
end

function dims = filtration_quotient_dims(R)
    m = length(R) - 1;
    dims = zeros(1,m);

    for i = 1:m
        dims(i) = gf2_rank(R{i+1}) - gf2_rank(R{i});
    end
end

function P = gf2_preimage(A, S)
    A = mod(full(A),2);
    S = mod(full(S),2);

    n = size(A,2);
    Big = [A S];
    Z = gf2_null(Big);
    X = Z(1:n,:);
    P = gf2_col_basis(X);
end

function B = gf2_col_basis(A)
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
    A = mod(full(A),2);

    if isempty(A)
        r = 0;
        return;
    end

    [~, pivots] = gf2_rref(A);
    r = length(pivots);
end

function Z = gf2_null(A)
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

        if pivot ~= row
            tmp = R(row,:);
            R(row,:) = R(pivot,:);
            R(pivot,:) = tmp;
        end

        for rr = 1:m
            if rr ~= row && R(rr,col) == 1
                R(rr,:) = mod(R(rr,:) + R(row,:), 2);
            end
        end

        pivots(end+1) = col; %#ok<AGROW>
        row = row + 1;
    end
end
