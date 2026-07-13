clear; clc; close all;

%% ============================================================
%  p_merge(t), p_break(t) ET p_disp(t) EMPIRIQUES
%  Orbites aleatoires a inclinaison fixe
%
%  Pour chaque transition G_k -> G_{k+1}, on suit les composantes
%  connexes par leurs ensembles de satellites.
%
%  Une composante C de G_k subit :
%    - une rupture si ses sommets se retrouvent dans au moins deux
%      composantes distinctes de G_{k+1};
%    - une fusion si la composante de G_{k+1} contenant ses sommets
%      contient aussi des sommets provenant d'une autre composante de G_k.
%
%  Une composante peut subir simultanement rupture et fusion.
%
%  Definitions :
%    p_break(t_k) = N_break(t_k) / N_comp(t_k)
%    p_merge(t_k) = N_merge(t_k) / N_comp(t_k)
%    p_disp,union(t_k) = N_(break OU merge)(t_k) / N_comp(t_k)
%
%  On compare aussi :
%    p_disp,indep = p_break + p_merge - p_break*p_merge
%
%  Ce dernier terme suppose l'independance des evenements rupture/fusion.
%% ============================================================

analysis_file = 'leo_zigzag_analysis_results_delta.mat';

if ~isfile(analysis_file)
    error(['Fichier introuvable : %s.\n' ...
           'Lance d''abord analysis_temp_random_init.m.'], analysis_file);
end

S = load(analysis_file, 'Adjacency', 'time_values', 'dt', 'inc_deg');

Adjacency   = S.Adjacency;
time_values = S.time_values(:);
dt          = S.dt;

if isfield(S, 'inc_deg')
    inc_deg = S.inc_deg;
else
    inc_deg = NaN;
end

Nt = numel(Adjacency);

if Nt < 2
    error('Il faut au moins deux graphes temporels.');
end

if numel(time_values) ~= Nt
    error('Le nombre de temps ne correspond pas au nombre de graphes.');
end

%% ============================================================
%  Calcul temporel
%% ============================================================

t_events = time_values(1:end-1);

n_components = zeros(Nt-1,1);
n_break      = zeros(Nt-1,1);
n_merge      = zeros(Nt-1,1);
n_both       = zeros(Nt-1,1);
n_disp_union = zeros(Nt-1,1);

p_break_emp     = NaN(Nt-1,1);
p_merge_emp     = NaN(Nt-1,1);
p_both_emp      = NaN(Nt-1,1);
p_disp_union    = NaN(Nt-1,1);
p_disp_indep    = NaN(Nt-1,1);

for k = 1:Nt-1
    A0 = Adjacency{k};
    A1 = Adjacency{k+1};

    stats = component_transition_stats(A0, A1);

    n_components(k) = stats.n_source;
    n_break(k)      = stats.n_break;
    n_merge(k)      = stats.n_merge;
    n_both(k)       = stats.n_both;
    n_disp_union(k) = stats.n_union;

    if stats.n_source > 0
        p_break_emp(k)  = stats.n_break / stats.n_source;
        p_merge_emp(k)  = stats.n_merge / stats.n_source;
        p_both_emp(k)   = stats.n_both  / stats.n_source;
        p_disp_union(k) = stats.n_union / stats.n_source;

        p_disp_indep(k) = p_break_emp(k) + p_merge_emp(k) ...
                        - p_break_emp(k)*p_merge_emp(k);
    end
end

%% Moyennes globales ponderees par le nombre de composantes exposees
den = sum(n_components);

if den > 0
    p_break_global  = sum(n_break)      / den;
    p_merge_global  = sum(n_merge)      / den;
    p_both_global   = sum(n_both)       / den;
    p_disp_global   = sum(n_disp_union) / den;
else
    p_break_global = NaN;
    p_merge_global = NaN;
    p_both_global  = NaN;
    p_disp_global  = NaN;
end

p_disp_indep_global = p_break_global + p_merge_global ...
                    - p_break_global*p_merge_global;

%% Lissage
smoothing_window = 5;

p_break_smooth  = movmean(p_break_emp,  smoothing_window, 'omitnan');
p_merge_smooth  = movmean(p_merge_emp,  smoothing_window, 'omitnan');
p_disp_smooth   = movmean(p_disp_union, smoothing_window, 'omitnan');
p_indep_smooth  = movmean(p_disp_indep, smoothing_window, 'omitnan');

%% ============================================================
%  Figures
%% ============================================================

figure;
hold on;
grid on;

plot(t_events, p_break_emp, '-', 'LineWidth', 0.8, ...
    'DisplayName', 'p_{break} brut');
plot(t_events, p_merge_emp, '-', 'LineWidth', 0.8, ...
    'DisplayName', 'p_{merge} brut');

plot(t_events, p_break_smooth, 'LineWidth', 2.0, ...
    'DisplayName', sprintf('p_{break} moyenne mobile (%d pas)', smoothing_window));
plot(t_events, p_merge_smooth, 'LineWidth', 2.0, ...
    'DisplayName', sprintf('p_{merge} moyenne mobile (%d pas)', smoothing_window));

yline(p_break_global, '--', ...
    sprintf('moyenne p_{break} = %.4f', p_break_global), ...
    'HandleVisibility', 'off');
yline(p_merge_global, '--', ...
    sprintf('moyenne p_{merge} = %.4f', p_merge_global), ...
    'HandleVisibility', 'off');

xlabel('Temps t_k (s)');
ylabel('Probabilite empirique');
title(sprintf('Ruptures et fusions des composantes - inclinaison %.1f deg', inc_deg));
legend('Location', 'best');
hold off;

figure;
hold on;
grid on;

plot(t_events, p_disp_union, '-', 'LineWidth', 0.9, ...
    'DisplayName', 'p_{disp} empirique : break \cup merge');
plot(t_events, p_disp_indep, '--', 'LineWidth', 0.9, ...
    'DisplayName', 'p_{break}+p_{merge}-p_{break}p_{merge}');

plot(t_events, p_disp_smooth, 'LineWidth', 2.2, ...
    'DisplayName', sprintf('p_{disp} empirique lisse (%d pas)', smoothing_window));
plot(t_events, p_indep_smooth, '--', 'LineWidth', 2.2, ...
    'DisplayName', 'Approximation independante lissee');

yline(p_disp_global, ':', ...
    sprintf('moyenne empirique = %.4f', p_disp_global), ...
    'LineWidth', 1.5, 'HandleVisibility', 'off');

xlabel('Temps t_k (s)');
ylabel('p_{disp}(t_k)');
title('Comparaison du taux de disparition des composantes');
legend('Location', 'best');
hold off;

figure;
hold on;
grid on;

plot(t_events, n_components, 'LineWidth', 1.5, ...
    'DisplayName', 'Composantes exposees');
plot(t_events, n_break, 'LineWidth', 1.2, ...
    'DisplayName', 'Composantes rompues');
plot(t_events, n_merge, 'LineWidth', 1.2, ...
    'DisplayName', 'Composantes fusionnees');
plot(t_events, n_both, 'LineWidth', 1.2, ...
    'DisplayName', 'Rupture et fusion simultanees');

xlabel('Temps t_k (s)');
ylabel('Nombre de composantes');
title('Effectifs servant au calcul des probabilites');
legend('Location', 'best');
hold off;

%% ============================================================
%  Console
%% ============================================================

fprintf('\n=== Transitions empiriques de composantes ===\n');
fprintf('Pas temporel dt                           : %.2f s\n', dt);
fprintf('Nombre total d''expositions               : %d\n', den);
fprintf('p_break global                            : %.6f\n', p_break_global);
fprintf('p_merge global                            : %.6f\n', p_merge_global);
fprintf('p_break ET p_merge global                 : %.6f\n', p_both_global);
fprintf('p_disp empirique (union exacte)           : %.6f\n', p_disp_global);
fprintf('p_disp sous hypothese d''independance      : %.6f\n', p_disp_indep_global);
fprintf('Ecart independence - union                : %.6f\n', ...
    p_disp_indep_global - p_disp_global);

%% Sauvegarde
save('pbreak_pmerge_emp_temp_delta_results.mat', ...
    't_events', ...
    'p_break_emp', 'p_merge_emp', 'p_both_emp', ...
    'p_disp_union', 'p_disp_indep', ...
    'p_break_smooth', 'p_merge_smooth', ...
    'p_disp_smooth', 'p_indep_smooth', ...
    'n_components', 'n_break', 'n_merge', 'n_both', 'n_disp_union', ...
    'p_break_global', 'p_merge_global', 'p_both_global', ...
    'p_disp_global', 'p_disp_indep_global', ...
    'smoothing_window', 'dt', 'inc_deg');

fprintf('Resultats sauvegardes dans pbreak_pmerge_emp_temp_delta_results.mat\n');

%% ============================================================
%  Fonction locale
%% ============================================================

function stats = component_transition_stats(A0, A1)
    % Analyse les transitions des composantes de G0 vers G1.
    %
    % La classification est faite du point de vue des composantes
    % presentes a l'instant initial.

    G0 = graph(A0);
    G1 = graph(A1);

    labels0 = conncomp(G0).';
    labels1 = conncomp(G1).';

    n0 = max(labels0);
    n1 = max(labels1);

    % Matrice de recouvrement :
    % overlap(a,b) = nombre de sommets appartenant simultanement
    % a la composante a de G0 et a la composante b de G1.
    overlap = accumarray([labels0, labels1], 1, [n0, n1]);

    is_break = false(n0,1);
    is_merge = false(n0,1);

    % Une composante source est rompue si elle intersecte plusieurs
    % composantes cibles.
    for a = 1:n0
        target_ids = find(overlap(a,:) > 0);
        is_break(a) = numel(target_ids) >= 2;
    end

    % Une composante source participe a une fusion si au moins une
    % composante cible qu'elle atteint recoit aussi des sommets provenant
    % d'une autre composante source.
    for a = 1:n0
        target_ids = find(overlap(a,:) > 0);

        for b = target_ids
            source_ids = find(overlap(:,b) > 0);

            if numel(source_ids) >= 2
                is_merge(a) = true;
                break;
            end
        end
    end

    stats.n_source = n0;
    stats.n_break  = sum(is_break);
    stats.n_merge  = sum(is_merge);
    stats.n_both   = sum(is_break & is_merge);
    stats.n_union  = sum(is_break | is_merge);
end
