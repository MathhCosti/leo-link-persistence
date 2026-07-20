clear; clc; close all;

%% ============================================================
%  PROBABILITES EMPIRIQUES p_merge, p_break ET p_change
%  A PARTIR DES RESULTATS DU BARCODE ZIGZAG H0
%
%  Zigzag utilise :
%      G_k -> G_k U G_{k+1} <- G_{k+1}
%
%  Convention des indices dans le barcode :
%  - indice impair  2k-1 : graphe reel G_k ;
%  - indice pair    2k   : union G_k U G_{k+1}.
%
%  Evenements :
%  - une fusion entre t_k et t_{k+1} correspond a une barre qui meurt
%    a l'indice 2k-1 lors du passage G_k -> G_k U G_{k+1} ;
%  - une rupture correspond a une barre qui nait a l'indice 2k+1 lors
%    du passage G_k U G_{k+1} <- G_{k+1}.
%
%  Normalisations choisies :
%      p_merge(k) = nb_fusions(k) / beta0(G_k)
%      p_break(k) = nb_ruptures(k) / beta0(G_{k+1})
%
%  Puis :
%      p_change(k) = 1 - (1-p_merge(k))*(1-p_break(k))
%
%  Cette derniere expression represente la probabilite qu'au moins un
%  des deux mecanismes de changement topologique se produise.
%% ============================================================

%% Fichiers d'entree
barcode_file = 'leo_H0_zigzag_random_vectors_barcodes.mat';
analysis_file = 'leo_zigzag_analysis_random_vectors_results.mat';

% Compatibilite avec l'ancien nom du fichier barcode.
if ~isfile(barcode_file)
    barcode_file = 'leo_H0_zigzag_barcodes_delta.mat';
end

if ~isfile(barcode_file)
    error(['Fichier de barcode introuvable. Lancez d''abord le script ', ...
           'barcodes_delta_uniformite_spatiale.m.']);
end

%% Chargement du barcode
Sbar = load(barcode_file, ...
    'birth_index', 'death_index', 'intervals', ...
    'ZigzagTime', 'ZigzagLabels', 'h0_dims');

birth_index = Sbar.birth_index(:);
death_index = Sbar.death_index(:);
ZigzagTime = Sbar.ZigzagTime(:);
ZigzagLabels = Sbar.ZigzagLabels(:);
h0_dims = Sbar.h0_dims(:);

Nz = length(h0_dims);

if mod(Nz,2) ~= 1
    error('Le nombre d''objets du zigzag doit etre impair : Nz = 2*Nt-1.');
end

Nt = (Nz+1)/2;
Ntrans = Nt-1;

%% Recuperation facultative des temps reels et parametres
if isfile(analysis_file)
    Sana = load(analysis_file, 'time_values', 'dt', 'inc_deg', 'N');
else
    Sana = struct();
end

if isfield(Sana, 'time_values') && length(Sana.time_values) == Nt
    time_values = Sana.time_values(:);
else
    % Les graphes reels sont situes aux indices impairs du zigzag.
    time_values = ZigzagTime(1:2:end);
end

if isfield(Sana, 'dt')
    dt = Sana.dt;
else
    dt = median(diff(time_values));
end

if isfield(Sana, 'inc_deg')
    inc_deg = Sana.inc_deg;
else
    inc_deg = NaN;
end

%% ============================================================
%  1. COMPTAGE DES FUSIONS ET DES RUPTURES
%% ============================================================

merge_count = zeros(Ntrans,1);
break_count = zeros(Ntrans,1);

beta0_before = zeros(Ntrans,1);
beta0_after = zeros(Ntrans,1);
beta0_union = zeros(Ntrans,1);

for k = 1:Ntrans
    idx_Gk = 2*k - 1;
    idx_union = 2*k;
    idx_Gkp1 = 2*k + 1;

    % Barres qui disparaissent lors de G_k -> union.
    merge_count(k) = sum(death_index == idx_Gk);

    % Barres qui apparaissent lors de union <- G_{k+1}.
    break_count(k) = sum(birth_index == idx_Gkp1);

    beta0_before(k) = h0_dims(idx_Gk);
    beta0_union(k) = h0_dims(idx_union);
    beta0_after(k) = h0_dims(idx_Gkp1);
end

%% Verifications topologiques directes
% Pour H0, les comptages issus des barres doivent coincider avec :
%   fusion  = beta0(G_k)     - beta0(union)
%   rupture = beta0(G_{k+1}) - beta0(union)
merge_count_from_beta0 = beta0_before - beta0_union;
break_count_from_beta0 = beta0_after - beta0_union;

max_merge_error = max(abs(merge_count - merge_count_from_beta0));
max_break_error = max(abs(break_count - break_count_from_beta0));

fprintf('Erreur max comptage fusion barcode/beta0  : %g\n', max_merge_error);
fprintf('Erreur max comptage rupture barcode/beta0 : %g\n', max_break_error);

if max_merge_error > 0 || max_break_error > 0
    warning(['Les comptages barcode et les differences de beta0 ne ', ...
             'coincident pas exactement. Verifiez la convention des ', ...
             'bornes d''intervalles du calcul de barcode.']);
end

%% ============================================================
%  2. PROBABILITES EMPIRIQUES
%% ============================================================

p_merge = merge_count ./ max(beta0_before, 1);
p_break = break_count ./ max(beta0_after, 1);

% Probabilite combinee de dispersion/changement topologique.
p_change = 1 - (1-p_merge).*(1-p_break);

% Version additive disponible pour comparaison.
p_change_additive = min(1, p_merge + p_break);

% Axe temporel place au milieu de chaque transition.
time_transition = 0.5*(time_values(1:end-1) + time_values(2:end));

%% Moyennes glissantes
% Nombre de transitions utilisees dans chaque fenetre.
% Par exemple, avec dt = 20 s et moving_window = 15, la tendance est
% lissee sur environ 300 s.
moving_window = 15;
moving_window = min(moving_window, Ntrans);

p_merge_moving = movmean(p_merge, moving_window, 'Endpoints', 'shrink');
p_break_moving = movmean(p_break, moving_window, 'Endpoints', 'shrink');
p_change_moving = movmean(p_change, moving_window, 'Endpoints', 'shrink');
p_change_additive_moving = movmean(p_change_additive, moving_window, ...
    'Endpoints', 'shrink');

% Moyennes glissantes des nombres bruts d'evenements.
merge_count_moving = movmean(merge_count, moving_window, ...
    'Endpoints', 'shrink');
break_count_moving = movmean(break_count, moving_window, ...
    'Endpoints', 'shrink');

%% Moyennes globales ponderees
% Les moyennes ponderees correspondent aux rapports entre le nombre total
% d'evenements et le nombre total de composantes exposees.
p_merge_mean = sum(merge_count) / max(sum(beta0_before), 1);
p_break_mean = sum(break_count) / max(sum(beta0_after), 1);
p_change_mean = 1 - (1-p_merge_mean)*(1-p_break_mean);

% Moyennes temporelles simples, utiles si chaque transition doit avoir le
% meme poids independamment du nombre de composantes presentes.
p_merge_time_mean = mean(p_merge);
p_break_time_mean = mean(p_break);
p_change_time_mean = mean(p_change);

%% ============================================================
%  3. TRACES TEMPORELS
%% ============================================================

%% Probabilites brutes et moyennes glissantes
figure;
hold on;
grid on;

% Donnees instantanees en traits fins.
plot(time_transition, p_merge, '-', 'LineWidth', 0.7, ...
    'DisplayName', 'p_{merge} instantane');
plot(time_transition, p_break, '-', 'LineWidth', 0.7, ...
    'DisplayName', 'p_{break} instantane');
plot(time_transition, p_change, '-', 'LineWidth', 0.8, ...
    'DisplayName', 'p_{change} instantane');

% Tendances lissees en traits epais.
plot(time_transition, p_merge_moving, 'LineWidth', 2.2, ...
    'DisplayName', sprintf('p_{merge} moyenne glissante (%d points)', ...
    moving_window));
plot(time_transition, p_break_moving, 'LineWidth', 2.2, ...
    'DisplayName', sprintf('p_{break} moyenne glissante (%d points)', ...
    moving_window));
plot(time_transition, p_change_moving, 'LineWidth', 2.6, ...
    'DisplayName', sprintf('p_{change} moyenne glissante (%d points)', ...
    moving_window));

xlabel('Temps au milieu de la transition (s)');
ylabel('Probabilite empirique');

if isnan(inc_deg)
    title('Probabilites instantanees et moyennes glissantes');
else
    title(sprintf(['Probabilites instantanees et moyennes glissantes ', ...
        '- Delta spatial, i = %.1f deg'], inc_deg));
end

legend('Location', 'best');
ylim([0, min(1, 1.05*max([p_disp; p_disp_moving; eps]))]);
hold off;

%% Figure plus lisible : moyennes glissantes seules
figure;
hold on;
grid on;
plot(time_transition, p_merge_moving, 'LineWidth', 2.0, ...
    'DisplayName', 'p_{merge} lisse');
plot(time_transition, p_break_moving, 'LineWidth', 2.0, ...
    'DisplayName', 'p_{break} lisse');
plot(time_transition, p_change_moving, 'LineWidth', 2.5, ...
    'DisplayName', 'p_{change} lisse');

yline(p_merge_mean, '--', 'LineWidth', 1.0, ...
    'DisplayName', sprintf('moyenne globale p_{merge} = %.3f', ...
    p_merge_mean));
yline(p_break_mean, '--', 'LineWidth', 1.0, ...
    'DisplayName', sprintf('moyenne globale p_{break} = %.3f', ...
    p_break_mean));
yline(p_change_mean, '--', 'LineWidth', 1.2, ...
    'DisplayName', sprintf('moyenne globale p_{disp} = %.3f', ...
    p_change_mean));

xlabel('Temps au milieu de la transition (s)');
ylabel('Probabilite empirique lissee');
title(sprintf('Moyennes glissantes sur %d transitions', moving_window));
legend('Location', 'best');
ylim([0, min(1, 1.05*max([p_change_moving; eps]))]);
hold off;

%% Comptages bruts et leurs moyennes glissantes
figure;
hold on;
grid on;
stairs(time_transition, merge_count, 'LineWidth', 0.7, ...
    'DisplayName', 'Fusions instantanees');
stairs(time_transition, break_count, 'LineWidth', 0.7, ...
    'DisplayName', 'Ruptures instantanees');
plot(time_transition, merge_count_moving, 'LineWidth', 2.0, ...
    'DisplayName', 'Fusions - moyenne glissante');
plot(time_transition, break_count_moving, 'LineWidth', 2.0, ...
    'DisplayName', 'Ruptures - moyenne glissante');
xlabel('Temps au milieu de la transition (s)');
ylabel('Nombre d''evenements');
title('Evenements topologiques et moyennes glissantes');
legend('Location', 'best');
hold off;

%% Comparaison des deux definitions possibles de p_change
figure;
hold on;
grid on;
plot(time_transition, p_change, '-', 'LineWidth', 0.7, ...
    'DisplayName', 'p_{change} combine instantane');
plot(time_transition, p_change_additive, '--', 'LineWidth', 0.7, ...
    'DisplayName', 'p_{change} additif instantane');
plot(time_transition, p_change_moving, 'LineWidth', 2.3, ...
    'DisplayName', 'p_{change} combine lisse');
plot(time_transition, p_change_additive_moving, '--', 'LineWidth', 2.3, ...
    'DisplayName', 'p_{change} additif lisse');
xlabel('Temps au milieu de la transition (s)');
ylabel('p_{change}');
title('Comparaison des definitions de p_{change} avec lissage');
legend('Location', 'best');
hold off;

%% ============================================================
%  4. AFFICHAGE DES RESULTATS
%% ============================================================

fprintf('\n--- Probabilites issues du barcode H0 ---\n');
fprintf('Nombre de transitions : %d\n', Ntrans);
fprintf('Pas temporel moyen     : %.3f s\n', dt);
fprintf('Nombre total de fusions  : %d\n', sum(merge_count));
fprintf('Nombre total de ruptures : %d\n', sum(break_count));

fprintf('\nMoyennes ponderees par le nombre de composantes :\n');
fprintf('p_merge = %.6f\n', p_merge_mean);
fprintf('p_break = %.6f\n', p_break_mean);
fprintf('p_change  = %.6f\n', p_change_mean);

fprintf('\nMoyennes temporelles simples :\n');
fprintf('mean[p_merge(t)] = %.6f\n', p_merge_time_mean);
fprintf('mean[p_break(t)] = %.6f\n', p_break_time_mean);
fprintf('mean[p_change(t)]  = %.6f\n', p_change_time_mean);

%% ============================================================
%  5. SAUVEGARDE
%% ============================================================

output_file = 'pmerge_pbreak_pdisp_barcodes_delta_spatial.mat';

save(output_file, ...
    'time_values', 'time_transition', 'dt', ...
    'merge_count', 'break_count', ...
    'merge_count_from_beta0', 'break_count_from_beta0', ...
    'beta0_before', 'beta0_union', 'beta0_after', ...
    'p_merge', 'p_break', 'p_change', 'p_change_additive', ...
    'moving_window', ...
    'p_merge_moving', 'p_break_moving', 'p_change_moving', ...
    'p_change_additive_moving', ...
    'merge_count_moving', 'break_count_moving', ...
    'p_merge_mean', 'p_break_mean', 'p_change_mean', ...
    'p_merge_time_mean', 'p_break_time_mean', 'p_change_time_mean');

fprintf('\nResultats sauvegardes dans %s\n', output_file);
