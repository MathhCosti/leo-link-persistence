clear; clc; close all;

%% ============================================================
%  PROBABILITE EMPIRIQUE DE FUSION p_merge
%  A PARTIR DES RESULTATS DU BARCODE ZIGZAG H0
%
%  Zigzag utilise :
%      G_k -> G_k U G_{k+1} <- G_{k+1}
%
%  Une fusion entre t_k et t_{k+1} correspond a une barre H0
%  qui meurt a l'indice 2k-1 lors du passage :
%      G_k -> G_k U G_{k+1}.
%
%  Normalisation :
%      p_merge(k) = nb_fusions(k) / beta0(G_k)
%% ============================================================

%% 1. Fichiers d'entree
script_dir = fileparts(mfilename('fullpath'));

barcode_file = fullfile(script_dir, ...
    '..', 'barcodes_results.mat');

analysis_file = fullfile(script_dir, ...
    '..', 'analysis_temp_results.mat');

% Compatibilite si le script est place directement dans le sous-dossier.
if ~isfile(barcode_file)
    barcode_file = fullfile(script_dir, 'barcodes_results.mat');
end
if ~isfile(analysis_file)
    analysis_file = fullfile(script_dir, 'analysis_temp_results.mat');
end

% Compatibilite avec l'ancien fichier de barcode.
if ~isfile(barcode_file)
    old_barcode_file = fullfile(script_dir, ...
        'leo_H0_zigzag_barcodes_delta.mat');
    if isfile(old_barcode_file)
        barcode_file = old_barcode_file;
    end
end

if ~isfile(barcode_file)
    error(['Fichier de barcode introuvable :\n%s\n' ...
        'Lancez d''abord barcodes.m.'], barcode_file);
end

%% 2. Chargement des donnees du barcode
Sbar = load(barcode_file);

required_fields = {'death_index', 'ZigzagTime', 'h0_dims'};
for i = 1:numel(required_fields)
    if ~isfield(Sbar, required_fields{i})
        error('Variable ''%s'' absente de %s.', ...
            required_fields{i}, barcode_file);
    end
end

death_index = Sbar.death_index(:);
ZigzagTime = Sbar.ZigzagTime(:);
h0_dims = Sbar.h0_dims(:);

Nz = numel(h0_dims);

if mod(Nz, 2) ~= 1
    error('Le nombre d''objets du zigzag doit verifier Nz = 2*Nt-1.');
end

Nt = (Nz + 1)/2;
Ntrans = Nt - 1;

%% 3. Recuperation des temps et parametres
if isfile(analysis_file)
    Sana = load(analysis_file, 'time_values', 'dt', 'inc_deg', 'N');
else
    Sana = struct();
end

if isfield(Sana, 'time_values') && numel(Sana.time_values) == Nt
    time_values = Sana.time_values(:);
else
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

if isfield(Sana, 'N')
    N = Sana.N;
else
    N = NaN;
end

time_transition = 0.5 * ...
    (time_values(1:end-1) + time_values(2:end));

%% 4. Comptage des fusions
merge_count = zeros(Ntrans, 1);
beta0_before = zeros(Ntrans, 1);
beta0_union = zeros(Ntrans, 1);

for k = 1:Ntrans
    idx_Gk = 2*k - 1;
    idx_union = 2*k;

    % Barres qui meurent lors du passage G_k -> union.
    merge_count(k) = sum(death_index == idx_Gk);

    beta0_before(k) = h0_dims(idx_Gk);
    beta0_union(k) = h0_dims(idx_union);
end

% Verification directe par les nombres de composantes connexes.
merge_count_from_beta0 = beta0_before - beta0_union;
max_merge_error = max(abs(merge_count - merge_count_from_beta0));

fprintf('Erreur max comptage fusion barcode/beta0 : %g\n', ...
    max_merge_error);

if max_merge_error > 0
    warning(['Le comptage des fusions issu du barcode ne coincide pas ', ...
        'exactement avec beta0(G_k)-beta0(union).']);
end

%% 5. Probabilite empirique de fusion
p_merge = merge_count ./ max(beta0_before, 1);

moving_window = min(15, Ntrans);
p_merge_moving = movmean(p_merge, moving_window, ...
    'Endpoints', 'shrink');
merge_count_moving = movmean(merge_count, moving_window, ...
    'Endpoints', 'shrink');

% Moyenne ponderee par le nombre de composantes exposees.
p_merge_mean = sum(merge_count) / max(sum(beta0_before), 1);

% Moyenne temporelle simple.
p_merge_time_mean = mean(p_merge);

%% 6. Traces
figure;
hold on;
grid on;

plot(time_transition, p_merge, '-', ...
    'LineWidth', 0.7, ...
    'DisplayName', 'p_{merge} instantane');

plot(time_transition, p_merge_moving, ...
    'LineWidth', 2.2, ...
    'DisplayName', sprintf( ...
    'p_{merge} moyenne glissante (%d points)', moving_window));

yline(p_merge_mean, '--', ...
    'LineWidth', 1.2, ...
    'DisplayName', sprintf( ...
    'moyenne globale = %.4f', p_merge_mean));

xlabel('Temps au milieu de la transition (s)');
ylabel('Probabilite empirique de fusion');

if isnan(inc_deg)
    title('Probabilite empirique de fusion');
else
    title(sprintf( ...
        'Probabilite empirique de fusion - i = %.1f deg', inc_deg));
end

legend('Location', 'best');
ylim([0, min(1, 1.05 * max([p_merge; p_merge_moving; eps]))]);
hold off;

figure;
hold on;
grid on;

stairs(time_transition, merge_count, ...
    'LineWidth', 0.8, ...
    'DisplayName', 'Fusions instantanees');

plot(time_transition, merge_count_moving, ...
    'LineWidth', 2.0, ...
    'DisplayName', 'Fusions - moyenne glissante');

xlabel('Temps au milieu de la transition (s)');
ylabel('Nombre de fusions');
title('Nombre de fusions entre deux graphes consecutifs');
legend('Location', 'best');
hold off;

%% 7. Affichage des resultats
fprintf('\n--- Probabilite empirique de fusion ---\n');
fprintf('Nombre de transitions            : %d\n', Ntrans);
fprintf('Pas temporel moyen               : %.3f s\n', dt);
fprintf('Nombre total de fusions          : %d\n', sum(merge_count));
fprintf('p_merge moyen pondere            : %.8f\n', p_merge_mean);
fprintf('moyenne temporelle de p_merge(t) : %.8f\n', ...
    p_merge_time_mean);

%% 8. Sauvegarde
output_file = fullfile(script_dir, 'pmerge_emp_results.mat');

save(output_file, ...
    'time_values', 'time_transition', 'dt', 'inc_deg', 'N', ...
    'merge_count', 'merge_count_from_beta0', ...
    'beta0_before', 'beta0_union', ...
    'p_merge', 'moving_window', 'p_merge_moving', ...
    'merge_count_moving', ...
    'p_merge_mean', 'p_merge_time_mean', ...
    'max_merge_error');

fprintf('Resultats sauvegardes dans :\n%s\n', output_file);
