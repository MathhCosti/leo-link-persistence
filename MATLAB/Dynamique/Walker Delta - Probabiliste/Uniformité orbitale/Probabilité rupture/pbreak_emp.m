clear; clc; close all;

%% ============================================================
%  PROBABILITE EMPIRIQUE DE RUPTURE p_break
%  A PARTIR DES RESULTATS DU BARCODE ZIGZAG H0
%
%  Zigzag utilise :
%      G_k -> G_k U G_{k+1} <- G_{k+1}
%
%  Une rupture entre t_k et t_{k+1} correspond a une barre H0
%  qui nait a l'indice 2k+1 lors du passage :
%      G_k U G_{k+1} <- G_{k+1}.
%
%  Normalisation :
%      p_break(k) = nb_ruptures(k) / beta0(G_{k+1})
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

required_fields = {'birth_index', 'ZigzagTime', 'h0_dims'};
for i = 1:numel(required_fields)
    if ~isfield(Sbar, required_fields{i})
        error('Variable ''%s'' absente de %s.', ...
            required_fields{i}, barcode_file);
    end
end

birth_index = Sbar.birth_index(:);
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

%% 4. Comptage des ruptures
break_count = zeros(Ntrans, 1);
beta0_union = zeros(Ntrans, 1);
beta0_after = zeros(Ntrans, 1);

for k = 1:Ntrans
    idx_union = 2*k;
    idx_Gkp1 = 2*k + 1;

    % Barres qui naissent lors du passage union <- G_{k+1}.
    break_count(k) = sum(birth_index == idx_Gkp1);

    beta0_union(k) = h0_dims(idx_union);
    beta0_after(k) = h0_dims(idx_Gkp1);
end

% Verification directe par les nombres de composantes connexes.
break_count_from_beta0 = beta0_after - beta0_union;
max_break_error = max(abs(break_count - break_count_from_beta0));

fprintf('Erreur max comptage rupture barcode/beta0 : %g\n', ...
    max_break_error);

if max_break_error > 0
    warning(['Le comptage des ruptures issu du barcode ne coincide pas ', ...
        'exactement avec beta0(G_{k+1})-beta0(union).']);
end

%% 5. Probabilite empirique de rupture
p_break = break_count ./ max(beta0_after, 1);

moving_window = min(15, Ntrans);
p_break_moving = movmean(p_break, moving_window, ...
    'Endpoints', 'shrink');
break_count_moving = movmean(break_count, moving_window, ...
    'Endpoints', 'shrink');

% Moyenne ponderee par le nombre de composantes exposees.
p_break_mean = sum(break_count) / max(sum(beta0_after), 1);

% Moyenne temporelle simple.
p_break_time_mean = mean(p_break);

%% 6. Traces
figure;
hold on;
grid on;

plot(time_transition, p_break, '-', ...
    'LineWidth', 0.7, ...
    'DisplayName', 'p_{break} instantane');

plot(time_transition, p_break_moving, ...
    'LineWidth', 2.2, ...
    'DisplayName', sprintf( ...
    'p_{break} moyenne glissante (%d points)', moving_window));

yline(p_break_mean, '--', ...
    'LineWidth', 1.2, ...
    'DisplayName', sprintf( ...
    'moyenne globale = %.4f', p_break_mean));

xlabel('Temps au milieu de la transition (s)');
ylabel('Probabilite empirique de rupture');

if isnan(inc_deg)
    title('Probabilite empirique de rupture');
else
    title(sprintf( ...
        'Probabilite empirique de rupture - i = %.1f deg', inc_deg));
end

legend('Location', 'best');
ylim([0, min(1, 1.05 * max([p_break; p_break_moving; eps]))]);
hold off;

figure;
hold on;
grid on;

stairs(time_transition, break_count, ...
    'LineWidth', 0.8, ...
    'DisplayName', 'Ruptures instantanees');

plot(time_transition, break_count_moving, ...
    'LineWidth', 2.0, ...
    'DisplayName', 'Ruptures - moyenne glissante');

xlabel('Temps au milieu de la transition (s)');
ylabel('Nombre de ruptures');
title('Nombre de ruptures entre deux graphes consecutifs');
legend('Location', 'best');
hold off;

%% 7. Affichage des resultats
fprintf('\n--- Probabilite empirique de rupture ---\n');
fprintf('Nombre de transitions            : %d\n', Ntrans);
fprintf('Pas temporel moyen               : %.3f s\n', dt);
fprintf('Nombre total de ruptures         : %d\n', sum(break_count));
fprintf('p_break moyen pondere            : %.8f\n', p_break_mean);
fprintf('moyenne temporelle de p_break(t) : %.8f\n', ...
    p_break_time_mean);

%% 8. Sauvegarde
output_file = fullfile(script_dir, 'pbreak_emp_results.mat');

save(output_file, ...
    'time_values', 'time_transition', 'dt', 'inc_deg', 'N', ...
    'break_count', 'break_count_from_beta0', ...
    'beta0_union', 'beta0_after', ...
    'p_break', 'moving_window', 'p_break_moving', ...
    'break_count_moving', ...
    'p_break_mean', 'p_break_time_mean', ...
    'max_break_error');

fprintf('Resultats sauvegardes dans :\n%s\n', output_file);
