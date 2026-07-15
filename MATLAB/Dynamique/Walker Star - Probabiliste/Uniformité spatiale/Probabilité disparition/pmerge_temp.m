clear; clc; close all;

%% ============================================================
% p_merge(t) EMPIRIQUE DEPUIS LE BARCODE H0 - WALKER STAR
% Comparaison au modele theorique temporel
%% ============================================================

script_dir = fileparts(mfilename('fullpath'));

barcode_file = fullfile(script_dir, '..', 'leo_H0_zigzag_barcodes.mat');
analysis_file = fullfile(script_dir, '..', 'leo_zigzag_analysis_results.mat');
theory_file = 'pdisp_modele_lambda_vrel_results.mat';

assert(isfile(barcode_file), 'Fichier introuvable : %s', barcode_file);
assert(isfile(analysis_file), 'Fichier introuvable : %s', analysis_file);
assert(isfile(theory_file), 'Fichier introuvable : %s. Lance d''abord le modele theorique.', theory_file);

B = load(barcode_file, 'birth_index', 'death_index');
A = load(analysis_file, 'time_values');
T = load(theory_file, 'time_values', 'p_merge_t');

birth_index = B.birth_index(:);
death_index = B.death_index(:);
time_values = A.time_values(:);

Nt = numel(time_values);
Nz = 2*Nt - 1;

% Suppression eventuelle de la classe globale persistante.
persist = (birth_index == 1) & (death_index == Nz);
birth_index(persist) = [];
death_index(persist) = [];

% La fusion G_k -> G_k U G_{k+1} est datee a t_k.
t_emp = time_values(1:end-1);
merge_count = zeros(Nt-1,1);
exposed_count = zeros(Nt-1,1);
p_merge_emp = NaN(Nt-1,1);

for k = 1:Nt-1
    idx_Gk = 2*k - 1;

    % Une barre mourant a G_k est absorbee par une autre composante
    % dans G_k U G_{k+1} : evenement de fusion en H0.
    merge_count(k) = nnz(death_index == idx_Gk);

    % Nombre de classes H0 presentes dans G_k.
    exposed_count(k) = nnz((birth_index <= idx_Gk) & ...
                           (death_index >= idx_Gk));

    if exposed_count(k) > 0
        p_merge_emp(k) = merge_count(k) / exposed_count(k);
    end
end

% Interpolation du modele theorique sur les instants du barcode.
t_th_raw = T.time_values(:);
p_th_raw = T.p_merge_t(:);
[t_th_raw, ia] = unique(t_th_raw, 'stable');
p_th_raw = p_th_raw(ia);
p_merge_th = interp1(t_th_raw, p_th_raw, t_emp, 'linear', 'extrap');
p_merge_th = min(max(p_merge_th,0),1);

p_merge_emp_global = sum(merge_count) / max(sum(exposed_count),1);
p_merge_th_mean = mean(p_merge_th, 'omitnan');

figure;
hold on; grid on;
plot(t_emp, p_merge_emp, 'o-', 'LineWidth', 1.0, 'MarkerSize', 3, ...
    'DisplayName', 'p_{merge}^{emp}(t)');
plot(t_emp, p_merge_th, '--', 'LineWidth', 2.0, ...
    'DisplayName', 'p_{merge}^{th}(t)');
yline(p_merge_emp_global, ':', 'LineWidth', 1.6, ...
    'DisplayName', sprintf('moyenne empirique = %.4g', p_merge_emp_global));
xlabel('Temps (s)');
ylabel('p_{merge}(t)');
title('Walker Star : p_{merge}(t) empirique et theorique');
legend('Location','best');
local_ylim([p_merge_emp; p_merge_th; p_merge_emp_global]);
hold off;

fprintf('\n=== p_merge(t) Walker Star ===\n');
fprintf('Nombre total de fusions       : %d\n', sum(merge_count));
fprintf('Nombre total d''expositions   : %d\n', sum(exposed_count));
fprintf('p_merge empirique global      : %.8f\n', p_merge_emp_global);
fprintf('p_merge theorique moyen       : %.8f\n', p_merge_th_mean);
fprintf('RMSE temporelle               : %.8f\n', ...
    sqrt(mean((p_merge_emp-p_merge_th).^2,'omitnan')));

save('pmerge_temp_barcodes_walker_star_results.mat', ...
    't_emp','p_merge_emp','p_merge_th','merge_count','exposed_count', ...
    'p_merge_emp_global','p_merge_th_mean');

function local_ylim(values)
    values = values(isfinite(values));
    if isempty(values), return; end
    ymax = max(values);
    if ymax <= 0, ymax = 1; end
    ylim([0, min(1, 1.15*ymax)]);
end
