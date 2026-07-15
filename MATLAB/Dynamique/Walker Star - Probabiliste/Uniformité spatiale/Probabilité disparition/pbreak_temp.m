clear; clc; close all;

%% ============================================================
% p_break(t) EMPIRIQUE DEPUIS LE BARCODE H0 - WALKER STAR
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
T = load(theory_file, 'time_values', 'p_break_t');

birth_index = B.birth_index(:);
death_index = B.death_index(:);
time_values = A.time_values(:);

Nt = numel(time_values);
Nz = 2*Nt - 1;

persist = (birth_index == 1) & (death_index == Nz);
birth_index(persist) = [];
death_index(persist) = [];

% La rupture G_k U G_{k+1} <- G_{k+1} est datee a t_{k+1}.
t_emp = time_values(2:end);
break_count = zeros(Nt-1,1);
exposed_count = zeros(Nt-1,1);
p_break_emp = NaN(Nt-1,1);

for k = 1:Nt-1
    idx_Gkp1 = 2*k + 1;

    % Une barre naissant dans G_{k+1} correspond a une nouvelle
    % composante creee par fragmentation.
    break_count(k) = nnz(birth_index == idx_Gkp1);

    % Nombre de classes H0 presentes dans G_{k+1}.
    exposed_count(k) = nnz((birth_index <= idx_Gkp1) & ...
                           (death_index >= idx_Gkp1));

    if exposed_count(k) > 0
        p_break_emp(k) = break_count(k) / exposed_count(k);
    end
end

% Interpolation du modele theorique sur les instants du barcode.
t_th_raw = T.time_values(:);
p_th_raw = T.p_break_t(:);
[t_th_raw, ia] = unique(t_th_raw, 'stable');
p_th_raw = p_th_raw(ia);
p_break_th = interp1(t_th_raw, p_th_raw, t_emp, 'linear', 'extrap');
p_break_th = min(max(p_break_th,0),1);

p_break_emp_global = sum(break_count) / max(sum(exposed_count),1);
p_break_th_mean = mean(p_break_th, 'omitnan');

figure;
hold on; grid on;
plot(t_emp, p_break_emp, 'o-', 'LineWidth', 1.0, 'MarkerSize', 3, ...
    'DisplayName', 'p_{break}^{emp}(t)');
plot(t_emp, p_break_th, '--', 'LineWidth', 2.0, ...
    'DisplayName', 'p_{break}^{th}(t)');
yline(p_break_emp_global, ':', 'LineWidth', 1.6, ...
    'DisplayName', sprintf('moyenne empirique = %.4g', p_break_emp_global));
xlabel('Temps (s)');
ylabel('p_{break}(t)');
title('Walker Star : p_{break}(t) empirique et theorique');
legend('Location','best');
local_ylim([p_break_emp; p_break_th; p_break_emp_global]);
hold off;

fprintf('\n=== p_break(t) Walker Star ===\n');
fprintf('Nombre total de ruptures      : %d\n', sum(break_count));
fprintf('Nombre total d''expositions   : %d\n', sum(exposed_count));
fprintf('p_break empirique global      : %.8f\n', p_break_emp_global);
fprintf('p_break theorique moyen       : %.8f\n', p_break_th_mean);
fprintf('RMSE temporelle               : %.8f\n', ...
    sqrt(mean((p_break_emp-p_break_th).^2,'omitnan')));

save('pbreak_temp_barcodes_walker_star_results.mat', ...
    't_emp','p_break_emp','p_break_th','break_count','exposed_count', ...
    'p_break_emp_global','p_break_th_mean');

function local_ylim(values)
    values = values(isfinite(values));
    if isempty(values), return; end
    ymax = max(values);
    if ymax <= 0, ymax = 1; end
    ylim([0, min(1, 1.15*ymax)]);
end
