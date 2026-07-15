%% comparaison_pdisp_th_emp_spectres.m
% Compare p_disp theorique et empirique, puis compare leurs spectres.
% Les artefacts de bord du barcode sont retires avant les calculs.

clear; clc; close all;

%% Chargement des resultats
load('pdisp_modele_lambda_vrel_results.mat', ...
    'time_values', 'p_disp_t', 'p_merge_t', 'p_break_t');

load('pdisp_moyenne_iterations_results.mat', ...
    'time_ref', 'p_disp_emp_mean_t', 'p_disp_emp_std_t', ...
    'p_disp_emp_global_t');

% Mise en colonnes double
t_th = double(time_values(:));
p_th = double(p_disp_t(:));
p_merge_t = double(p_merge_t(:));
p_break_t = double(p_break_t(:));

t_emp = double(time_ref(:));
p_emp = double(p_disp_emp_mean_t(:));
p_emp_std = double(p_disp_emp_std_t(:));
p_emp_global_t = double(p_disp_emp_global_t(:));

%% Interpolation de la theorie sur les temps empiriques
p_th_i = interp1(t_th, p_th, t_emp, 'linear', 'extrap');
p_merge_i = interp1(t_th, p_merge_t, t_emp, 'linear', 'extrap');
p_break_i = interp1(t_th, p_break_t, t_emp, 'linear', 'extrap');

%% Suppression des effets de bord
edge_trim_steps = 2;

idx_valid = true(size(t_emp));
idx_valid(1:edge_trim_steps) = false;
idx_valid(end-edge_trim_steps+1:end) = false;

% Securite : points finis et pas de pic artificiel de fin de barcode.
idx_valid = idx_valid & isfinite(p_emp) & isfinite(p_emp_std) & isfinite(p_th_i);
idx_valid = idx_valid & (p_emp < 0.95);

% Donnees filtrees
tv = t_emp(idx_valid);
p_emp_v = p_emp(idx_valid);
p_emp_v = movmean(p_emp_v, 10, 'omitnan');
p_emp_std_v = p_emp_std(idx_valid);
p_emp_global_v = p_emp_global_t(idx_valid);
p_th_v = p_th_i(idx_valid);
p_merge_v = p_merge_i(idx_valid);
p_break_v = p_break_i(idx_valid);

mean_emp = mean(p_emp_v, 'omitnan');
mean_emp_global = mean(p_emp_global_v, 'omitnan');
mean_th = mean(p_th_v, 'omitnan');

fprintf('\n=== Comparaison p_disp sans artefacts de bord ===\n');
fprintf('Points conserves : %d / %d\n', nnz(idx_valid), numel(idx_valid));
fprintf('Moyenne empirique temporelle : %.6f\n', mean_emp);
fprintf('Moyenne empirique globale    : %.6f\n', mean_emp_global);
fprintf('Moyenne theorique            : %.6f\n', mean_th);

%% Figure 1 : comparaison principale
figure;
hold on; grid on; box on;

fill([tv; flipud(tv)], ...
     [p_emp_v + p_emp_std_v; flipud(p_emp_v - p_emp_std_v)], ...
     [0.75 0.75 1], 'FaceAlpha', 0.25, 'EdgeColor', 'none');

plot(tv, p_emp_v, 'b-', 'LineWidth', 1.5);
plot(tv, p_emp_global_v, 'c--', 'LineWidth', 1.0);
plot(tv, p_th_v, 'r-', 'LineWidth', 1.7);

yline(mean_emp, 'b:', sprintf('moyenne emp = %.3f', mean_emp), 'LineWidth', 1.2);
yline(mean_th, 'r:', sprintf('moyenne th = %.3f', mean_th), 'LineWidth', 1.2);

xlabel('Temps (s)');
ylabel('p_{disp}(t)');
title('Comparaison de p_{disp}(t) theorique et empirique sans effets de bord');
legend('Empirique \pm ecart-type', ...
       'p_{disp}^{emp}(t) moyen', ...
       'p_{disp}^{emp}(t) global', ...
       'p_{disp}^{th}(t)', ...
       'Location', 'best');

%% Figure 2 : decomposition theorique
figure;
hold on; grid on; box on;
plot(tv, p_merge_v, 'LineWidth', 1.4);
plot(tv, p_break_v, 'LineWidth', 1.4);
plot(tv, p_th_v, 'LineWidth', 1.8);
plot(tv, p_emp_v, 'k-', 'LineWidth', 1.2);
xlabel('Temps (s)');
ylabel('Probabilite par pas');
title('Decomposition theorique et comparaison empirique');
legend('p_{merge}^{th}(t)', 'p_{break}^{th}(t)', ...
       'p_{disp}^{th}(t)', 'p_{disp}^{emp}(t)', ...
       'Location', 'best');

%% ============================================================
%  SPECTRES THEORIQUE ET EMPIRIQUE
%% ============================================================
% On retire la moyenne pour comparer les oscillations, pas les offsets.
% On utilise seulement les points valides deja filtres.

% Verification du pas temporel
if numel(tv) < 4
    error('Pas assez de points valides pour calculer un spectre.');
end

dt_vec = diff(tv);
dt_spec = median(dt_vec);
if max(abs(dt_vec - dt_spec)) > 1e-6 * max(1, dt_spec)
    warning('Le pas temporel n''est pas parfaitement uniforme. Interpolation sur une grille uniforme.');
    t_uniform = (tv(1):dt_spec:tv(end)).';
    p_emp_spec_sig = interp1(tv, p_emp_v, t_uniform, 'linear');
    p_th_spec_sig  = interp1(tv, p_th_v,  t_uniform, 'linear');
else
    t_uniform = tv;
    p_emp_spec_sig = p_emp_v;
    p_th_spec_sig  = p_th_v;
end

% Suppression NaN eventuels apres interpolation
idx_spec = isfinite(p_emp_spec_sig) & isfinite(p_th_spec_sig);
t_uniform = t_uniform(idx_spec);
p_emp_spec_sig = p_emp_spec_sig(idx_spec);
p_th_spec_sig  = p_th_spec_sig(idx_spec);

Nspec = numel(t_uniform);
Fs = 1 / dt_spec;

% Signaux centres
x_emp = p_emp_spec_sig - mean(p_emp_spec_sig, 'omitnan');
x_th  = p_th_spec_sig  - mean(p_th_spec_sig,  'omitnan');

% Fenetre de Hann simple pour limiter la fuite spectrale
if Nspec >= 4
    n = (0:Nspec-1).';
    win = 0.5 * (1 - cos(2*pi*n/(Nspec-1)));
else
    win = ones(Nspec,1);
end

x_emp_w = x_emp .* win;
x_th_w  = x_th  .* win;

Y_emp = fft(x_emp_w);
Y_th  = fft(x_th_w);

P2_emp = abs(Y_emp / Nspec);
P2_th  = abs(Y_th  / Nspec);

n_half = floor(Nspec/2) + 1;
freq = Fs * (0:n_half-1).' / Nspec;

Amp_emp = 2 * P2_emp(1:n_half);
Amp_th  = 2 * P2_th(1:n_half);
Amp_emp(1) = P2_emp(1); % DC, normalement proche de 0 apres centrage
Amp_th(1)  = P2_th(1);

% Frequences dominantes hors DC
if numel(freq) >= 2
    [~, idx_dom_emp_rel] = max(Amp_emp(2:end));
    [~, idx_dom_th_rel]  = max(Amp_th(2:end));
    idx_dom_emp = idx_dom_emp_rel + 1;
    idx_dom_th  = idx_dom_th_rel  + 1;
    f_dom_emp = freq(idx_dom_emp);
    f_dom_th  = freq(idx_dom_th);
else
    f_dom_emp = NaN;
    f_dom_th = NaN;
end

fprintf('\n=== Spectres de p_disp ===\n');
fprintf('Frequence dominante empirique : %.6g Hz  (periode %.2f s)\n', ...
    f_dom_emp, 1/f_dom_emp);
fprintf('Frequence dominante theorique : %.6g Hz  (periode %.2f s)\n', ...
    f_dom_th, 1/f_dom_th);

%% Figure 3 : comparaison des spectres
figure;
hold on; grid on; box on;
plot(freq, Amp_emp, 'b-', 'LineWidth', 1.5);
plot(freq, Amp_th,  'r-', 'LineWidth', 1.5);

if isfinite(f_dom_emp)
    xline(f_dom_emp, 'b--', sprintf('f_{dom}^{emp}=%.3g Hz', f_dom_emp), ...
        'LabelOrientation', 'aligned', 'LineWidth', 1.1);
end
if isfinite(f_dom_th)
    xline(f_dom_th, 'r--', sprintf('f_{dom}^{th}=%.3g Hz', f_dom_th), ...
        'LabelOrientation', 'aligned', 'LineWidth', 1.1);
end

xlabel('Frequence (Hz)');
ylabel('Amplitude');
title('Comparaison des spectres de p_{disp}^{emp}(t) et p_{disp}^{th}(t)');
legend('Spectre empirique', 'Spectre theorique', 'Location', 'best');

%% Figure 4 : spectres normalises pour comparer uniquement les pics
Amp_emp_norm = Amp_emp / max(Amp_emp(2:end), [], 'omitnan');
Amp_th_norm  = Amp_th  / max(Amp_th(2:end),  [], 'omitnan');

figure;
hold on; grid on; box on;
plot(freq, Amp_emp_norm, 'b-', 'LineWidth', 1.5);
plot(freq, Amp_th_norm,  'r-', 'LineWidth', 1.5);
if isfinite(f_dom_emp)
    xline(f_dom_emp, 'b--', 'f_{dom}^{emp}', 'LineWidth', 1.0);
end
if isfinite(f_dom_th)
    xline(f_dom_th, 'r--', 'f_{dom}^{th}', 'LineWidth', 1.0);
end
xlabel('Frequence (Hz)');
ylabel('Amplitude normalisee');
title('Spectres normalises de p_{disp}^{emp}(t) et p_{disp}^{th}(t)');
legend('Spectre empirique normalise', 'Spectre theorique normalise', ...
       'Location', 'best');

%% Sauvegarde
save('comparaison_pdisp_th_emp_spectres_results.mat', ...
     'tv', 'p_emp_v', 'p_emp_std_v', 'p_emp_global_v', ...
     'p_th_v', 'p_merge_v', 'p_break_v', ...
     'mean_emp', 'mean_emp_global', 'mean_th', 'idx_valid', ...
     't_uniform', 'freq', 'Amp_emp', 'Amp_th', ...
     'Amp_emp_norm', 'Amp_th_norm', 'f_dom_emp', 'f_dom_th');
