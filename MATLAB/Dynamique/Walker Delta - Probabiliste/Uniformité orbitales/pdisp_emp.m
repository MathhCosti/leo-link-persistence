clear; clc; close all;

%% ============================================================
%  p_disp EMPIRIQUE EN FONCTION DU TEMPS
%  Orbites aleatoires a inclinaison fixe
%
%  Definition utilisee sur chaque intervalle [t_k,t_{k+1}] :
%
%      p_disp(t_k) = nombre de barres H0 vivantes a t_k
%                    qui meurent avant ou a t_{k+1}
%                    ---------------------------------------
%                    nombre de barres H0 vivantes a t_k
%
%  Les barres qui atteignent la fin de la simulation sont traitees
%  comme censurees a droite et ne sont pas comptees comme des morts.
%% ============================================================

%% Fichiers d'entree
analysis_file = 'leo_zigzag_analysis_results_delta.mat';
barcode_file  = 'leo_H0_zigzag_barcodes_delta.mat';

if ~isfile(analysis_file)
    error('Fichier introuvable : %s. Lance d''abord analysis_temp_random_init.m.', ...
        analysis_file);
end

if ~isfile(barcode_file)
    error('Fichier introuvable : %s. Lance d''abord barcodes.m.', ...
        barcode_file);
end

Sanalysis = load(analysis_file, 'time_values', 'dt', 'inc_deg');
Sbarcode  = load(barcode_file, 'birth_time', 'death_time', 'lifetimes');

time_values = Sanalysis.time_values(:);
dt = Sanalysis.dt;

birth_time = Sbarcode.birth_time(:);
death_time = Sbarcode.death_time(:);

if isfield(Sanalysis, 'inc_deg')
    inc_deg = Sanalysis.inc_deg;
else
    inc_deg = NaN;
end

%% Verification
if numel(time_values) < 2
    error('La grille temporelle doit contenir au moins deux instants.');
end

T_end = time_values(end);
tol = 1e-10 * max(1, abs(T_end));

%% ============================================================
%  Calcul de p_disp(t)
%% ============================================================

Nt = numel(time_values);
t_pdisp = time_values(1:end-1);

alive_counts = zeros(Nt-1,1);
death_counts = zeros(Nt-1,1);
p_disp_emp_t = NaN(Nt-1,1);

% Une barre qui se termine au dernier instant peut simplement etre encore
% vivante a la fin de la fenetre d'observation. On la considere censuree.
is_right_censored = abs(death_time - T_end) <= tol;

for k = 1:Nt-1
    t0 = time_values(k);
    t1 = time_values(k+1);

    % Barres deja nees et encore vivantes au debut de l'intervalle.
    alive = (birth_time <= t0 + tol) & (death_time > t0 + tol);

    % Parmi ces barres, morts observees pendant l'intervalle.
    dying = alive & ...
            (death_time <= t1 + tol) & ...
            ~is_right_censored;

    alive_counts(k) = sum(alive);
    death_counts(k) = sum(dying);

    if alive_counts(k) > 0
        p_disp_emp_t(k) = death_counts(k) / alive_counts(k);
    end
end

%% Moyennes
p_disp_mean_unweighted = mean(p_disp_emp_t, 'omitnan');

if sum(alive_counts) > 0
    % Estimateur global pondere par le nombre de barres exposees au risque.
    p_disp_global = sum(death_counts) / sum(alive_counts);
else
    p_disp_global = NaN;
end

% Lissage uniquement pour rendre la tendance temporelle lisible.
smoothing_window = 5;
p_disp_smooth = movmean(p_disp_emp_t, smoothing_window, 'omitnan');

%% ============================================================
%  Affichage
%% ============================================================

figure;
hold on;
grid on;

plot(t_pdisp, p_disp_emp_t, 'o-', ...
    'LineWidth', 1.0, 'MarkerSize', 3, ...
    'DisplayName', 'p_{disp} empirique brut');

plot(t_pdisp, p_disp_smooth, '-', ...
    'LineWidth', 2.0, ...
    'DisplayName', sprintf('Moyenne mobile (%d pas)', smoothing_window));

yline(p_disp_global, '--', ...
    sprintf('Moyenne globale ponderee = %.4g', p_disp_global), ...
    'LineWidth', 1.7, ...
    'LabelHorizontalAlignment', 'left', ...
    'DisplayName', 'Moyenne globale ponderee');

xlabel('Temps t_k (s)');
ylabel('p_{disp}(t_k)');
title(sprintf(['Probabilite empirique de disparition des barres H_0 ' ...
               '- inclinaison %.1f deg'], inc_deg));
legend('Location', 'best');

valid_values = [p_disp_emp_t(isfinite(p_disp_emp_t)); p_disp_smooth(isfinite(p_disp_smooth))];
if ~isempty(valid_values)
    ymax = max(valid_values);
    ylim([0, max(0.01, 1.15*ymax)]);
end

hold off;

%% Graphe des effectifs servant au calcul
figure;
plot(t_pdisp, alive_counts, 'LineWidth', 1.5); hold on;
plot(t_pdisp, death_counts, 'LineWidth', 1.5);
grid on;
xlabel('Temps t_k (s)');
ylabel('Nombre de barres');
title('Barres exposees au risque et disparitions observees');
legend('Barres vivantes a t_k', 'Disparitions sur [t_k,t_{k+1}]', ...
    'Location', 'best');
hold off;

%% Console
fprintf('\n=== p_disp empirique temporel ===\n');
fprintf('Pas temporel dt                         : %.2f s\n', dt);
fprintf('Nombre total de barres H0               : %d\n', numel(birth_time));
fprintf('Barres censurees a droite               : %d\n', sum(is_right_censored));
fprintf('Moyenne temporelle non ponderee         : %.6g\n', p_disp_mean_unweighted);
fprintf('Moyenne globale ponderee                : %.6g\n', p_disp_global);
fprintf('Nombre total de disparitions observees  : %d\n', sum(death_counts));
fprintf('Nombre total d''expositions             : %d\n', sum(alive_counts));

%% Sauvegarde
save('pdisp_emp_temp_delta_results.mat', ...
    't_pdisp', 'p_disp_emp_t', 'p_disp_smooth', ...
    'alive_counts', 'death_counts', ...
    'p_disp_mean_unweighted', 'p_disp_global', ...
    'smoothing_window', 'dt', 'inc_deg');

fprintf('Resultats sauvegardes dans pdisp_emp_temp_delta_results.mat\n');
