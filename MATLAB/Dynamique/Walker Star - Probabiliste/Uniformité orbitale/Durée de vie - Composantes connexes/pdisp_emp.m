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
analysis_file = 'leo_zigzag_analysis_results.mat';
barcode_file  = 'leo_H0_zigzag_barcodes_walker_star_deux_parts.mat';

if ~isfile(analysis_file)
    error('Fichier introuvable : %s. Lance d''abord analysis_temp_random_init.m.', ...
        analysis_file);
end

if ~isfile(barcode_file)
    error('Fichier introuvable : %s. Lance d''abord barcodes.m.', ...
        barcode_file);
end

Sanalysis = load(analysis_file, 'time_values', 'dt', 'inc_deg', 'R');
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


% Periode orbitale theorique
mu = 398600;  % parametre gravitationnel terrestre, km^3/s^2

if isfield(Sanalysis, 'R')
    R = Sanalysis.R;
else
    % Valeur par defaut correspondant a une altitude de 550 km
    R = 6371 + 550;
end

omega = sqrt(mu / R^3);
T_orb = 2*pi / omega;
T_half_orb = T_orb / 2;

f_orb = 1 / T_orb;
f_half_orb = 1 / T_half_orb;  % = 2*f_orb

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
smoothing_window = 20;
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

%% ============================================================
%  SPECTRE FREQUENTIEL DE LA MOYENNE GLISSANTE
%% ============================================================

% La FFT exige une serie sans NaN. Les eventuelles valeurs manquantes sont
% interpolees lineairement, puis la composante continue est retiree afin
% qu'elle ne masque pas les frequences non nulles.
signal_spectral = p_disp_smooth(:);

if any(~isfinite(signal_spectral))
    signal_spectral = fillmissing(signal_spectral, 'linear', ...
        'EndValues', 'nearest');
end

signal_spectral = signal_spectral - mean(signal_spectral);

Ns = numel(signal_spectral);
Fs = 1 / dt;                         % frequence d'echantillonnage en Hz

Y = fft(signal_spectral);
P2 = abs(Y / Ns);
P1 = P2(1:floor(Ns/2) + 1);

if Ns > 2
    P1(2:end-1) = 2 * P1(2:end-1);
end

frequencies = Fs * (0:floor(Ns/2))' / Ns;
periods = NaN(size(frequencies));
periods(frequencies > 0) = 1 ./ frequencies(frequencies > 0);

% Recherche du pic dominant hors composante continue.
if numel(P1) >= 2
    [dominant_amplitude, idx_peak_local] = max(P1(2:end));
    idx_peak = idx_peak_local + 1;
    dominant_frequency = frequencies(idx_peak);
    dominant_period = periods(idx_peak);
else
    dominant_amplitude = NaN;
    dominant_frequency = NaN;
    dominant_period = NaN;
end

figure;
plot(frequencies(2:end), P1(2:end), 'LineWidth', 1.7);
grid on;
xlabel('Frequence (Hz)');
ylabel('Amplitude');
title(sprintf(['Spectre de la moyenne mobile de p_{disp}^{emp}(t) ' ...
               '- fenetre = %d pas'], smoothing_window));

% Une echelle logarithmique rend les basses frequences plus lisibles.
if numel(frequencies) > 2 && all(frequencies(2:end) > 0)
    set(gca, 'XScale', 'log');
end

hold on;

if isfinite(dominant_frequency)
    xline(dominant_frequency, '--', ...
        sprintf('Pic : T = %.1f s', dominant_period), ...
        'LineWidth', 1.4, ...
        'LabelHorizontalAlignment', 'left');
end

% Frequence correspondant a la periode orbitale complete
xline(f_orb, '-.', ...
    sprintf('T_{orb} = %.1f s', T_orb), ...
    'LineWidth', 1.5, ...
    'LabelHorizontalAlignment', 'left');

% Frequence correspondant a la demi-periode orbitale
xline(f_half_orb, ':', ...
    sprintf('T_{orb}/2 = %.1f s', T_half_orb), ...
    'LineWidth', 1.7, ...
    'LabelHorizontalAlignment', 'left');

legend('Spectre de p_{disp}^{emp} lisse', ...
       'Pic dominant', ...
       'Frequence orbitale', ...
       'Frequence de demi-periode', ...
       'Location', 'best');

hold off;

fprintf('\n=== Analyse frequentielle de la moyenne mobile ===\n');
fprintf('Frequence dominante                     : %.6g Hz\n', dominant_frequency);
fprintf('Periode dominante                       : %.6g s\n', dominant_period);
fprintf('Amplitude du pic dominant               : %.6g\n', dominant_amplitude);
fprintf('Periode orbitale theorique              : %.6g s\n', T_orb);
fprintf('Frequence orbitale                      : %.6g Hz\n', f_orb);
fprintf('Demi-periode orbitale                   : %.6g s\n', T_half_orb);
fprintf('Frequence associee a la demi-periode    : %.6g Hz\n', f_half_orb);

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
    'smoothing_window', 'dt', 'inc_deg', ...
    'frequencies', 'periods', 'P1', ...
    'dominant_frequency', 'dominant_period', ...
    'dominant_amplitude', ...
    'R', 'mu', 'omega', 'T_orb', 'T_half_orb', ...
    'f_orb', 'f_half_orb');

fprintf('Resultats sauvegardes dans pdisp_emp_temp_delta_results.mat\n');
