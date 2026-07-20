clear; clc; close all;

%% ============================================================
%  p_disp EMPIRIQUE — WALKER DELTA A UNIFORMITE SPATIALE
%
%  Sur chaque intervalle [t_k,t_{k+1}] :
%
%      p_disp^emp(t_k)
%      = nombre de barres H0 vivantes à t_k qui meurent sur l'intervalle
%        ---------------------------------------------------------------
%        nombre de barres H0 vivantes à t_k
%
%  Les barres atteignant la fin de la fenêtre sont censurées à droite.
%% ============================================================

script_dir = fileparts(mfilename('fullpath'));

%% Candidats pour les fichiers Delta
barcode_candidates = { ...
    fullfile(script_dir,'leo_H0_zigzag_barcodes_delta_spatial.mat'), ...
    fullfile(script_dir,'leo_H0_zigzag_barcodes_delta.mat'), ...
    fullfile(script_dir,'leo_H0_zigzag_barcodes.mat'), ...
    fullfile(script_dir,'..','leo_H0_zigzag_barcodes_delta_spatial.mat'), ...
    fullfile(script_dir,'..','leo_H0_zigzag_barcodes_delta.mat')};

analysis_candidates = { ...
    fullfile(script_dir,'leo_zigzag_analysis_results_delta.mat'), ...
    fullfile(script_dir,'leo_zigzag_analysis_results_delta_spatial.mat'), ...
    fullfile(script_dir,'leo_zigzag_analysis_results.mat'), ...
    fullfile(script_dir,'..','leo_zigzag_analysis_results_delta.mat'), ...
    fullfile(script_dir,'..','leo_zigzag_analysis_results_delta_spatial.mat')};

barcode_file = first_existing_file(barcode_candidates);
analysis_file = first_existing_file(analysis_candidates);

if isempty(barcode_file)
    error('Aucun fichier barcode Delta trouvé.');
end

if isempty(analysis_file)
    error('Aucun fichier d''analyse Delta trouvé.');
end

fprintf('Barcode chargé : %s\n',barcode_file);
fprintf('Analyse chargée : %s\n',analysis_file);

Sanalysis = load(analysis_file);
Sbarcode = load(barcode_file);

%% Variables temporelles
if ~isfield(Sanalysis,'time_values')
    error('time_values absent de %s.',analysis_file);
end

time_values = double(Sanalysis.time_values(:));

if isfield(Sanalysis,'dt')
    dt = double(Sanalysis.dt);
else
    dt = median(diff(time_values));
end

if isfield(Sanalysis,'inc_deg')
    inc_deg = double(Sanalysis.inc_deg);
else
    inc_deg = NaN;
end

if isfield(Sanalysis,'R')
    R = double(Sanalysis.R);
else
    R = 6371+550;
end

%% Barcodes
required_barcode = {'birth_time','death_time'};

for q = 1:numel(required_barcode)
    if ~isfield(Sbarcode,required_barcode{q})
        error('Variable %s absente de %s.', ...
            required_barcode{q},barcode_file);
    end
end

birth_time = double(Sbarcode.birth_time(:));
death_time = double(Sbarcode.death_time(:));

%% Périodes orbitales
mu = 398600;
omega = sqrt(mu/R^3);
T_orb = 2*pi/omega;
T_half_orb = T_orb/2;

f_orb = 1/T_orb;
f_half_orb = 1/T_half_orb;

%% Calcul
Nt = numel(time_values);

if Nt<2
    error('La grille temporelle doit contenir au moins deux instants.');
end

t_pdisp = time_values(1:end-1);
T_end = time_values(end);
tol = 1e-10*max(1,abs(T_end));

alive_counts = zeros(Nt-1,1);
death_counts = zeros(Nt-1,1);
p_disp_emp_t = NaN(Nt-1,1);

is_right_censored = abs(death_time-T_end)<=tol;

for k = 1:Nt-1
    t0 = time_values(k);
    t1 = time_values(k+1);

    alive = birth_time<=t0+tol & death_time>t0+tol;

    dying = alive & ...
        death_time<=t1+tol & ...
        ~is_right_censored;

    alive_counts(k) = nnz(alive);
    death_counts(k) = nnz(dying);

    if alive_counts(k)>0
        p_disp_emp_t(k) = death_counts(k)/alive_counts(k);
    end
end

p_disp_mean_unweighted = mean(p_disp_emp_t,'omitnan');

if sum(alive_counts)>0
    p_disp_global = sum(death_counts)/sum(alive_counts);
else
    p_disp_global = NaN;
end

smoothing_window = 15;
p_disp_smooth = movmean(p_disp_emp_t,smoothing_window,'omitnan');

%% Figures
figure;
hold on; grid on;
plot(t_pdisp,p_disp_emp_t,'o-','LineWidth',0.9, ...
    'MarkerSize',3,'DisplayName','p_{disp}^{emp} brut');
plot(t_pdisp,p_disp_smooth,'LineWidth',2, ...
    'DisplayName',sprintf('Moyenne glissante (%d pas)',smoothing_window));
yline(p_disp_global,'--', ...
    sprintf('Moyenne globale = %.4g',p_disp_global), ...
    'DisplayName','Moyenne globale pondérée');
xlabel('Temps (s)');
ylabel('p_{disp}^{emp}(t)');
title(sprintf(['Walker-Delta spatial : disparition empirique ', ...
    'des barres H_0 — i = %.1f deg'],inc_deg));
legend('Location','best');

figure;
hold on; grid on;
plot(t_pdisp,alive_counts,'LineWidth',1.5);
plot(t_pdisp,death_counts,'LineWidth',1.5);
xlabel('Temps (s)');
ylabel('Nombre de barres');
title('Barres exposées au risque et disparitions');
legend('Barres vivantes','Disparitions','Location','best');

%% Spectre
signal_spectral = p_disp_smooth(:);

if any(~isfinite(signal_spectral))
    signal_spectral = fillmissing(signal_spectral,'linear', ...
        'EndValues','nearest');
end

signal_spectral = signal_spectral-mean(signal_spectral);

Ns = numel(signal_spectral);
Fs = 1/dt;

Y = fft(signal_spectral);
P2 = abs(Y/Ns);
P1 = P2(1:floor(Ns/2)+1);

if Ns>2
    P1(2:end-1) = 2*P1(2:end-1);
end

frequencies = Fs*(0:floor(Ns/2)).'/Ns;
periods = NaN(size(frequencies));
periods(frequencies>0) = 1./frequencies(frequencies>0);

if numel(P1)>=2
    [dominant_amplitude,idx_rel] = max(P1(2:end));
    idx_peak = idx_rel+1;
    dominant_frequency = frequencies(idx_peak);
    dominant_period = periods(idx_peak);
else
    dominant_amplitude = NaN;
    dominant_frequency = NaN;
    dominant_period = NaN;
end

figure;
hold on; grid on;
plot(frequencies(2:end),P1(2:end),'LineWidth',1.7);
xline(f_orb,'-.',sprintf('T_{orb}=%.0f s',T_orb));
xline(f_half_orb,':',sprintf('T_{orb}/2=%.0f s',T_half_orb));
if isfinite(dominant_frequency)
    xline(dominant_frequency,'--', ...
        sprintf('Pic T=%.0f s',dominant_period));
end
xlabel('Fréquence (Hz)');
ylabel('Amplitude');
title('Spectre de p_{disp}^{emp}(t) lissé');
legend('Spectre','Fréquence orbitale','Demi-période', ...
    'Pic dominant','Location','best');

%% Console
fprintf('\n=== p_disp empirique Delta spatial ===\n');
fprintf('dt                              : %.2f s\n',dt);
fprintf('Nombre de barres                : %d\n',numel(birth_time));
fprintf('Barres censurées à droite       : %d\n',nnz(is_right_censored));
fprintf('Moyenne temporelle              : %.6g\n',p_disp_mean_unweighted);
fprintf('Moyenne globale pondérée        : %.6g\n',p_disp_global);
fprintf('Période dominante               : %.2f s\n',dominant_period);
fprintf('Demi-période orbitale théorique : %.2f s\n',T_half_orb);

%% Sauvegarde
save('pdisp_emp_delta_spatial_results.mat', ...
    't_pdisp','p_disp_emp_t','p_disp_smooth', ...
    'alive_counts','death_counts', ...
    'p_disp_mean_unweighted','p_disp_global', ...
    'smoothing_window','dt','inc_deg', ...
    'frequencies','periods','P1', ...
    'dominant_frequency','dominant_period','dominant_amplitude', ...
    'R','mu','omega','T_orb','T_half_orb','f_orb','f_half_orb', ...
    'barcode_file','analysis_file');

fprintf('\nRésultats sauvegardés dans ');
fprintf('pdisp_emp_delta_spatial_results.mat\n');

function file = first_existing_file(candidates)
    file = '';
    for k = 1:numel(candidates)
        if isfile(candidates{k})
            file = candidates{k};
            return;
        end
    end
end
