%% comparaison_pdisp_delta_spatial.m
% Comparaison de p_disp théorique et empirique pour le Walker-Delta
% à uniformité spatiale, avec comparaison fréquentielle.

clear; clc; close all;

script_dir = fileparts(mfilename('fullpath'));

th_candidates = { ...
    fullfile(script_dir,'pdisp_modele_delta_spatial_results.mat'), ...
    fullfile(script_dir,'..','pdisp_modele_delta_spatial_results.mat')};

emp_candidates = { ...
    fullfile(script_dir,'pdisp_emp_delta_spatial_results.mat'), ...
    fullfile(script_dir,'..','pdisp_emp_delta_spatial_results.mat')};

th_file = first_existing_file(th_candidates);
emp_file = first_existing_file(emp_candidates);

if isempty(th_file)
    error('pdisp_modele_delta_spatial_results.mat introuvable.');
end

if isempty(emp_file)
    error('pdisp_emp_delta_spatial_results.mat introuvable.');
end

Sth = load(th_file);
Semp = load(emp_file);

required_th = {'t_transition','p_disp_t','p_merge_t','p_break_t'};
required_emp = {'t_pdisp','p_disp_emp_t','p_disp_smooth'};

for q = 1:numel(required_th)
    if ~isfield(Sth,required_th{q})
        error('Variable %s absente du fichier théorique.',required_th{q});
    end
end

for q = 1:numel(required_emp)
    if ~isfield(Semp,required_emp{q})
        error('Variable %s absente du fichier empirique.',required_emp{q});
    end
end

t_th = double(Sth.t_transition(:));
p_th = double(Sth.p_disp_t(:));
p_merge = double(Sth.p_merge_t(:));
p_break = double(Sth.p_break_t(:));

t_emp = double(Semp.t_pdisp(:));
p_emp_raw = double(Semp.p_disp_emp_t(:));
p_emp_smooth = double(Semp.p_disp_smooth(:));

%% Interpolation de la théorie
p_th_i = interp1(t_th,p_th,t_emp,'linear','extrap');
p_merge_i = interp1(t_th,p_merge,t_emp,'linear','extrap');
p_break_i = interp1(t_th,p_break,t_emp,'linear','extrap');

%% Suppression des bords
edge_trim_steps = 2;

idx_valid = true(size(t_emp));

if numel(idx_valid)>2*edge_trim_steps
    idx_valid(1:edge_trim_steps) = false;
    idx_valid(end-edge_trim_steps+1:end) = false;
end

idx_valid = idx_valid & ...
    isfinite(p_emp_raw) & ...
    isfinite(p_emp_smooth) & ...
    isfinite(p_th_i);

tv = t_emp(idx_valid);
p_emp_v = p_emp_smooth(idx_valid);
p_emp_raw_v = p_emp_raw(idx_valid);
p_th_v = p_th_i(idx_valid);
p_merge_v = p_merge_i(idx_valid);
p_break_v = p_break_i(idx_valid);

mean_emp = mean(p_emp_v,'omitnan');
mean_th = mean(p_th_v,'omitnan');

rmse = sqrt(mean((p_emp_v-p_th_v).^2,'omitnan'));
mae = mean(abs(p_emp_v-p_th_v),'omitnan');

fprintf('\n=== Comparaison Delta spatial ===\n');
fprintf('Points conservés : %d / %d\n',nnz(idx_valid),numel(idx_valid));
fprintf('Moyenne empirique : %.6f\n',mean_emp);
fprintf('Moyenne théorique : %.6f\n',mean_th);
fprintf('RMSE              : %.6f\n',rmse);
fprintf('MAE               : %.6f\n',mae);

%% Courbe principale
figure;
hold on; grid on; box on;
plot(tv,p_emp_raw_v,'Color',[0.7 0.7 0.7], ...
    'LineWidth',0.8,'DisplayName','Empirique brut');
plot(tv,p_emp_v,'k','LineWidth',1.8, ...
    'DisplayName','Empirique lissé');
plot(tv,p_th_v,'--','LineWidth',2, ...
    'DisplayName','Théorie Delta spatial');
yline(mean_emp,':',sprintf('Moyenne emp = %.3f',mean_emp));
yline(mean_th,':',sprintf('Moyenne th = %.3f',mean_th));
xlabel('Temps (s)');
ylabel('p_{disp}(t)');
title('Walker-Delta spatial : comparaison théorique / empirique');
legend('Location','best');

%% Décomposition théorique
figure;
hold on; grid on; box on;
plot(tv,p_merge_v,'LineWidth',1.4);
plot(tv,p_break_v,'LineWidth',1.4);
plot(tv,p_th_v,'LineWidth',1.8);
plot(tv,p_emp_v,'k','LineWidth',1.2);
xlabel('Temps (s)');
ylabel('Probabilité par pas');
title('Décomposition du modèle et comparaison empirique');
legend('p_{merge}^{th}','p_{break}^{th}', ...
    'p_{disp}^{th}','p_{disp}^{emp}', ...
    'Location','best');

%% Spectres
if numel(tv)<4
    error('Pas assez de points pour calculer les spectres.');
end

dt_spec = median(diff(tv));
t_uniform = (tv(1):dt_spec:tv(end)).';

x_emp = interp1(tv,p_emp_v,t_uniform,'linear');
x_th = interp1(tv,p_th_v,t_uniform,'linear');

x_emp = x_emp-mean(x_emp,'omitnan');
x_th = x_th-mean(x_th,'omitnan');

Nspec = numel(t_uniform);
n = (0:Nspec-1).';
win = 0.5*(1-cos(2*pi*n/max(Nspec-1,1)));

Y_emp = fft(x_emp.*win);
Y_th = fft(x_th.*win);

P_emp = abs(Y_emp/Nspec);
P_th = abs(Y_th/Nspec);

n_half = floor(Nspec/2)+1;
freq = (1/dt_spec)*(0:n_half-1).'/Nspec;

Amp_emp = 2*P_emp(1:n_half);
Amp_th = 2*P_th(1:n_half);
Amp_emp(1) = P_emp(1);
Amp_th(1) = P_th(1);

if numel(freq)>=2
    [~,ie] = max(Amp_emp(2:end));
    [~,it] = max(Amp_th(2:end));
    ie = ie+1;
    it = it+1;
    f_dom_emp = freq(ie);
    f_dom_th = freq(it);
else
    f_dom_emp = NaN;
    f_dom_th = NaN;
end

fprintf('\nFréquence dominante empirique : %.6g Hz, T = %.2f s\n', ...
    f_dom_emp,1/f_dom_emp);
fprintf('Fréquence dominante théorique : %.6g Hz, T = %.2f s\n', ...
    f_dom_th,1/f_dom_th);

figure;
hold on; grid on; box on;
plot(freq,Amp_emp,'LineWidth',1.5);
plot(freq,Amp_th,'LineWidth',1.5);
if isfinite(f_dom_emp)
    xline(f_dom_emp,'--',sprintf('T_{emp}=%.0f s',1/f_dom_emp));
end
if isfinite(f_dom_th)
    xline(f_dom_th,'--',sprintf('T_{th}=%.0f s',1/f_dom_th));
end
xlabel('Fréquence (Hz)');
ylabel('Amplitude');
title('Spectres de p_{disp}^{emp} et p_{disp}^{th}');
legend('Empirique','Théorique','Location','best');

%% Sauvegarde
save('comparaison_pdisp_delta_spatial_results.mat', ...
    'tv','p_emp_raw_v','p_emp_v','p_th_v', ...
    'p_merge_v','p_break_v','mean_emp','mean_th', ...
    'rmse','mae','idx_valid', ...
    't_uniform','freq','Amp_emp','Amp_th', ...
    'f_dom_emp','f_dom_th','th_file','emp_file');

fprintf('\nRésultats sauvegardés dans ');
fprintf('comparaison_pdisp_delta_spatial_results.mat\n');

function file = first_existing_file(candidates)
    file = '';
    for k = 1:numel(candidates)
        if isfile(candidates{k})
            file = candidates{k};
            return;
        end
    end
end
