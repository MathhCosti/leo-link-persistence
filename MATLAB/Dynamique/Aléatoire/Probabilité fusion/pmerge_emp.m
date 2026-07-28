clear; clc; close all;

%% ============================================================
%  p_merge EMPIRIQUE A PARTIR DU BARCODE ZIGZAG H0
%
%  Zigzag :
%      G_k -> U_k <- G_{k+1},   U_k = G_k union G_{k+1}
%
%  Une fusion correspond a une barre qui meurt a l'indice 2k-1.
%
%  Definition :
%      p_merge(k) = nb_fusions(k) / beta0(G_k)
%% ============================================================
script_dir = fileparts(mfilename('fullpath'));
analysis_file = fullfile(script_dir, '..', 'analysis_temp_results.mat');
barcode_file = fullfile(script_dir, '..', 'barcodes_results.mat');

if ~isfile(barcode_file)
    error('Fichier introuvable : %s', barcode_file);
end

Sbar = load(barcode_file, ...
    'death_index', 'ZigzagTime', 'h0_dims');

death_index = Sbar.death_index(:);
ZigzagTime  = Sbar.ZigzagTime(:);
h0_dims     = Sbar.h0_dims(:);

Nz = numel(h0_dims);

if mod(Nz,2) ~= 1
    error('Le zigzag doit contenir Nz = 2*Nt-1 objets.');
end

Nt = (Nz+1)/2;
Ntrans = Nt-1;

%% Temps physiques
if isfile(analysis_file)
    Sana = load(analysis_file, 'time_values', 'dt');
else
    Sana = struct();
end

if isfield(Sana,'time_values') && numel(Sana.time_values) == Nt
    time_values = Sana.time_values(:);
else
    time_values = ZigzagTime(1:2:end);
end

if isfield(Sana,'dt')
    dt = Sana.dt;
else
    dt = median(diff(time_values));
end

time_transition = 0.5*(time_values(1:end-1) + time_values(2:end));

%% Comptage des fusions
merge_count = zeros(Ntrans,1);
beta0_before = zeros(Ntrans,1);
beta0_union = zeros(Ntrans,1);

for k = 1:Ntrans
    idx_Gk    = 2*k-1;
    idx_union = 2*k;

    merge_count(k) = sum(death_index == idx_Gk);
    beta0_before(k) = h0_dims(idx_Gk);
    beta0_union(k)  = h0_dims(idx_union);
end

%% Verification topologique
merge_count_from_beta0 = beta0_before - beta0_union;
max_merge_error = max(abs(merge_count - merge_count_from_beta0));

fprintf('Erreur max comptage fusion barcode/beta0 : %g\n', ...
    max_merge_error);

if max_merge_error > 0
    warning(['Les comptages du barcode et les differences de beta0 ', ...
             'ne coincident pas exactement.']);
end

%% Probabilite empirique
p_merge = merge_count ./ max(beta0_before,1);

moving_window = min(15,Ntrans);
p_merge_moving = movmean(p_merge,moving_window,'Endpoints','shrink');
merge_count_moving = movmean(merge_count,moving_window,'Endpoints','shrink');

p_merge_mean = sum(merge_count) / max(sum(beta0_before),1);
p_merge_time_mean = mean(p_merge);

%% Traces
figure;
plot(time_transition,p_merge,'-','LineWidth',0.8); hold on;
plot(time_transition,p_merge_moving,'LineWidth',2.2);
yline(p_merge_mean,'--','LineWidth',1.5, ...
    'Label',sprintf('moyenne globale = %.4f',p_merge_mean));
grid on;
xlabel('Temps au milieu de la transition (s)');
ylabel('p_{merge}^{emp}(t)');
title('Probabilite empirique de fusion');
legend('Instantanee', ...
       sprintf('Moyenne glissante (%d points)',moving_window), ...
       'Moyenne globale','Location','best');
ylim([0,min(1,1.05*max([p_merge; p_merge_moving; p_merge_mean; eps]))]);
hold off;

figure;
stairs(time_transition,merge_count,'LineWidth',0.8); hold on;
plot(time_transition,merge_count_moving,'LineWidth',2);
grid on;
xlabel('Temps au milieu de la transition (s)');
ylabel('Nombre de fusions');
title('Fusions topologiques par transition');
legend('Comptage instantane','Moyenne glissante','Location','best');
hold off;

%% Affichage
fprintf('\n--- p_merge empirique ---\n');
fprintf('Nombre de transitions        : %d\n',Ntrans);
fprintf('Pas temporel moyen           : %.6f s\n',dt);
fprintf('Nombre total de fusions      : %d\n',sum(merge_count));
fprintf('Expositions dans G_k         : %d\n',sum(beta0_before));
fprintf('p_merge global pondere       : %.8f\n',p_merge_mean);
fprintf('moyenne temporelle simple    : %.8f\n',p_merge_time_mean);

%% Sauvegarde
output_file = 'pmerge_emp_results.mat';

save(output_file, ...
    'time_values','time_transition','dt', ...
    'merge_count','merge_count_from_beta0', ...
    'beta0_before','beta0_union', ...
    'p_merge','moving_window','p_merge_moving', ...
    'merge_count_moving','p_merge_mean','p_merge_time_mean');

fprintf('Resultats sauvegardes dans %s\n',output_file);
