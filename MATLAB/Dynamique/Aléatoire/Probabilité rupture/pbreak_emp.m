clear; clc; close all;

%% ============================================================
%  p_break EMPIRIQUE A PARTIR DU BARCODE ZIGZAG H0
%
%  Zigzag :
%      G_k -> U_k <- G_{k+1},   U_k = G_k union G_{k+1}
%
%  Une rupture correspond a une barre qui nait a l'indice 2k+1.
%
%  Definition :
%      p_break(k) = nb_ruptures(k) / beta0(U_k)
%% ============================================================
script_dir = fileparts(mfilename('fullpath'));
analysis_file = fullfile(script_dir, '..', 'analysis_temp_results.mat');
barcode_file = fullfile(script_dir, '..', 'barcodes_results.mat');

if ~isfile(barcode_file)
    error('Fichier introuvable : %s', barcode_file);
end

Sbar = load(barcode_file, ...
    'birth_index', 'ZigzagTime', 'h0_dims');

birth_index = Sbar.birth_index(:);
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

%% Comptage des ruptures
break_count = zeros(Ntrans,1);
beta0_union = zeros(Ntrans,1);
beta0_after = zeros(Ntrans,1);

for k = 1:Ntrans
    idx_union = 2*k;
    idx_Gkp1  = 2*k+1;

    break_count(k) = sum(birth_index == idx_Gkp1);
    beta0_union(k) = h0_dims(idx_union);
    beta0_after(k) = h0_dims(idx_Gkp1);
end

%% Verification topologique
break_count_from_beta0 = beta0_after - beta0_union;
max_break_error = max(abs(break_count - break_count_from_beta0));

fprintf('Erreur max comptage rupture barcode/beta0 : %g\n', ...
    max_break_error);

if max_break_error > 0
    warning(['Les comptages du barcode et les differences de beta0 ', ...
             'ne coincident pas exactement.']);
end

%% Probabilite empirique
p_break = break_count ./ max(beta0_union,1);

moving_window = min(15,Ntrans);
p_break_moving = movmean(p_break,moving_window,'Endpoints','shrink');
break_count_moving = movmean(break_count,moving_window,'Endpoints','shrink');

p_break_mean = sum(break_count) / max(sum(beta0_union),1);
p_break_time_mean = mean(p_break);

%% Traces
figure;
plot(time_transition,p_break,'-','LineWidth',0.8); hold on;
plot(time_transition,p_break_moving,'LineWidth',2.2);
yline(p_break_mean,'--','LineWidth',1.5, ...
    'Label',sprintf('moyenne globale = %.4f',p_break_mean));
grid on;
xlabel('Temps au milieu de la transition (s)');
ylabel('p_{break}^{emp}(t)');
title('Probabilite empirique de rupture');
legend('Instantanee', ...
       sprintf('Moyenne glissante (%d points)',moving_window), ...
       'Moyenne globale','Location','best');
ylim([0,min(1,1.05*max([p_break; p_break_moving; p_break_mean; eps]))]);
hold off;

figure;
stairs(time_transition,break_count,'LineWidth',0.8); hold on;
plot(time_transition,break_count_moving,'LineWidth',2);
grid on;
xlabel('Temps au milieu de la transition (s)');
ylabel('Nombre de ruptures');
title('Ruptures topologiques par transition');
legend('Comptage instantane','Moyenne glissante','Location','best');
hold off;

%% Affichage
fprintf('\n--- p_break empirique ---\n');
fprintf('Nombre de transitions        : %d\n',Ntrans);
fprintf('Pas temporel moyen           : %.6f s\n',dt);
fprintf('Nombre total de ruptures     : %d\n',sum(break_count));
fprintf('Expositions dans les unions  : %d\n',sum(beta0_union));
fprintf('p_break global pondere       : %.8f\n',p_break_mean);
fprintf('moyenne temporelle simple    : %.8f\n',p_break_time_mean);

%% Sauvegarde
output_file = 'pbreak_emp_results.mat';

save(output_file, ...
    'time_values','time_transition','dt', ...
    'break_count','break_count_from_beta0', ...
    'beta0_union','beta0_after', ...
    'p_break','moving_window','p_break_moving', ...
    'break_count_moving','p_break_mean','p_break_time_mean');

fprintf('Resultats sauvegardes dans %s\n',output_file);
