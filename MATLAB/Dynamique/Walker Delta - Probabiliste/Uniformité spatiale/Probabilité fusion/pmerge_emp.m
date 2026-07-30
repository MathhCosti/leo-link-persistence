clear; clc; close all;

%% ============================================================
%  p_merge EMPIRIQUE A PARTIR DU BARCODE ZIGZAG H0
%
%  Zigzag :
%      G_k -> G_k U G_{k+1} <- G_{k+1}
%
%  Convention :
%      indice 2k-1 : G_k
%      indice 2k   : G_k U G_{k+1}
%
%  Une fusion entre t_k et t_{k+1} correspond a une barre qui meurt
%  a l'indice 2k-1.
%
%  Normalisation :
%      p_merge(k) = nb_fusions(k) / beta0(G_k)
%% ============================================================

%% Fichiers d'entree
script_dir = fileparts(mfilename('fullpath'));
barcode_file = fullfile(script_dir, '..', 'barcodes_results.mat');
analysis_file = fullfile(script_dir, '..', 'analysis_temp_results.mat');

if isempty(barcode_file)
    error(['Fichier barcode introuvable. Lancez d''abord le script ', ...
        'barcodes.m.']);
end

%% Chargement du barcode
Sbar = load(barcode_file, ...
    'birth_index', 'death_index', ...
    'ZigzagTime', 'ZigzagLabels', 'h0_dims');

death_index = Sbar.death_index(:);
ZigzagTime = Sbar.ZigzagTime(:);
h0_dims = Sbar.h0_dims(:);

Nz = numel(h0_dims);

if mod(Nz,2) ~= 1
    error('Nz doit etre impair : Nz = 2*Nt-1.');
end

Nt = (Nz+1)/2;
Ntrans = Nt-1;

%% Temps et parametres
if ~isempty(analysis_file)
    Sana = load(analysis_file, 'time_values', 'dt', 'inc_deg');
else
    Sana = struct();
end

if isfield(Sana,'time_values') && numel(Sana.time_values) == Nt
    time_values = double(Sana.time_values(:));
else
    time_values = ZigzagTime(1:2:end);
end

if isfield(Sana,'dt')
    dt = double(Sana.dt);
else
    dt = median(diff(time_values));
end

if isfield(Sana,'inc_deg')
    inc_deg = double(Sana.inc_deg);
else
    inc_deg = NaN;
end

%% Comptage des fusions
merge_count = zeros(Ntrans,1);
beta0_before = zeros(Ntrans,1);
beta0_union = zeros(Ntrans,1);

for k = 1:Ntrans
    idx_Gk = 2*k-1;
    idx_union = 2*k;

    merge_count(k) = sum(death_index == idx_Gk);

    beta0_before(k) = h0_dims(idx_Gk);
    beta0_union(k) = h0_dims(idx_union);
end

%% Verification topologique
merge_count_from_beta0 = beta0_before-beta0_union;
max_merge_error = max(abs(merge_count-merge_count_from_beta0));

fprintf('Erreur max comptage fusion barcode/beta0 : %g\n', ...
    max_merge_error);

if max_merge_error > 0
    warning(['Le comptage barcode ne coincide pas exactement avec ', ...
        'beta0(G_k)-beta0(union).']);
end

%% Probabilite empirique
p_merge = merge_count ./ max(beta0_before,1);

time_transition = ...
    0.5*(time_values(1:end-1)+time_values(2:end));

moving_window = min(15,Ntrans);
p_merge_moving = movmean( ...
    p_merge,moving_window,'Endpoints','shrink');

merge_count_moving = movmean( ...
    merge_count,moving_window,'Endpoints','shrink');

%% Moyennes globales
p_merge_mean = ...
    sum(merge_count)/max(sum(beta0_before),1);

p_merge_time_mean = mean(p_merge,'omitnan');

%% Figures
figure;
hold on; grid on;
plot(time_transition,p_merge,'LineWidth',0.8, ...
    'DisplayName','p_{merge}^{emp} instantane');
plot(time_transition,p_merge_moving,'LineWidth',2.2, ...
    'DisplayName',sprintf( ...
        'Moyenne glissante (%d transitions)',moving_window));
yline(p_merge_mean,'--', ...
    sprintf('Moyenne globale = %.4f',p_merge_mean));
xlabel('Temps au milieu de la transition (s)');
ylabel('p_{merge}^{emp}(t)');
title(sprintf( ...
    'Walker-Delta spatial : p_{merge}^{emp}(t), i = %.1f deg', ...
    inc_deg));
legend('Location','best');
hold off;

figure;
hold on; grid on;
stairs(time_transition,merge_count,'LineWidth',0.8, ...
    'DisplayName','Fusions instantanees');
plot(time_transition,merge_count_moving,'LineWidth',2, ...
    'DisplayName','Moyenne glissante');
xlabel('Temps au milieu de la transition (s)');
ylabel('Nombre de fusions');
title('Nombre de fusions H_0 par transition');
legend('Location','best');
hold off;

%% Console
fprintf('\n=== p_merge empirique ===\n');
fprintf('Fichier barcode                  : %s\n',barcode_file);
fprintf('Nombre de transitions            : %d\n',Ntrans);
fprintf('Nombre total de fusions          : %d\n',sum(merge_count));
fprintf('p_merge moyen pondere            : %.6f\n',p_merge_mean);
fprintf('Moyenne temporelle de p_merge(t) : %.6f\n', ...
    p_merge_time_mean);

%% Sauvegarde
save('pmerge_emp_results.mat', ...
    'time_values','time_transition','dt','inc_deg', ...
    'merge_count','merge_count_from_beta0', ...
    'beta0_before','beta0_union', ...
    'p_merge','moving_window','p_merge_moving', ...
    'merge_count_moving','p_merge_mean', ...
    'p_merge_time_mean','barcode_file','analysis_file');

fprintf('\nResultats sauvegardes dans pmerge_emp_results.mat\n');

function file = first_existing_file(candidates)
    file = '';
    for k = 1:numel(candidates)
        if isfile(candidates{k})
            file = candidates{k};
            return;
        end
    end
end
