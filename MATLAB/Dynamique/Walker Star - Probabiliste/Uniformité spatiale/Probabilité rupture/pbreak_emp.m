clear; clc; close all;

%% ============================================================
%  p_break EMPIRIQUE A PARTIR DU BARCODE ZIGZAG H0
%
%  Zigzag :
%      G_k -> G_k U G_{k+1} <- G_{k+1}
%
%  Convention :
%      indice 2k   : G_k U G_{k+1}
%      indice 2k+1 : G_{k+1}
%
%  Une rupture entre t_k et t_{k+1} correspond a une barre qui nait
%  a l'indice 2k+1.
%
%  Normalisation :
%      p_break(k) = nb_ruptures(k) / beta0(G_k U G_{k+1})
%
%  Le denominateur est donc le nombre de composantes de l'union,
%  exposees a une eventuelle separation lors du passage vers G_{k+1}.
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

birth_index = Sbar.birth_index(:);
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

%% Comptage des ruptures
break_count = zeros(Ntrans,1);
beta0_union = zeros(Ntrans,1);
beta0_after = zeros(Ntrans,1);

for k = 1:Ntrans
    idx_union = 2*k;
    idx_Gkp1 = 2*k+1;

    break_count(k) = sum(birth_index == idx_Gkp1);

    beta0_union(k) = h0_dims(idx_union);
    beta0_after(k) = h0_dims(idx_Gkp1);
end

%% Verification topologique
break_count_from_beta0 = beta0_after-beta0_union;
max_break_error = max(abs(break_count-break_count_from_beta0));

fprintf('Erreur max comptage rupture barcode/beta0 : %g\n', ...
    max_break_error);

if max_break_error > 0
    warning(['Le comptage barcode ne coincide pas exactement avec ', ...
        'beta0(G_{k+1})-beta0(union).']);
end

%% Probabilite empirique
p_break = break_count ./ max(beta0_union,1);

time_transition = ...
    0.5*(time_values(1:end-1)+time_values(2:end));

moving_window = min(15,Ntrans);
p_break_moving = movmean( ...
    p_break,moving_window,'Endpoints','shrink');

break_count_moving = movmean( ...
    break_count,moving_window,'Endpoints','shrink');

%% Moyennes globales
p_break_mean = ...
    sum(break_count)/max(sum(beta0_union),1);

p_break_time_mean = mean(p_break,'omitnan');

%% Figures
figure;
hold on; grid on;
plot(time_transition,p_break,'LineWidth',0.8, ...
    'DisplayName','p_{break}^{emp} instantane');
plot(time_transition,p_break_moving,'LineWidth',2.2, ...
    'DisplayName',sprintf( ...
        'Moyenne glissante (%d transitions)',moving_window));
yline(p_break_mean,'--', ...
    sprintf('Moyenne globale = %.4f',p_break_mean));
xlabel('Temps au milieu de la transition (s)');
ylabel('p_{break}^{emp}(t)');
title(sprintf( ...
    'Walker-Delta spatial : p_{break}^{emp}(t), i = %.1f deg', ...
    inc_deg));
legend('Location','best');
hold off;

figure;
hold on; grid on;
stairs(time_transition,break_count,'LineWidth',0.8, ...
    'DisplayName','Ruptures instantanees');
plot(time_transition,break_count_moving,'LineWidth',2, ...
    'DisplayName','Moyenne glissante');
xlabel('Temps au milieu de la transition (s)');
ylabel('Nombre de ruptures');
title('Nombre de ruptures H_0 par transition');
legend('Location','best');
hold off;

%% Console
fprintf('\n=== p_break empirique ===\n');
fprintf('Fichier barcode                  : %s\n',barcode_file);
fprintf('Nombre de transitions            : %d\n',Ntrans);
fprintf('Nombre total de ruptures         : %d\n',sum(break_count));
fprintf('p_break moyen pondere            : %.6f\n',p_break_mean);
fprintf('Moyenne temporelle de p_break(t) : %.6f\n', ...
    p_break_time_mean);

%% Sauvegarde
save('pbreak_emp_results.mat', ...
    'time_values','time_transition','dt','inc_deg', ...
    'break_count','break_count_from_beta0', ...
    'beta0_union','beta0_after', ...
    'p_break','moving_window','p_break_moving', ...
    'break_count_moving','p_break_mean', ...
    'p_break_time_mean','barcode_file','analysis_file');

fprintf('\nResultats sauvegardes dans pbreak_emp_results.mat\n');

function file = first_existing_file(candidates)
    file = '';
    for k = 1:numel(candidates)
        if isfile(candidates{k})
            file = candidates{k};
            return;
        end
    end
end
