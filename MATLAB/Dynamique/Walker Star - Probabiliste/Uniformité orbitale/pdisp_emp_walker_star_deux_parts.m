clear; clc; close all;

%% ============================================================
% p_disp^emp(t) DEPUIS LE BARCODE H0 - WALKER STAR DEUX PARTIES
%% ============================================================

barcode_file = 'leo_H0_zigzag_barcodes_walker_star_deux_parts.mat';
analysis_file = 'leo_zigzag_analysis_results.mat';
assert(isfile(barcode_file),'Fichier introuvable : %s',barcode_file);
assert(isfile(analysis_file),'Fichier introuvable : %s',analysis_file);

B = load(barcode_file,'birth_index','death_index','birth_time','death_time');
A = load(analysis_file,'time_values');

birth_index = B.birth_index(:);
death_index = B.death_index(:);
birth_time = B.birth_time(:);
death_time = B.death_time(:);
time_values = A.time_values(:);

Nt = numel(time_values);
t_emp = time_values(1:end-1);
alive_count = zeros(Nt-1,1);
death_count = zeros(Nt-1,1);
p_disp_emp = NaN(Nt-1,1);

for k=1:Nt-1
    t0=time_values(k);
    t1=time_values(k+1);

    % Classes vivantes au debut de l'intervalle [t_k,t_{k+1}].
    alive = (birth_time<=t0) & (death_time>t0);

    % Parmi ces classes, celles qui disparaissent avant ou a t_{k+1}.
    dying = alive & (death_time<=t1);

    alive_count(k)=nnz(alive);
    death_count(k)=nnz(dying);
    if alive_count(k)>0
        p_disp_emp(k)=death_count(k)/alive_count(k);
    end
end

p_disp_emp_global = sum(death_count)/max(sum(alive_count),1);

figure; hold on; grid on;
plot(t_emp,p_disp_emp,'o-','LineWidth',1.2,'MarkerSize',4, ...
    'DisplayName','p_{disp}^{emp}(t)');
yline(p_disp_emp_global,':','LineWidth',1.8, ...
    'DisplayName',sprintf('moyenne globale = %.4g',p_disp_emp_global));
xlabel('Temps (s)'); ylabel('p_{disp}^{emp}(t)');
title('Walker Star deux parties : probabilite empirique de disparition');
legend('Location','best');
valid=p_disp_emp(isfinite(p_disp_emp));
if ~isempty(valid)
    ymax=max([valid;p_disp_emp_global]);
    ylim([0,min(1,max(1e-12,1.15*ymax))]);
end
hold off;

fprintf('\n=== p_disp empirique Walker Star deux parties ===\n');
fprintf('Nombre total de disparitions : %d\n',sum(death_count));
fprintf('Nombre total d''expositions  : %d\n',sum(alive_count));
fprintf('p_disp empirique global      : %.8f\n',p_disp_emp_global);

save('pdisp_emp_walker_star_deux_parts_results.mat', ...
    't_emp','p_disp_emp','death_count','alive_count','p_disp_emp_global');

window = 15;
pdisp_smooth = movmean(p_disp_emp, window, 'omitnan');

figure;
plot(time_values(1:end-1), p_disp_emp, ...
    'Color', [0.7 0.7 0.7]); hold on;
plot(time_values(1:end-1), pdisp_smooth, ...
    'LineWidth', 2);
yline(mean(p_disp_emp, 'omitnan'), '--', 'Moyenne');
grid on;
xlabel('Temps (s)');
ylabel('p_{\mathrm{disp}}^{emp}(t)');
legend('Empirique', 'Moyenne glissante', 'Moyenne globale');