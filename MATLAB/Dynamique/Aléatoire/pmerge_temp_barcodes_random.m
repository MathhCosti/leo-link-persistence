clear; clc; close all;

%% p_merge(t) extrait des morts du barcode aux étapes G_k -> G_k U G_{k+1}
load('leo_H0_zigzag_random_vectors_barcodes.mat', ...
    'birth_index','death_index','h0_dims');
load('leo_zigzag_analysis_random_vectors_results.mat', ...
    'time_values','N','R','lambda','dmax','dt');

Nt = numel(time_values);
Nz = 2*Nt-1;
persist = (birth_index == 1) & (death_index == Nz);
birth_index(persist) = [];
death_index(persist) = [];

t = time_values(1:end-1).';
merge_count = zeros(Nt-1,1);
alive_before = zeros(Nt-1,1);
p_merge_emp = NaN(Nt-1,1);

for k = 1:Nt-1
    idx_G = 2*k-1;

    % Une classe qui meurt à G_k lors de la flèche vers l'union est absorbée
    % par une autre composante : c'est un événement de fusion H0.
    merge_count(k) = nnz(death_index == idx_G);

    % Nombre de barres du barcode présentes à l'indice G_k.
    alive_before(k) = nnz((birth_index <= idx_G) & (death_index >= idx_G));

    if alive_before(k) > 0
        p_merge_emp(k) = merge_count(k)/alive_before(k);
    end
end

[p_merge_th,~,~,chi] = random_model_probabilities(N,R,lambda,dmax,dt);
p_merge_global = sum(merge_count)/max(sum(alive_before),1);

figure;
plot(t,p_merge_emp,'o-','LineWidth',1.2,'MarkerSize',4); hold on; grid on;
yline(p_merge_th,'--','LineWidth',2, ...
    'Label',sprintf('théorie = %.4g',p_merge_th));
yline(p_merge_global,':','LineWidth',1.8, ...
    'Label',sprintf('moyenne empirique = %.4g',p_merge_global));
xlabel('Temps (s)');
ylabel('p_{merge}(t)');
title('p_{merge}(t) issue des barcodes H_0 — mouvement aléatoire');
legend('Empirique','Modèle théorique aléatoire','Moyenne empirique', ...
    'Location','best');
valid = [p_merge_emp(isfinite(p_merge_emp));p_merge_th;p_merge_global];
if ~isempty(valid), ylim([0,1.15*max(valid)]); end

fprintf('\n--- p_merge temporelle ---\n');
fprintf('p_merge théorique       : %.8f\n',p_merge_th);
fprintf('p_merge empirique global: %.8f\n',p_merge_global);
fprintf('Nombre total de fusions : %d\n',sum(merge_count));
fprintf('chi topologique         : %.8f\n',chi);

function [p_merge,p_break,p_disp,chi] = random_model_probabilities(N,R,lambda,dmax,dt)
    mu = 398600; omega = sqrt(mu/R^3);
    v_rel = (4/pi)*(R*omega);
    alpha_eff = 2*asin(min(dmax/(2*R),1));
    p_link = (1-cos(alpha_eff))/2;
    E_th = N*(N-1)/2*p_link;
    c2 = 1+3*sqrt(3)/(4*pi); c3c = 1+3*sqrt(3)/(2*pi); c3u = 1.80;
    N1 = N*max(1-p_link,0)^(N-1);
    N2 = 0; N3 = 0;
    if N>=2, N2=nchoosek(N,2)*p_link*max(1-c2*p_link,0)^(N-2); end
    if N>=3
        N3=nchoosek(N,3)*min(max(c3c*p_link^2,0),1)*max(1-c3u*p_link,0)^(N-3);
    end
    beta0_th=min(max(1+N1+N2+N3,1),N);
    if E_th>0, chi=min(max((N-beta0_th)/E_th,0),1); else, chi=0; end
    p_merge=1-exp(-lambda*(2*dmax*v_rel*dt)*chi);
    p_break=min(max((2/pi)*(v_rel*dt/dmax)*chi,0),1-eps);
    p_disp=1-(1-p_merge)*(1-p_break);
end
