clear; clc; close all;

%% p_break(t) extrait des naissances du barcode aux étapes G_k U G_{k+1} <- G_{k+1}
load('leo_H0_zigzag_random_vectors_barcodes.mat', ...
    'birth_index','death_index','h0_dims');
load('leo_zigzag_analysis_random_vectors_results.mat', ...
    'time_values','N','R','lambda','dmax','dt');

Nt = numel(time_values);
Nz = 2*Nt-1;
persist = (birth_index == 1) & (death_index == Nz);
birth_index(persist) = [];
death_index(persist) = [];

t = time_values(2:end).';
break_count = zeros(Nt-1,1);
alive_after = zeros(Nt-1,1);
p_break_emp = NaN(Nt-1,1);

for k = 1:Nt-1
    idx_next_G = 2*k+1;

    % Une classe née lors du passage de l'union à G_{k+1} correspond à une
    % nouvelle composante créée par fragmentation/rupture.
    break_count(k) = nnz(birth_index == idx_next_G);

    % Normalisation par les classes effectivement présentes dans G_{k+1}.
    alive_after(k) = nnz((birth_index <= idx_next_G) & ...
                         (death_index >= idx_next_G));

    if alive_after(k) > 0
        p_break_emp(k) = break_count(k)/alive_after(k);
    end
end

[~,p_break_th,~,chi] = random_model_probabilities(N,R,lambda,dmax,dt);
p_break_global = sum(break_count)/max(sum(alive_after),1);

figure;
plot(t,p_break_emp,'o-','LineWidth',1.2,'MarkerSize',4); hold on; grid on;
yline(p_break_th,'--','LineWidth',2, ...
    'Label',sprintf('théorie = %.4g',p_break_th));
yline(p_break_global,':','LineWidth',1.8, ...
    'Label',sprintf('moyenne empirique = %.4g',p_break_global));
xlabel('Temps (s)');
ylabel('p_{break}(t)');
title('p_{break}(t) issue des barcodes H_0 — mouvement aléatoire');
legend('Empirique','Modèle théorique aléatoire','Moyenne empirique', ...
    'Location','best');
valid = [p_break_emp(isfinite(p_break_emp));p_break_th;p_break_global];
if ~isempty(valid), ylim([0,1.15*max(valid)]); end

fprintf('\n--- p_break temporelle ---\n');
fprintf('p_break théorique       : %.8f\n',p_break_th);
fprintf('p_break empirique global: %.8f\n',p_break_global);
fprintf('Nombre total de ruptures: %d\n',sum(break_count));
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
