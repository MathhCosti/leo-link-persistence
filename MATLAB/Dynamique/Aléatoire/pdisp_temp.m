clear; clc; close all;

%% p_disp(t) extrait des barcodes H0 et comparaison au modèle aléatoire
load('leo_H0_zigzag_random_vectors_barcodes.mat', ...
    'birth_index','death_index','birth_time','death_time','lifetimes');
load('leo_zigzag_analysis_random_vectors_results.mat', ...
    'time_values','N','R','lambda','dmax','dt');

% La composante globale persistante doit être retirée si elle est encore présente.
Nz = 2*numel(time_values)-1;
persist = (birth_index == 1) & (death_index == Nz);
birth_index(persist) = [];
death_index(persist) = [];
birth_time(persist) = [];
death_time(persist) = [];
lifetimes(persist) = [];

Nt = numel(time_values);
t = time_values(1:end-1).';
alive_count = zeros(Nt-1,1);
death_count = zeros(Nt-1,1);
p_disp_emp = NaN(Nt-1,1);

for k = 1:Nt-1
    t0 = time_values(k);
    t1 = time_values(k+1);

    alive = (birth_time <= t0) & (death_time > t0);
    dying = alive & (death_time <= t1);

    alive_count(k) = nnz(alive);
    death_count(k) = nnz(dying);

    if alive_count(k) > 0
        p_disp_emp(k) = death_count(k)/alive_count(k);
    end
end

[p_merge_th,p_break_th,p_disp_th,chi] = random_model_probabilities(N,R,lambda,dmax,dt);

p_disp_global = sum(death_count)/max(sum(alive_count),1);

figure;
plot(t,p_disp_emp,'o-','LineWidth',1.2,'MarkerSize',4); hold on; grid on;
yline(p_disp_th,'--','LineWidth',2, ...
    'Label',sprintf('théorie = %.4g',p_disp_th));
yline(p_disp_global,':','LineWidth',1.8, ...
    'Label',sprintf('moyenne empirique = %.4g',p_disp_global));
xlabel('Temps (s)');
ylabel('p_{disp}(t)');
title('p_{disp}(t) issue des barcodes H_0 — mouvement aléatoire');
legend('Empirique','Modèle théorique aléatoire','Moyenne empirique', ...
    'Location','best');
valid = [p_disp_emp(isfinite(p_disp_emp));p_disp_th;p_disp_global];
if ~isempty(valid), ylim([0,1.15*max(valid)]); end

fprintf('\n--- p_disp temporelle ---\n');
fprintf('p_disp théorique       : %.8f\n',p_disp_th);
fprintf('p_disp empirique global: %.8f\n',p_disp_global);
fprintf('p_merge théorique      : %.8f\n',p_merge_th);
fprintf('p_break théorique      : %.8f\n',p_break_th);
fprintf('chi topologique        : %.8f\n',chi);

function [p_merge,p_break,p_disp,chi] = random_model_probabilities(N,R,lambda,dmax,dt)
    mu = 398600;
    omega = sqrt(mu/R^3);
    v_orb = R*omega;
    v_rel = (4/pi)*v_orb;

    alpha_eff = 2*asin(min(dmax/(2*R),1));
    p_link = (1-cos(alpha_eff))/2;
    E_th = N*(N-1)/2*p_link;

    c2_union = 1 + 3*sqrt(3)/(4*pi);
    c3_conn = 1 + 3*sqrt(3)/(2*pi);
    c3_union = 1.80;

    N1 = N*max(1-p_link,0)^(N-1);
    if N >= 2
        N2 = nchoosek(N,2)*p_link*max(1-c2_union*p_link,0)^(N-2);
    else
        N2 = 0;
    end
    if N >= 3
        p_conn3 = min(max(c3_conn*p_link^2,0),1);
        N3 = nchoosek(N,3)*p_conn3*max(1-c3_union*p_link,0)^(N-3);
    else
        N3 = 0;
    end

    beta0_th = min(max(1+N1+N2+N3,1),N);
    if E_th > 0
        chi = min(max((N-beta0_th)/E_th,0),1);
    else
        chi = 0;
    end

    p_merge = 1-exp(-lambda*(2*dmax*v_rel*dt)*chi);
    p_break = min(max((2/pi)*(v_rel*dt/dmax)*chi,0),1-eps);
    p_disp = 1-(1-p_merge)*(1-p_break);
end
