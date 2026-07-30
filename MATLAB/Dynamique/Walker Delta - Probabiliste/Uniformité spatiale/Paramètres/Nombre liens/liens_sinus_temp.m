clear; clc; close all;

%% ============================================================
%  WALKER-DELTA A UNIFORMITE SPATIALE
%  COMPARAISON EMPIRIQUE / SINUSOIDE MIN-MAX ANALYTIQUE
%
%  Minimum analytique :
%
%      L_min ~= E[C(N,2)] * (1-cos(alpha_max))/(2*sin(i))
%
%  Maximum analytique approché par strates :
%
%      mu_k^max = N/sin(i) * [
%          sqrt(sin(i)^2-sin(phi_k)^2)
%        - sqrt(sin(i)^2-sin(phi_{k+1})^2)
%      ]
%
%      L_max ~= sum_k 1/2 * mu_k^max * k_k^max
%
%  Sinusoïde :
%
%      L_sin(t) = (L_max+L_min)/2
%               - (L_max-L_min)/2*cos(2*omega*t)
%
%  La période est T_orb/2.
%% ============================================================

%% Paramètres physiques
R_earth = 6371;          % km
h = 550;                 % km
R = R_earth + h;         % km

mu_earth = 398600;       % km^3/s^2
omega = sqrt(mu_earth/R^3);

%% Paramètres Walker-Delta
inc_deg = 58;
inc = deg2rad(inc_deg);

lambda = 4e-7;                       % sat/km^2 dans la bande
surface_band = 4*pi*R^2*sin(inc);
N_mean = lambda*surface_band;

%% Paramètres de liaison
dmax = 1500;             % km

d_LOS = 2*sqrt(R^2-R_earth^2);
d_eff = min(dmax,d_LOS);

alpha_max = 2*asin(d_eff/(2*R));
A_link = 2*pi*R^2*(1-cos(alpha_max));

%% Paramètres temporels
dt = 60;
Tmax = 12000;
time_values = 0:dt:Tmax;
Nt = numel(time_values);

T_orb = 2*pi/omega;
T_links = T_orb/2;

%% Paramètres Monte-Carlo
n_iter = 100;

%% Paramètres des strates
phi_step = alpha_max/16;

phi_edges = 0:phi_step:inc;

if phi_edges(end) < inc-1e-12
    phi_edges(end+1) = inc;
else
    phi_edges(end) = inc;
end

phi_edges = unique(phi_edges,'stable');

phi_in = phi_edges(1:end-1).';
phi_out = phi_edges(2:end).';
phi_mid = 0.5*(phi_in+phi_out);

n_strates = numel(phi_in);

%% ============================================================
%  MINIMUM THEORIQUE
%% ============================================================

p_link_min = min(1,(1-cos(alpha_max))/(2*sin(inc)));

% Pour N ~ Poisson(N_mean), E[N(N-1)] = N_mean^2.
L_min_theory = 0.5*N_mean^2*p_link_min;

%% ============================================================
%  MAXIMUM THEORIQUE PAR STRATES
%% ============================================================

sin_i = sin(inc);

term_in = sqrt(max(0,sin_i^2-sin(phi_in).^2));
term_out = sqrt(max(0,sin_i^2-sin(phi_out).^2));

mu_max = N_mean/sin_i .* (term_in-term_out);

A_strate = 4*pi*R^2 .* (sin(phi_out)-sin(phi_in));
lambda_local_max = mu_max ./ max(A_strate,eps);

%% Correction simple de bord
distance_to_edge = inc-phi_mid;

eta_edge = min(1,0.5+distance_to_edge/(2*alpha_max));
eta_edge = max(0.5,eta_edge);

A_link_eff = eta_edge.*A_link;

k_mean_max = lambda_local_max.*A_link_eff;
k_mean_max = min(k_mean_max,max(mu_max-1,0));

L_strate_max = 0.5.*mu_max.*k_mean_max;
L_max_theory = sum(L_strate_max);

%% ============================================================
%  SINUSOIDE THEORIQUE
%% ============================================================

L_center = 0.5*(L_max_theory+L_min_theory);
L_amplitude = 0.5*(L_max_theory-L_min_theory);

L_sin_theory = L_center ...
    - L_amplitude*cos(2*pi*time_values/T_links);

%% ============================================================
%  SIMULATION MONTE-CARLO
%% ============================================================

num_edges_all = zeros(n_iter,Nt);
N_all = zeros(n_iter,1);

links_strate_internal_all = zeros(n_iter,Nt,n_strates);
n_sat_strate_all = zeros(n_iter,Nt,n_strates);

for it = 1:n_iter

    N = poissrnd(lambda*surface_band);
    N_all(it) = N;

    fprintf('Itération %d / %d : N = %d\n',it,n_iter,N);

    %% Uniformité spatiale initiale dans la bande
    longitude0 = 2*pi*rand(N,1);
    sin_latitude0 = sin(inc)*(2*rand(N,1)-1);

    %% Reconstruction de u0 et Omega
    u_principal = asin(sin_latitude0/sin(inc));
    ascending_branch = rand(N,1)<0.5;

    u0 = u_principal;
    u0(~ascending_branch) = pi-u_principal(~ascending_branch);
    u0 = mod(u0,2*pi);

    argument_longitude = atan2(sin(u0)*cos(inc),cos(u0));
    Omega = mod(longitude0-argument_longitude,2*pi);

    for k = 1:Nt

        t = time_values(k);
        u_t = u0+omega*t;

        positions_t = walker_delta_positions(R,inc,Omega,u_t);

        D = squareform(pdist(positions_t));
        A = (D<=d_eff) & (D>0);

        num_edges_all(it,k) = nnz(triu(A,1));

        [L_internal,n_sat] = ...
            liens_empiriques_par_strate(A,positions_t,R,phi_in,phi_out);

        links_strate_internal_all(it,k,:) = reshape(L_internal,1,1,[]);
        n_sat_strate_all(it,k,:) = reshape(n_sat,1,1,[]);
    end
end

%% ============================================================
%  STATISTIQUES EMPIRIQUES
%% ============================================================

mean_edges = mean(num_edges_all,1);
mean_edges_global = mean(mean_edges);

mean_links_internal_strate_t = ...
    reshape(mean(links_strate_internal_all,1),Nt,n_strates);

mean_n_sat_strate_t = ...
    reshape(mean(n_sat_strate_all,1),Nt,n_strates);

[mean_edges_peak,idx_peak_edges] = max(mean_edges);
t_peak_edges = time_values(idx_peak_edges);

[mean_edges_min,idx_min_edges] = min(mean_edges);
t_min_edges = time_values(idx_min_edges);

L_internal_emp_peak = mean_links_internal_strate_t(idx_peak_edges,:).';
N_emp_peak = mean_n_sat_strate_t(idx_peak_edges,:).';

mean_N_emp = mean(N_all);

%% ============================================================
%  ERREURS
%% ============================================================

rmse_sin = sqrt(mean((mean_edges-L_sin_theory).^2));
mae_sin = mean(abs(mean_edges-L_sin_theory));

relative_error_min = abs(L_min_theory-mean_edges_min)/max(mean_edges_min,eps);
relative_error_max = abs(L_max_theory-mean_edges_peak)/max(mean_edges_peak,eps);

%% ============================================================
%  TABLE DE COMPARAISON PAR STRATE
%% ============================================================

strates_comparison = table( ...
    (1:n_strates).', ...
    phi_in,phi_out,phi_mid, ...
    rad2deg(phi_in),rad2deg(phi_out),rad2deg(phi_mid), ...
    A_strate,mu_max,N_emp_peak, ...
    lambda_local_max,eta_edge,A_link_eff, ...
    k_mean_max,L_strate_max,L_internal_emp_peak, ...
    'VariableNames',{ ...
    'index', ...
    'phi_in','phi_out','phi_mid', ...
    'phi_in_deg','phi_out_deg','phi_mid_deg', ...
    'A_strate','N_theory_max','N_emp_peak', ...
    'lambda_local_max','eta_edge','A_link_eff', ...
    'k_mean_max','L_internal_theory_max','L_internal_emp_peak'});

strates_comparison.ratio_emp_theory = ...
    strates_comparison.L_internal_emp_peak ./ ...
    max(strates_comparison.L_internal_theory_max,eps);

%% ============================================================
%  AFFICHAGE TEXTE
%% ============================================================

fprintf('\n=== Walker-Delta spatial : comparaison empirique / analytique ===\n');
fprintf('Inclinaison : %.2f deg\n',inc_deg);
fprintf('N théorique moyen : %.3f satellites\n',N_mean);
fprintf('N empirique moyen : %.3f satellites\n',mean_N_emp);
fprintf('alpha_max : %.5f rad = %.3f deg\n', ...
    alpha_max,rad2deg(alpha_max));
fprintf('Période orbitale : %.3f s\n',T_orb);
fprintf('Période des liens : %.3f s\n',T_links);

fprintf('\nMinimum théorique : %.3f liens\n',L_min_theory);
fprintf('Minimum empirique moyen : %.3f liens à t = %.1f s\n', ...
    mean_edges_min,t_min_edges);
fprintf('Erreur relative minimum : %.2f %%\n',100*relative_error_min);

fprintf('\nMaximum théorique : %.3f liens\n',L_max_theory);
fprintf('Maximum empirique moyen : %.3f liens à t = %.1f s\n', ...
    mean_edges_peak,t_peak_edges);
fprintf('Erreur relative maximum : %.2f %%\n',100*relative_error_max);

fprintf('\nCentre théorique : %.3f liens\n',L_center);
fprintf('Amplitude théorique : %.3f liens\n',L_amplitude);
fprintf('RMSE sinusoïde : %.3f liens\n',rmse_sin);
fprintf('MAE sinusoïde : %.3f liens\n\n',mae_sin);

disp(strates_comparison);

%% ============================================================
%  GRAPHE PRINCIPAL
%% ============================================================

figure;
hold on;
grid on;

plot(time_values,mean_edges,'k','LineWidth',2.2, ...
    'DisplayName','Moyenne empirique');

plot(time_values,L_sin_theory,'--','LineWidth',2.2, ...
    'DisplayName','Sinusoïde analytique min/max');

yline(L_min_theory,':','LineWidth',1.4, ...
    'DisplayName',sprintf('Minimum th. = %.1f',L_min_theory));

yline(L_max_theory,':','LineWidth',1.4, ...
    'DisplayName',sprintf('Maximum th. = %.1f',L_max_theory));

xlabel('Temps (s)');
ylabel('Nombre de liens inter-satellites');

title(sprintf(['Walker-Delta spatial : moyenne empirique et ', ...
    'sinusoïde analytique — %d itérations'],n_iter));

legend('Location','best');
hold off;

%% ============================================================
%  COMPARAISON DES SATELLITES PAR STRATE AU PIC
%% ============================================================

figure;
bar(strates_comparison.index, ...
    [strates_comparison.N_emp_peak, ...
     strates_comparison.N_theory_max]);

grid on;
xlabel('Indice de strate');
ylabel('Nombre moyen de satellites');

title(sprintf('Satellites par strate au pic empirique, t = %.0f s', ...
    t_peak_edges));

legend('Empirique','Théorie analytique','Location','best');

%% ============================================================
%  COMPARAISON DES LIENS INTERNES PAR STRATE
%% ============================================================

figure;
bar(strates_comparison.index, ...
    [strates_comparison.L_internal_emp_peak, ...
     strates_comparison.L_internal_theory_max]);

grid on;
xlabel('Indice de strate');
ylabel('Nombre de liens internes');

title(sprintf('Liens internes par strate au pic, t = %.0f s', ...
    t_peak_edges));

legend('Empirique','Théorie analytique','Location','best');

%% ============================================================
%  SAUVEGARDE
%% ============================================================

save('liens_sinus_results.mat', ...
    'R_earth','h','R','mu_earth','omega', ...
    'inc_deg','inc','lambda','surface_band','N_mean', ...
    'dmax','d_LOS','d_eff','alpha_max','A_link', ...
    'dt','Tmax','time_values','T_orb','T_links','n_iter', ...
    'phi_step','phi_edges','strates_comparison', ...
    'p_link_min','L_min_theory','L_max_theory', ...
    'L_center','L_amplitude','L_sin_theory', ...
    'N_all','num_edges_all','mean_edges','mean_edges_global', ...
    'mean_edges_peak','t_peak_edges','mean_edges_min','t_min_edges', ...
    'mean_links_internal_strate_t','mean_n_sat_strate_t', ...
    'rmse_sin','mae_sin','relative_error_min','relative_error_max');

fprintf('\nRésultats sauvegardés dans ');
fprintf('liens_inter_satellites_sinus.mat\n');

%% ============================================================
%  FONCTIONS LOCALES
%% ============================================================

function positions = walker_delta_positions(R,inc,Omega,u)

    x = R*(cos(Omega).*cos(u) ...
        - sin(Omega).*sin(u).*cos(inc));

    y = R*(sin(Omega).*cos(u) ...
        + cos(Omega).*sin(u).*cos(inc));

    z = R*(sin(u).*sin(inc));

    positions = [x y z];
end

function [L_internal,n_sat] = ...
    liens_empiriques_par_strate(A,positions,R,phi_in,phi_out)

    n_strates = numel(phi_in);

    L_internal = zeros(n_strates,1);
    n_sat = zeros(n_strates,1);

    latitude = asin(max(-1,min(1,positions(:,3)/R)));
    abs_latitude = abs(latitude);

    [row,col] = find(triu(A,1));

    for s = 1:n_strates

        if s<n_strates
            in_s = abs_latitude>=phi_in(s) & ...
                abs_latitude<phi_out(s);
        else
            in_s = abs_latitude>=phi_in(s) & ...
                abs_latitude<=phi_out(s)+1e-12;
        end

        n_sat(s) = nnz(in_s);

        if isempty(row)
            continue;
        end

        L_internal(s) = nnz(in_s(row) & in_s(col));
    end
end
