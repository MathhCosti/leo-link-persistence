clear; clc; close all;

%% ============================================================
%  NOMBRE DE LIENS INTER-SATELLITES TEMPOREL
%  WALKER-DELTA A UNIFORMITE SPATIALE INITIALE
%
%  Correction principale :
%  - le maximum théorique n'est plus obtenu en transposant artificiellement
%    le modèle de transport par strates du Walker-Star ;
%  - la probabilité de lien temporelle est calculée directement à partir
%    de la loi exacte de phase induite par l'uniformité spatiale initiale :
%
%        f_u(u,t) = |cos(u-omega*t)|/4.
%
%  L'intégration sur la différence de RAAN est effectuée analytiquement,
%  puis l'intégrale double en phase est évaluée par Gauss-Legendre.
%
%  Les comparaisons par strate mettent toujours en regard des quantités
%  de même nature :
%  - liens internes empiriques / liens internes théoriques ;
%  - liens incidents empiriques / liens incidents théoriques.
%% ============================================================

%% Paramètres physiques
R_earth = 6371;          % km
h = 550;                 % km
R = R_earth + h;         % rayon orbital [km]

mu = 398600;             % km^3/s^2
omega = sqrt(mu/R^3);    % rad/s

%% Walker-Delta et processus de Poisson
inc_deg = 58;
inc = deg2rad(inc_deg);

lambda = 4e-7;                       % sat/km^2 dans la bande
surface_band = 4*pi*R^2*sin(inc);
N_mean_theory = lambda*surface_band;

%% Liens et temps
dmax = 1500;             % km
dt = 60;                 % s
Tmax = 12000;            % s

time_values = 0:dt:Tmax;
Nt = numel(time_values);

%% Monte-Carlo réseau
n_iter = 100;

%% Quadrature théorique
n_quad_u = 120;          % augmenter pour plus de précision

%% ============================================================
%  GRANDEURS GEOMETRIQUES
%% ============================================================

d_LOS = 2*sqrt(R^2-R_earth^2);
d_eff = min(dmax,d_LOS);

alpha_max = 2*asin(d_eff/(2*R));
cmax = cos(alpha_max);

T_orb = 2*pi/omega;
T_links = T_orb/2;

fprintf('N moyen théorique dans la bande : %.3f\n',N_mean_theory);
fprintf('Inclinaison : %.2f deg\n',inc_deg);
fprintf('alpha_max : %.5f rad = %.3f deg\n', ...
    alpha_max,rad2deg(alpha_max));
fprintf('Période orbitale : %.2f s\n',T_orb);
fprintf('Demi-période orbitale : %.2f s\n',T_links);

%% ============================================================
%  STRATES DE LATITUDE DU DELTA
%% ============================================================

beta_step_strates = alpha_max;
beta_max_strates = inc;

strates_delta = strates_delta_spatiales( ...
    lambda,R,alpha_max,inc, ...
    'beta_step',beta_step_strates, ...
    'beta_max',beta_max_strates, ...
    'verbose',true);

n_strates = height(strates_delta.active_table);

%% ============================================================
%  MODELE THEORIQUE TEMPOREL EXACT SOUS FORME D'INTEGRALE
%% ============================================================

% Noeuds et poids sur [0,2*pi].
[u_nodes,w_nodes] = gauss_legendre(n_quad_u,0,2*pi);

% g(u1,u2) est la fraction des différences de RAAN donnant un lien.
G_link = fraction_raan_liee(u_nodes,u_nodes,inc,cmax);

p_link_theory_time = zeros(Nt,1);
L_theory_time = zeros(Nt,1);

% Nombre théorique moyen de satellites dans chaque strate.
N_strate_theory_time = zeros(Nt,n_strates);

for k = 1:Nt
    t = time_values(k);

    % Loi de la phase courante issue de l'uniformité spatiale à t=0.
    f_u = abs(cos(u_nodes-omega*t))/4;
    wf = w_nodes.*f_u;

    % Sécurité numérique : renormalisation de la quadrature.
    wf = wf/sum(wf);

    p_link_theory_time(k) = wf.'*G_link*wf;
    L_theory_time(k) = N_mean_theory*(N_mean_theory-1)/2 ...
        *p_link_theory_time(k);

    abs_lat = abs(asin(sin(inc).*sin(u_nodes)));

    for s = 1:n_strates
        lat_min = strates_delta.active_table.latitude_inner(s);
        lat_max = strates_delta.active_table.latitude_outer(s);

        if s < n_strates
            mask = abs_lat >= lat_min & abs_lat < lat_max;
        else
            mask = abs_lat >= lat_min & abs_lat <= lat_max+1e-12;
        end

        N_strate_theory_time(k,s) = N_mean_theory*sum(wf(mask));
    end
end

[L_min_theory,idx_min_theory] = min(L_theory_time);
[L_max_theory,idx_max_theory] = max(L_theory_time);
t_min_theory = time_values(idx_min_theory);
t_max_theory = time_values(idx_max_theory);

fprintf('\n--- Modèle théorique temporel par quadrature ---\n');
fprintf('p_link théorique min : %.8f\n',min(p_link_theory_time));
fprintf('p_link théorique max : %.8f\n',max(p_link_theory_time));
fprintf('L théorique min : %.3f liens à t = %.1f s\n', ...
    L_min_theory,t_min_theory);
fprintf('L théorique max : %.3f liens à t = %.1f s\n', ...
    L_max_theory,t_max_theory);

%% ============================================================
%  STOCKAGE MONTE-CARLO
%% ============================================================

num_edges_all = zeros(n_iter,Nt);
N_all = zeros(n_iter,1);

links_strate_internal_all = zeros(n_iter,Nt,n_strates);
links_strate_incident_all = zeros(n_iter,Nt,n_strates);
n_sat_strate_all = zeros(n_iter,Nt,n_strates);

%% ============================================================
%  BOUCLE MONTE-CARLO
%% ============================================================

for it = 1:n_iter

    N = poissrnd(lambda*surface_band);
    N_all(it) = N;

    fprintf('Itération %d/%d : N = %d\n',it,n_iter,N);

    %% Uniformité spatiale initiale dans |latitude| <= inc
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

        [L_internal,L_incident,n_sat] = ...
            liens_empiriques_delta(A,positions_t,R,strates_delta);

        links_strate_internal_all(it,k,:) = reshape(L_internal,1,1,[]);
        links_strate_incident_all(it,k,:) = reshape(L_incident,1,1,[]);
        n_sat_strate_all(it,k,:) = reshape(n_sat,1,1,[]);
    end
end

%% ============================================================
%  STATISTIQUES EMPIRIQUES
%% ============================================================

mean_edges = mean(num_edges_all,1).';
mean_edges_global = mean(mean_edges);

mean_links_internal_strate_t = ...
    reshape(mean(links_strate_internal_all,1),Nt,n_strates);
mean_links_incident_strate_t = ...
    reshape(mean(links_strate_incident_all,1),Nt,n_strates);
mean_n_sat_strate_t = ...
    reshape(mean(n_sat_strate_all,1),Nt,n_strates);

[mean_edges_peak,idx_peak_edges] = max(mean_edges);
t_peak_edges = time_values(idx_peak_edges);

mean_N = mean(N_all);

fprintf('\n--- Résultats Monte-Carlo ---\n');
fprintf('N empirique moyen : %.3f\n',mean_N);
fprintf('Nombre moyen global de liens : %.3f\n',mean_edges_global);
fprintf('Pic empirique moyen : %.3f liens à t = %.1f s\n', ...
    mean_edges_peak,t_peak_edges);

%% ============================================================
%  THEORIE PAR STRATE AU PIC EMPIRIQUE
%% ============================================================

t_ref = t_peak_edges;
f_ref = abs(cos(u_nodes-omega*t_ref))/4;
wf_ref = w_nodes.*f_ref;
wf_ref = wf_ref/sum(wf_ref);

abs_lat_nodes = abs(asin(sin(inc).*sin(u_nodes)));
masks = false(n_quad_u,n_strates);

for s = 1:n_strates
    lat_min = strates_delta.active_table.latitude_inner(s);
    lat_max = strates_delta.active_table.latitude_outer(s);

    if s < n_strates
        masks(:,s) = abs_lat_nodes>=lat_min & abs_lat_nodes<lat_max;
    else
        masks(:,s) = abs_lat_nodes>=lat_min & ...
            abs_lat_nodes<=lat_max+1e-12;
    end
end

L_internal_theory_peak = zeros(n_strates,1);
L_incident_theory_peak = zeros(n_strates,1);

pair_factor = N_mean_theory*(N_mean_theory-1)/2;

for s = 1:n_strates
    m = double(masks(:,s));

    % Deux extrémités dans la strate.
    W_internal = (wf_ref.*m)*(wf_ref.*m).';

    % Au moins une extrémité dans la strate.
    M_incident = (m+m.')>0;
    W_all = wf_ref*wf_ref.';

    p_internal_s = sum(G_link.*W_internal,'all');
    p_incident_s = sum(G_link.*W_all.*M_incident,'all');

    L_internal_theory_peak(s) = pair_factor*p_internal_s;
    L_incident_theory_peak(s) = pair_factor*p_incident_s;
end

L_internal_emp_peak = mean_links_internal_strate_t(idx_peak_edges,:).';
L_incident_emp_peak = mean_links_incident_strate_t(idx_peak_edges,:).';
N_emp_peak = mean_n_sat_strate_t(idx_peak_edges,:).';
N_theory_peak = N_strate_theory_time(idx_peak_edges,:).';

comparison_strates = strates_delta.active_table(:, ...
    {'index','latitude_inner_deg','latitude_outer_deg'});

comparison_strates.N_emp_peak = N_emp_peak;
comparison_strates.N_theory_peak = N_theory_peak;
comparison_strates.L_internal_emp_peak = L_internal_emp_peak;
comparison_strates.L_internal_theory_peak = L_internal_theory_peak;
comparison_strates.L_incident_emp_peak = L_incident_emp_peak;
comparison_strates.L_incident_theory_peak = L_incident_theory_peak;

disp(comparison_strates);

%% ============================================================
%  QUALITE DU MODELE
%% ============================================================

rmse_links = sqrt(mean((mean_edges-L_theory_time).^2));
mae_links = mean(abs(mean_edges-L_theory_time));

fprintf('\n--- Écart théorie / simulation ---\n');
fprintf('RMSE nombre de liens : %.4f\n',rmse_links);
fprintf('MAE nombre de liens  : %.4f\n',mae_links);

%% ============================================================
%  GRAPHE TEMPOREL PRINCIPAL
%% ============================================================

figure;
hold on;
grid on;

plot(time_values,mean_edges,'k','LineWidth',2.2, ...
    'DisplayName','Moyenne empirique');
plot(time_values,L_theory_time,'--','LineWidth',2.0, ...
    'DisplayName','Théorie par quadrature');

yline(L_min_theory,':','LineWidth',1.4, ...
    'DisplayName',sprintf('Minimum théorique = %.1f',L_min_theory));
yline(L_max_theory,':','LineWidth',1.4, ...
    'DisplayName',sprintf('Maximum théorique = %.1f',L_max_theory));

xlabel('Temps (s)');
ylabel('Nombre de liens inter-satellites');
title(sprintf(['Walker-Delta spatial : liens empiriques et modèle ', ...
    'temporel — %d itérations'],n_iter));
legend('Location','best');
hold off;

%% Probabilité de lien temporelle
figure;
plot(time_values,p_link_theory_time,'LineWidth',2);
grid on;
xlabel('Temps (s)');
ylabel('p_{link}^{\Delta,sp}(t)');
title('Probabilité théorique temporelle de lien');

%% Satellites par strate
figure;
bar(comparison_strates.index, ...
    [comparison_strates.N_emp_peak,comparison_strates.N_theory_peak]);
grid on;
xlabel('Indice de strate');
ylabel('Nombre de satellites');
title(sprintf('Satellites par strate au pic empirique, t = %.0f s',t_ref));
legend('Empirique','Théorie','Location','best');

%% Liens internes : comparaison cohérente
figure;
bar(comparison_strates.index, ...
    [comparison_strates.L_internal_emp_peak, ...
     comparison_strates.L_internal_theory_peak]);
grid on;
xlabel('Indice de strate');
ylabel('Nombre de liens internes');
title(sprintf('Liens internes par strate, t = %.0f s',t_ref));
legend('Empirique','Théorie','Location','best');

%% Liens incidents : comparaison cohérente
figure;
bar(comparison_strates.index, ...
    [comparison_strates.L_incident_emp_peak, ...
     comparison_strates.L_incident_theory_peak]);
grid on;
xlabel('Indice de strate');
ylabel('Nombre de liens incidents');
title(sprintf('Liens incidents par strate, t = %.0f s',t_ref));
legend('Empirique','Théorie','Location','best');

%% ============================================================
%  SAUVEGARDE
%% ============================================================

save('liens_inter_satellites_quadrature.mat', ...
    'R_earth','h','R','mu','omega', ...
    'inc_deg','inc','lambda','surface_band','N_mean_theory', ...
    'dmax','d_LOS','d_eff','alpha_max','cmax', ...
    'dt','Tmax','time_values','T_orb','T_links','n_iter', ...
    'strates_delta','n_quad_u', ...
    'p_link_theory_time','L_theory_time', ...
    'L_min_theory','L_max_theory','t_min_theory','t_max_theory', ...
    'N_strate_theory_time', ...
    'N_all','num_edges_all','mean_edges','mean_edges_global', ...
    'links_strate_internal_all','links_strate_incident_all', ...
    'n_sat_strate_all','mean_links_internal_strate_t', ...
    'mean_links_incident_strate_t','mean_n_sat_strate_t', ...
    'idx_peak_edges','t_peak_edges','mean_edges_peak', ...
    'comparison_strates','rmse_links','mae_links');

fprintf('\nRésultats sauvegardés dans liens_inter_satellites_quadrature.mat\n');

%% ============================================================
%  FONCTIONS LOCALES
%% ============================================================

function positions = walker_delta_positions(R,inc,Omega,u)
    x = R*(cos(Omega).*cos(u)-sin(Omega).*sin(u).*cos(inc));
    y = R*(sin(Omega).*cos(u)+cos(Omega).*sin(u).*cos(inc));
    z = R*(sin(u).*sin(inc));
    positions = [x y z];
end

function G = fraction_raan_liee(u1,u2,inc,cmax)
% Fraction de DeltaOmega uniforme dans [0,2*pi) donnant un lien.

    c1 = cos(u1(:));
    s1 = sin(u1(:));
    c2 = cos(u2(:)).';
    s2 = sin(u2(:)).';

    A = c1*c2 + cos(inc)^2*(s1*s2);
    B = cos(inc)*(s1*c2-c1*s2);
    C = sin(inc)^2*(s1*s2);

    rho = sqrt(A.^2+B.^2);
    G = zeros(size(rho));

    regular = rho>1e-13;
    q = zeros(size(rho));
    q(regular) = (cmax-C(regular))./rho(regular);

    G(regular & q<=-1) = 1;

    middle = regular & q>-1 & q<1;
    G(middle) = acos(q(middle))/pi;

    degenerate = ~regular;
    G(degenerate) = double(C(degenerate)>=cmax);
end

function [L_internal,L_incident,n_sat] = ...
    liens_empiriques_delta(A,positions,R,strates)

    T = strates.active_table;
    ns = height(T);

    L_internal = zeros(ns,1);
    L_incident = zeros(ns,1);
    n_sat = zeros(ns,1);

    latitude = asin(max(-1,min(1,positions(:,3)/R)));
    abs_latitude = abs(latitude);

    [row,col] = find(triu(A,1));

    for s = 1:ns
        lat_min = T.latitude_inner(s);
        lat_max = T.latitude_outer(s);

        if s<ns
            in_s = abs_latitude>=lat_min & abs_latitude<lat_max;
        else
            in_s = abs_latitude>=lat_min & ...
                abs_latitude<=lat_max+1e-12;
        end

        n_sat(s) = nnz(in_s);

        if isempty(row)
            continue;
        end

        row_in = in_s(row);
        col_in = in_s(col);

        L_internal(s) = nnz(row_in & col_in);
        L_incident(s) = nnz(row_in | col_in);
    end
end

function [x,w] = gauss_legendre(n,a,b)
% Quadrature de Gauss-Legendre par l'algorithme de Golub-Welsch.

    k = 1:n-1;
    beta = k./sqrt(4*k.^2-1);

    J = diag(beta,1)+diag(beta,-1);
    [V,D] = eig(J);

    x = diag(D);
    [x,order] = sort(x);
    V = V(:,order);

    w = 2*(V(1,:).^2).';

    x = (b-a)/2*x+(a+b)/2;
    w = (b-a)/2*w;
end
