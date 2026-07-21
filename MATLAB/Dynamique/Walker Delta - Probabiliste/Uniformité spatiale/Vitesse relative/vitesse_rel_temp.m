clear; clc; close all;

%% ============================================================
%  VITESSE RELATIVE TEMPORELLE — WALKER DELTA SPATIAL
%  COMPARAISON EMPIRIQUE / MODELE THEORIQUE PAR STRATES
%
%  Le modèle théorique n'utilise aucun angle ni lien empirique.
%
%  Pour chaque paire de phases orbitales (u1,u2), on intègre sur la
%  différence de RAAN DeltaOmega, supposée uniforme. On calcule :
%
%    G_m(u1,u2) =
%       P(lien et milieu du lien dans la strate m | u1,u2)
%
%    H_m(u1,u2) =
%       E[sin(gamma/2) 1_{lien,strate m} | u1,u2]
%
%  Puis, avec la densité de phase temporelle
%
%       f_u(u,t) = |cos(u-omega*t)|/4,
%
%  on obtient :
%
%       L_m^th(t) ∝ ∫∫ G_m(u1,u2) f_u(u1,t)f_u(u2,t) du1du2
%
%       s_m^th(t) =
%           ∫∫ H_m f_u f_u du1du2
%           --------------------------------
%           ∫∫ G_m f_u f_u du1du2
%
%  Enfin :
%
%       v_rel^th(t)
%       = 2 v_orb * sum_m L_m^th(t)s_m^th(t) / sum_m L_m^th(t)
%
%  C'est donc une prédiction théorique par quadrature, indépendante
%  des réalisations Monte-Carlo.
%% ============================================================

%% Paramètres physiques
R_earth = 6371;          % km
h = 550;                 % km
R = R_earth + h;         % km

mu = 398600;             % km^3/s^2
omega = sqrt(mu/R^3);    % rad/s
v_orb = sqrt(mu/R);      % km/s

%% Walker-Delta spatial
inc_deg = 90;
inc = deg2rad(inc_deg);

lambda = 4e-7;                       % sat/km^2 dans la bande
surface_band = 4*pi*R^2*sin(inc);
N_mean = lambda*surface_band;

%% Liens et temps
dmax = 1500;             % km
dt = 20;                 % s
Tmax = 12000;            % s

d_LOS = 2*sqrt(R^2-R_earth^2);
d_eff = min(dmax,d_LOS);

time_values = (0:dt:Tmax).';
Nt = numel(time_values);

%% Monte-Carlo empirique
n_iter = 30;

%% Strates de latitude absolue du milieu des liens
n_strates = 12;
alt_edges_deg = linspace(0,inc_deg,n_strates+1);
alt_centers_deg = 0.5*(alt_edges_deg(1:end-1)+alt_edges_deg(2:end));

%% Quadrature théorique
n_u = 48;                % noeuds pour chaque phase orbitale
n_dOmega = 120;          % points uniformes pour DeltaOmega

fprintf('Pré-calcul théorique : n_u = %d, n_dOmega = %d\n', ...
    n_u,n_dOmega);

%% ============================================================
%  PRE-CALCUL THEORIQUE DES NOYAUX PAR STRATE
%% ============================================================

[u_nodes,w_u] = gauss_legendre(n_u,0,2*pi);
dOmega_nodes = linspace(0,2*pi,n_dOmega+1);
dOmega_nodes(end) = [];
w_dOmega = 1/n_dOmega;

% G(:,:,m) : probabilité de lien dans la strate m
% H(:,:,m) : moyenne de sin(gamma/2)*1_lien,strate sur DeltaOmega
G = zeros(n_u,n_u,n_strates);
H = zeros(n_u,n_u,n_strates);

for a = 1:n_u
    u1 = u_nodes(a);

    % Satellite 1 avec Omega1 = 0
    r1 = walker_delta_positions(R,inc,0,u1);
    v1 = walker_delta_velocities(R,inc,0,u1,omega);

    for b = 1:n_u
        u2 = u_nodes(b);

        % Satellite 2 pour toutes les différences de RAAN
        Omega2 = dOmega_nodes(:);
        u2_vec = u2*ones(n_dOmega,1);

        r2 = walker_delta_positions(R,inc,Omega2,u2_vec);
        v2 = walker_delta_velocities(R,inc,Omega2,u2_vec,omega);

        dpos = r2-r1;
        distances = sqrt(sum(dpos.^2,2));
        linked = distances<=d_eff;

        if ~any(linked)
            continue;
        end

        v1_rep = repmat(v1,n_dOmega,1);
        cos_gamma = sum(v1_rep.*v2,2)/(v_orb^2);
        cos_gamma = max(-1,min(1,cos_gamma));
        sin_half = sqrt(max(0,(1-cos_gamma)/2));

        rmid = 0.5*(r2+r1);
        rmid_norm = sqrt(sum(rmid.^2,2));
        abs_lat_mid_deg = rad2deg(abs(asin(max(-1,min(1, ...
            rmid(:,3)./rmid_norm)))));

        bin_idx = discretize(abs_lat_mid_deg,alt_edges_deg);

        for m = 1:n_strates
            mask = linked & (bin_idx==m);

            if any(mask)
                G(a,b,m) = sum(mask)*w_dOmega;
                H(a,b,m) = sum(sin_half(mask))*w_dOmega;
            end
        end
    end
end

%% ============================================================
%  COURBES THEORIQUES TEMPORELLES
%% ============================================================

L_strate_theory = zeros(Nt,n_strates);
s_strate_theory = NaN(Nt,n_strates);
E_sin_half_theory = NaN(Nt,1);
gamma_eff_theory = NaN(Nt,1);
vrel_theory = NaN(Nt,1);

pair_factor = 0.5*N_mean^2;  % Poisson : E[N(N-1)]/2 = N_mean^2/2

for k = 1:Nt
    t = time_values(k);

    f = abs(cos(u_nodes-omega*t))/4;
    q = w_u.*f;

    % Renormalisation de la quadrature
    q = q/sum(q);
    W = q*q.';

    total_G = 0;
    total_H = 0;

    for m = 1:n_strates
        Gm = G(:,:,m);
        Hm = H(:,:,m);

        Pg_m = sum(W.*Gm,'all');
        Ph_m = sum(W.*Hm,'all');

        L_strate_theory(k,m) = pair_factor*Pg_m;

        if Pg_m>0
            s_strate_theory(k,m) = Ph_m/Pg_m;
        end

        total_G = total_G+Pg_m;
        total_H = total_H+Ph_m;
    end

    if total_G>0
        E_sin_half_theory(k) = total_H/total_G;
        gamma_eff_theory(k) = 2*asin(E_sin_half_theory(k));
        vrel_theory(k) = 2*v_orb*E_sin_half_theory(k);
    end
end

%% ============================================================
%  SIMULATION EMPIRIQUE MONTE-CARLO
%% ============================================================

vrel_direct_all = NaN(n_iter,Nt);
E_sin_half_emp_all = NaN(n_iter,Nt);
gamma_eff_emp_all = NaN(n_iter,Nt);
nb_links_all = zeros(n_iter,Nt);

L_strate_emp_all = zeros(n_iter,Nt,n_strates);

for it = 1:n_iter
    N = poissrnd(lambda*surface_band);

    fprintf('Simulation %d/%d : N = %d\n',it,n_iter,N);

    longitude0 = 2*pi*rand(N,1);
    sin_latitude0 = sin(inc)*(2*rand(N,1)-1);

    u_principal = asin(sin_latitude0/sin(inc));
    ascending = rand(N,1)<0.5;

    u0 = u_principal;
    u0(~ascending) = pi-u_principal(~ascending);
    u0 = mod(u0,2*pi);

    argument_longitude = atan2(sin(u0)*cos(inc),cos(u0));
    Omega = mod(longitude0-argument_longitude,2*pi);

    [pair_i,pair_j] = find(triu(true(N),1));

    for k = 1:Nt
        u_t = u0+omega*time_values(k);

        positions = walker_delta_positions(R,inc,Omega,u_t);
        velocities = walker_delta_velocities(R,inc,Omega,u_t,omega);

        dpos = positions(pair_i,:)-positions(pair_j,:);
        distances = sqrt(sum(dpos.^2,2));

        linked = distances<=d_eff;
        nb_links_all(it,k) = nnz(linked);

        if ~any(linked)
            continue;
        end

        dvel = velocities(pair_i(linked),:)-velocities(pair_j(linked),:);
        vrel_links = sqrt(sum(dvel.^2,2));

        vi = velocities(pair_i(linked),:);
        vj = velocities(pair_j(linked),:);

        cos_gamma = sum(vi.*vj,2)/(v_orb^2);
        cos_gamma = max(-1,min(1,cos_gamma));
        sin_half = sqrt(max(0,(1-cos_gamma)/2));

        vrel_direct_all(it,k) = mean(vrel_links);
        E_sin_half_emp_all(it,k) = mean(sin_half);
        gamma_eff_emp_all(it,k) = 2*asin(E_sin_half_emp_all(it,k));

        ri = positions(pair_i(linked),:);
        rj = positions(pair_j(linked),:);
        rmid = 0.5*(ri+rj);
        rmid_norm = sqrt(sum(rmid.^2,2));

        abs_lat_mid_deg = rad2deg(abs(asin(max(-1,min(1, ...
            rmid(:,3)./rmid_norm)))));

        bins = discretize(abs_lat_mid_deg,alt_edges_deg);

        for m = 1:n_strates
            L_strate_emp_all(it,k,m) = nnz(bins==m);
        end
    end
end

vrel_emp = mean(vrel_direct_all,1,'omitnan').';
E_sin_half_emp = mean(E_sin_half_emp_all,1,'omitnan').';
gamma_eff_emp = mean(gamma_eff_emp_all,1,'omitnan').';
nb_links_emp = mean(nb_links_all,1).';

L_strate_emp = reshape(mean(L_strate_emp_all,1),Nt,n_strates);

%% ============================================================
%  STATISTIQUES
%% ============================================================

rmse_vrel = sqrt(mean((vrel_emp-vrel_theory).^2,'omitnan'));
mae_vrel = mean(abs(vrel_emp-vrel_theory),'omitnan');

fprintf('\n--- Comparaison empirique / théorie ---\n');
fprintf('Vitesse empirique moyenne : %.4f km/s\n', ...
    mean(vrel_emp,'omitnan'));
fprintf('Vitesse théorique moyenne : %.4f km/s\n', ...
    mean(vrel_theory,'omitnan'));
fprintf('RMSE : %.4f km/s\n',rmse_vrel);
fprintf('MAE  : %.4f km/s\n',mae_vrel);

%% ============================================================
%  GRAPHES
%% ============================================================

figure;
hold on;
grid on;

plot(time_values,vrel_emp,'k','LineWidth',2.0, ...
    'DisplayName','Empirique Monte-Carlo');

plot(time_values,vrel_theory,'--','LineWidth',2.0, ...
    'DisplayName','Théorie par quadrature et strates');

xlabel('Temps (s)');
ylabel('Vitesse relative moyenne conditionnée aux liens (km/s)');
title(sprintf(['Walker-Delta spatial : vitesse relative empirique ', ...
    'et théorique — %d itérations'],n_iter));
legend('Location','best');
hold off;

figure;
hold on;
grid on;

plot(time_values,rad2deg(gamma_eff_emp),'k','LineWidth',1.8, ...
    'DisplayName','Empirique');

plot(time_values,rad2deg(gamma_eff_theory),'--','LineWidth',1.8, ...
    'DisplayName','Théorie');

xlabel('Temps (s)');
ylabel('\gamma_{eff}(t) (degrés)');
title('Angle effectif conditionné aux liens');
legend('Location','best');
hold off;

% Comparaison des liens par strate au pic théorique
[~,k_peak] = max(sum(L_strate_theory,2));

figure;
bar(alt_centers_deg, ...
    [L_strate_emp(k_peak,:).',L_strate_theory(k_peak,:).']);
grid on;
xlabel('Latitude absolue du milieu du lien (degrés)');
ylabel('Nombre moyen de liens');
title(sprintf('Liens par strate à t = %.0f s',time_values(k_peak)));
legend('Empirique','Théorie','Location','best');

% Facteur théorique s_m au pic
figure;
plot(alt_centers_deg,s_strate_theory(k_peak,:),'-o','LineWidth',1.6);
grid on;
xlabel('Latitude absolue du milieu du lien (degrés)');
ylabel('s_m^{th}=E[sin(\gamma/2)|lien,strate m]');
title(sprintf('Facteur angulaire théorique par strate à t = %.0f s', ...
    time_values(k_peak)));

%% ============================================================
%  SAUVEGARDE
%% ============================================================

save('vitesse_rel_temp.mat', ...
    'R_earth','h','R','mu','omega','v_orb', ...
    'inc_deg','inc','lambda','surface_band','N_mean', ...
    'dmax','d_eff','dt','Tmax','time_values', ...
    'n_iter','n_strates','alt_edges_deg','alt_centers_deg', ...
    'n_u','n_dOmega','G','H', ...
    'L_strate_theory','s_strate_theory', ...
    'E_sin_half_theory','gamma_eff_theory','vrel_theory', ...
    'L_strate_emp','E_sin_half_emp','gamma_eff_emp', ...
    'vrel_emp','nb_links_emp','rmse_vrel','mae_vrel');

fprintf('\nRésultats sauvegardés dans ');
fprintf('vitesse_rel_temp.mat\n');

%% ============================================================
%  FONCTIONS LOCALES
%% ============================================================

function positions = walker_delta_positions(R,inc,Omega,u)
    x = R*(cos(Omega).*cos(u)-sin(Omega).*sin(u).*cos(inc));
    y = R*(sin(Omega).*cos(u)+cos(Omega).*sin(u).*cos(inc));
    z = R*(sin(u).*sin(inc));
    positions = [x y z];
end

function velocities = walker_delta_velocities(R,inc,Omega,u,omega)
    vx = R*omega*(-cos(Omega).*sin(u)-sin(Omega).*cos(u).*cos(inc));
    vy = R*omega*(-sin(Omega).*sin(u)+cos(Omega).*cos(u).*cos(inc));
    vz = R*omega*cos(u).*sin(inc);
    velocities = [vx vy vz];
end

function [x,w] = gauss_legendre(n,a,b)
    k = 1:n-1;
    beta = k./sqrt(4*k.^2-1);
    J = diag(beta,1)+diag(beta,-1);
    [V,D] = eig(J);
    x = diag(D);
    [x,idx] = sort(x);
    V = V(:,idx);
    w = 2*(V(1,:).^2).';
    x = (b-a)/2*x+(a+b)/2;
    w = (b-a)/2*w;
end
