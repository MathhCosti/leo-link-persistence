clear; clc; close all; rng(4);

%% ============================================================
%  ESTIMATION AMELIOREE DE N1, N2 ET N3
%  Walker Delta a uniformite orbitale
%
%  N1 : correction locale en latitude via p_link(phi)
%  N2 : esperance geometrique exacte, evaluee par Monte-Carlo imbrique
%  N3 : esperance geometrique exacte, evaluee par Monte-Carlo imbrique
%
%  Comparaison avec des graphes statiques simules.
%% ============================================================

%% Parametres physiques
R = 6371 + 550;      % km
inc_deg = 53;
inc = deg2rad(inc_deg);
dmax = 1500;         % km
N = 250;

%% Parametres numeriques
% Pour N1
n_phi = 180;
nQuad = 220;

% Pour N2
n_pair_configs = 2e4;   % configurations de paires
n_probe_pair = 2e4;     % satellites tests pour q_union,2

% Pour N3
n_triplet_configs = 2e4;
n_probe_triplet = 2e4;

% Validation par graphes complets
n_graph_realizations = 300;

%% Seuil de lien
alpha_max = 2*asin(min(dmax/(2*R),1));
cos_alpha_max = cos(alpha_max);

%% ============================================================
%  1. N1 : approximation locale en latitude
%% ============================================================

eps_phi = 1e-6;
phi_vals = linspace(-inc+eps_phi,inc-eps_phi,n_phi);
p_link_phi = zeros(size(phi_vals));

[u2,w2] = gauss_legendre_interval(nQuad,0,2*pi);
u2 = u2(:);
w2 = w2(:);

c2 = cos(u2);
s2 = sin(u2);
ci = cos(inc);
si = sin(inc);

for k = 1:numel(phi_vals)
    phi = phi_vals(k);

    s1 = sin(phi)/si;
    s1 = max(min(s1,1),-1);

    u1a = asin(s1);
    u1b = pi-u1a;

    branches = [u1a,u1b];
    p_branch = zeros(2,1);

    for b = 1:2
        u1 = branches(b);
        c1 = cos(u1);
        s1b = sin(u1);

        A = c1.*c2 + ci^2.*s1b.*s2;
        B = ci.*(s1b.*c2 - c1.*s2);
        C = si^2.*s1b.*s2;

        rho = sqrt(A.^2+B.^2);

        g = zeros(size(rho));
        mask = rho > 1e-14;

        q = zeros(size(rho));
        q(mask) = (cos_alpha_max-C(mask))./rho(mask);

        g(mask & q <= -1) = 1;

        middle = mask & q > -1 & q < 1;
        g(middle) = acos(q(middle))/pi;

        g(~mask) = double(C(~mask) >= cos_alpha_max);

        p_branch(b) = sum(w2.*g)/(2*pi);
    end

    p_link_phi(k) = mean(p_branch);
end

%% Calcul exact de N1 par integration sur la phase orbitale u1

nQuad_N1 = 300;
[u1_nodes, w1_nodes] = gauss_legendre_interval(nQuad_N1, 0, 2*pi);

p_link_u1 = zeros(nQuad_N1,1);

for k = 1:nQuad_N1

    u1 = u1_nodes(k);

    c1 = cos(u1);
    s1 = sin(u1);

    % Produit scalaire :
    % dot = A*cos(DeltaOmega) + B*sin(DeltaOmega) + C
    A = c1.*c2 + ci^2.*s1.*s2;
    B = ci.*(s1.*c2 - c1.*s2);
    C = si^2.*s1.*s2;

    rho = sqrt(A.^2 + B.^2);

    g = zeros(size(rho));
    mask = rho > 1e-14;

    q = zeros(size(rho));
    q(mask) = (cos_alpha_max - C(mask))./rho(mask);

    g(mask & q <= -1) = 1;

    middle = mask & q > -1 & q < 1;
    g(middle) = acos(q(middle))/pi;

    g(~mask) = double(C(~mask) >= cos_alpha_max);

    % Moyenne sur la phase u2
    p_link_u1(k) = sum(w2 .* g)/(2*pi);
end

% Phase u1 uniforme sur [0,2*pi)
N1_th_local = N/(2*pi) * ...
    sum(w1_nodes .* (1-p_link_u1).^(N-1));

%% ============================================================
%  2. N2 : dimeres isoles
%
%  N2 = C(N,2) E[ 1_{1~2} (1-q_union,2)^(N-2) ]
%
%  q_union,2(x1,x2) = P(X3 lie a x1 ou x2)
%% ============================================================

P1 = sample_walker_delta_unit(n_pair_configs,inc);
P2 = sample_walker_delta_unit(n_pair_configs,inc);

linked_pair = sum(P1.*P2,2) >= cos_alpha_max;

% Echantillon de satellites tests partage par toutes les configurations
Qpair = sample_walker_delta_unit(n_probe_pair,inc);

pair_terms = zeros(n_pair_configs,1);

linked_idx = find(linked_pair);

for ii = 1:numel(linked_idx)
    j = linked_idx(ii);

    link_to_1 = Qpair*P1(j,:).' >= cos_alpha_max;
    link_to_2 = Qpair*P2(j,:).' >= cos_alpha_max;

    q_union_2 = mean(link_to_1 | link_to_2);

    pair_terms(j) = (1-q_union_2)^(N-2);
end

N2_th_geom = nchoosek(N,2)*mean(pair_terms);

%% ============================================================
%  3. N3 : trimeres isoles
%
%  N3 = C(N,3) E[
%          1_{G3 connexe}
%          (1-q_union,3)^(N-3)
%       ]
%
%  q_union,3 = P(X4 lie a au moins un des trois satellites)
%% ============================================================

T1 = sample_walker_delta_unit(n_triplet_configs,inc);
T2 = sample_walker_delta_unit(n_triplet_configs,inc);
T3 = sample_walker_delta_unit(n_triplet_configs,inc);

e12 = sum(T1.*T2,2) >= cos_alpha_max;
e13 = sum(T1.*T3,2) >= cos_alpha_max;
e23 = sum(T2.*T3,2) >= cos_alpha_max;

% Un graphe a 3 sommets est connexe s'il possede au moins 2 aretes
connected_triplet = (e12 + e13 + e23) >= 2;

Qtrip = sample_walker_delta_unit(n_probe_triplet,inc);

triplet_terms = zeros(n_triplet_configs,1);

trip_idx = find(connected_triplet);

for ii = 1:numel(trip_idx)
    j = trip_idx(ii);

    l1 = Qtrip*T1(j,:).' >= cos_alpha_max;
    l2 = Qtrip*T2(j,:).' >= cos_alpha_max;
    l3 = Qtrip*T3(j,:).' >= cos_alpha_max;

    q_union_3 = mean(l1 | l2 | l3);

    triplet_terms(j) = (1-q_union_3)^(N-3);
end

N3_th_geom = nchoosek(N,3)*mean(triplet_terms);

%% ============================================================
%  4. Validation par graphes statiques complets
%% ============================================================

N1_emp = zeros(n_graph_realizations,1);
N2_emp = zeros(n_graph_realizations,1);
N3_emp = zeros(n_graph_realizations,1);
beta0_emp = zeros(n_graph_realizations,1);

for r = 1:n_graph_realizations
    positions = R*sample_walker_delta_unit(N,inc);

    D = squareform(pdist(positions));
    A = sparse((D <= dmax) & (D > 0));

    comp = conncomp(graph(A));
    comp_sizes = accumarray(comp',1);

    N1_emp(r) = sum(comp_sizes==1);
    N2_emp(r) = sum(comp_sizes==2);
    N3_emp(r) = sum(comp_sizes==3);
    beta0_emp(r) = numel(comp_sizes);
end

%% Approximation de beta0 par petites composantes
C_macro = 1;
beta0_th_123 = C_macro + N1_th_local + N2_th_geom + N3_th_geom;

%% Console
fprintf('\n=== Estimation amelioree de N1, N2, N3 ===\n');
fprintf('N = %d, i = %.1f deg, dmax = %.0f km\n',N,inc_deg,dmax);

fprintf('\nN1 satellites isoles\n');
fprintf('Theorie locale latitude : %.4f\n',N1_th_local);
fprintf('Simulation               : %.4f +/- %.4f\n', ...
    mean(N1_emp),std(N1_emp));

fprintf('\nN2 dimeres isoles\n');
fprintf('Theorie geometrique MC   : %.4f\n',N2_th_geom);
fprintf('Simulation               : %.4f +/- %.4f\n', ...
    mean(N2_emp),std(N2_emp));

fprintf('\nN3 trimeres isoles\n');
fprintf('Theorie geometrique MC   : %.4f\n',N3_th_geom);
fprintf('Simulation               : %.4f +/- %.4f\n', ...
    mean(N3_emp),std(N3_emp));

fprintf('\nbeta0\n');
fprintf('Approximation 1+N1+N2+N3 : %.4f\n',beta0_th_123);
fprintf('Simulation               : %.4f +/- %.4f\n', ...
    mean(beta0_emp),std(beta0_emp));

%% Graphiques
figure;
bar([N1_th_local,mean(N1_emp); ...
     N2_th_geom,mean(N2_emp); ...
     N3_th_geom,mean(N3_emp)]);
grid on;
set(gca,'XTickLabel',{'N_1','N_2','N_3'});
ylabel('Nombre moyen de composantes');
legend('Theorie amelioree','Simulation','Location','best');
title('Petites composantes : theorie et simulation');

figure;
histogram(beta0_emp,25,'Normalization','pdf'); hold on;
xline(beta0_th_123,'--','1+N_1+N_2+N_3','LineWidth',1.8);
grid on;
xlabel('\beta_0');
ylabel('Densite');
title('\beta_0 : simulation et approximation par petites composantes');
legend('Simulation','Approximation','Location','best');

figure;
plot(rad2deg(phi_vals),p_link_phi,'LineWidth',2);
grid on;
xlabel('Latitude \phi (deg)');
ylabel('p_{link}(\phi)');
title('Probabilite locale de lien utilisee pour N_1');

%% Sauvegarde
save('N1_N2_N3_walker_delta_results.mat', ...
    'R','inc_deg','inc','dmax','N','alpha_max', ...
    'phi_vals','p_link_phi','N1_th_local', ...
    'N2_th_geom','N3_th_geom','C_macro','beta0_th_123', ...
    'N1_emp','N2_emp','N3_emp','beta0_emp', ...
    'n_pair_configs','n_probe_pair', ...
    'n_triplet_configs','n_probe_triplet', ...
    'n_graph_realizations');

%% ============================================================
%  Fonctions locales
%% ============================================================

function P = sample_walker_delta_unit(M,inc)

Omega = 2*pi*rand(M,1);
u = 2*pi*rand(M,1);

cO = cos(Omega);
sO = sin(Omega);
cu = cos(u);
su = sin(u);
ci = cos(inc);
si = sin(inc);

x = cO.*cu - sO.*su.*ci;
y = sO.*cu + cO.*su.*ci;
z = su.*si;

P = [x y z];
end

function [x,w] = gauss_legendre_interval(n,a,b)

k = (1:n-1)';
beta = k./sqrt(4*k.^2-1);

J = diag(beta,1)+diag(beta,-1);

[V,D] = eig(J);
x0 = diag(D);

[x0,idx] = sort(x0);
V = V(:,idx);

w0 = 2*(V(1,:).^2)';

x = (b-a)/2*x0+(a+b)/2;
w = (b-a)/2*w0;
end
