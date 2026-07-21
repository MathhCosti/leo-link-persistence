clear; clc; close all; rng(4);

%% ============================================================
%  ESTIMATION AMELIOREE DE N1, N2 ET N3
%  Walker Delta a uniformite orbitale
%
%  N1 : quadrature deterministe sur les phases orbitales
%  N2 : quadrature deterministe imbriquee 3D + 2D
%  N3 : quadrature deterministe imbriquee 5D + 2D
%
%  Comparaison avec des graphes statiques simules.
%% ============================================================

%% Parametres physiques
R = 6371 + 550;      % km
inc_deg = 90;
inc = deg2rad(inc_deg);
dmax = 1500;         % km
N = 204;

%% Parametres numeriques
% Pour N1
n_phi = 180;
nQuad = 220;

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
%  2. N2 : dimeres isoles par quadrature deterministe
%
%  N2 = C(N,2)/(2*pi)^3 int du1 du2 dOmega2
%       1_{1~2} [1-q_union,2]^(N-2)
%
%  q_union,2 = 1/(2*pi)^2 int du3 dOmega3
%              1_{3~1 ou 3~2}
%% ============================================================

nQuad_N2_outer = 30;
nQuad_N2_probe = 50;

[uN2,wN2] = gauss_legendre_interval(nQuad_N2_outer,0,2*pi);
[oN2,woN2] = gauss_legendre_interval(nQuad_N2_outer,0,2*pi);
[uP2,wUP2] = gauss_legendre_interval(nQuad_N2_probe,0,2*pi);
[oP2,wOP2] = gauss_legendre_interval(nQuad_N2_probe,0,2*pi);

[U3,O3] = ndgrid(uP2,oP2);
[WU3,WO3] = ndgrid(wUP2,wOP2);
P3probe = walker_delta_position(U3(:),O3(:),inc);
W3 = WU3(:).*WO3(:);

integral_N2 = 0;

for a = 1:nQuad_N2_outer
    P1 = walker_delta_position(uN2(a),0,inc);

    for b = 1:nQuad_N2_outer
        for c = 1:nQuad_N2_outer
            P2 = walker_delta_position(uN2(b),oN2(c),inc);

            if dot(P1,P2) >= cos_alpha_max
                l31 = P3probe*P1.' >= cos_alpha_max;
                l32 = P3probe*P2.' >= cos_alpha_max;

                q_union_2 = sum(W3.*double(l31 | l32))/(2*pi)^2;
                term = (1-q_union_2)^(N-2);
            else
                term = 0;
            end

            integral_N2 = integral_N2 + ...
                wN2(a)*wN2(b)*woN2(c)*term;
        end
    end
end

N2_th_geom = nchoosek(N,2)*integral_N2/(2*pi)^3;

%% ============================================================
%  3. N3 : trimeres isoles par quadrature deterministe
%
%  N3 = C(N,3)/(2*pi)^5 int du1 du2 du3 dOmega2 dOmega3
%       1_{G3 connexe} [1-q_union,3]^(N-3)
%
%  q_union,3 = 1/(2*pi)^2 int du4 dOmega4
%              1_{4~1 ou 4~2 ou 4~3}
%% ============================================================

nQuad_N3_outer = 30;
nQuad_N3_probe = 80;

[uN3,wN3] = gauss_legendre_interval(nQuad_N3_outer,0,2*pi);
[oN3,woN3] = gauss_legendre_interval(nQuad_N3_outer,0,2*pi);
[uP3,wUP3] = gauss_legendre_interval(nQuad_N3_probe,0,2*pi);
[oP3,wOP3] = gauss_legendre_interval(nQuad_N3_probe,0,2*pi);

[U4,O4] = ndgrid(uP3,oP3);
[WU4,WO4] = ndgrid(wUP3,wOP3);
P4probe = walker_delta_position(U4(:),O4(:),inc);
W4 = WU4(:).*WO4(:);

integral_N3 = 0;

for a = 1:nQuad_N3_outer
    P1 = walker_delta_position(uN3(a),0,inc);

    for b = 1:nQuad_N3_outer
        for c = 1:nQuad_N3_outer
            P2 = walker_delta_position(uN3(b),oN3(c),inc);
            e12 = dot(P1,P2) >= cos_alpha_max;

            for d = 1:nQuad_N3_outer
                for e = 1:nQuad_N3_outer
                    P3 = walker_delta_position(uN3(d),oN3(e),inc);

                    e13 = dot(P1,P3) >= cos_alpha_max;
                    e23 = dot(P2,P3) >= cos_alpha_max;

                    connected3 = (e12 + e13 + e23) >= 2;

                    if connected3
                        l41 = P4probe*P1.' >= cos_alpha_max;
                        l42 = P4probe*P2.' >= cos_alpha_max;
                        l43 = P4probe*P3.' >= cos_alpha_max;

                        q_union_3 = sum(W4.*double(l41 | l42 | l43))/(2*pi)^2;
                        term = (1-q_union_3)^(N-3);
                    else
                        term = 0;
                    end

                    integral_N3 = integral_N3 + ...
                        wN3(a)*wN3(b)*woN3(c)*wN3(d)*woN3(e)*term;
                end
            end
        end
    end
end

N3_th_geom = nchoosek(N,3)*integral_N3/(2*pi)^5;

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
fprintf('Quadrature deterministe : %.4f\n',N2_th_geom);
fprintf('Simulation               : %.4f +/- %.4f\n', ...
    mean(N2_emp),std(N2_emp));

fprintf('\nN3 trimeres isoles\n');
fprintf('Quadrature deterministe : %.4f\n',N3_th_geom);
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
    'nQuad_N2_outer','nQuad_N2_probe', ...
    'nQuad_N3_outer','nQuad_N3_probe', ...
    'n_graph_realizations');

%% ============================================================
%  Fonctions locales
%% ============================================================


function P = walker_delta_position(u,Omega,inc)

u = u(:);
Omega = Omega(:);
if isscalar(u) && numel(Omega)>1, u = repmat(u,size(Omega)); end
if isscalar(Omega) && numel(u)>1, Omega = repmat(Omega,size(u)); end

cO = cos(Omega); sO = sin(Omega);
cu = cos(u); su = sin(u);
ci = cos(inc); si = sin(inc);

x = cO.*cu - sO.*su.*ci;
y = sO.*cu + cO.*su.*ci;
z = su.*si;
P = [x y z];
end

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
