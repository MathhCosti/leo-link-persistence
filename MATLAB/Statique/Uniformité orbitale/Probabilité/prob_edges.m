clear; clc; close all;

%% ============================================================
%  NOMBRE MOYEN D'ARETES - UNIFORMITE ORBITALE
%% ============================================================

%% Parametres
N = 200;
numTests = 1000;
alpha_vals = linspace(0.01,pi/2,60);

inc_deg = 90;
inc = deg2rad(inc_deg);

nQuad = 180;

E_edges_sim = zeros(size(alpha_vals));
E_edges_theo = zeros(size(alpha_vals));
E_edges_sphere = zeros(size(alpha_vals));

%% Boucle sur alpha_max
for k = 1:length(alpha_vals)

    alpha_max = alpha_vals(k);
    edges_tests = zeros(numTests,1);

    for t = 1:numTests

        %% Generation
        positions = sample_walker_delta(N,inc);

        %% Angles entre satellites
        cosAlpha = positions*positions';
        cosAlpha = max(min(cosAlpha,1),-1);
        alpha = acos(cosAlpha);

        %% Matrice d'adjacence
        A = (alpha <= alpha_max);
        A = A & ~eye(N);

        %% Nombre d'aretes
        edges_tests(t) = sum(A(:))/2;
    end

    %% Moyenne simulee
    E_edges_sim(k) = mean(edges_tests);

    %% Nouvelle formule
    p_link_delta = plink_delta_quadrature(alpha_max,inc,nQuad);
    E_edges_theo(k) = nchoosek(N,2)*p_link_delta;

    %% Ancienne formule sphere
    p_link_sphere = (1-cos(alpha_max))/2;
    E_edges_sphere(k) = nchoosek(N,2)*p_link_sphere;
end

%% Affichage
figure;
plot(alpha_vals,E_edges_theo,'LineWidth',2); hold on;
plot(alpha_vals,E_edges_sim,'o','MarkerSize',4);
plot(alpha_vals,E_edges_sphere,'--','LineWidth',1.5);
grid on;

xlabel('\alpha_{max} en radians');
ylabel('Nombre moyen d''aretes');
legend('Theorie Uniforme orbite','Simulation Monte-Carlo', ...
       'Uniforme sphere','Location','northwest');
title(sprintf('E[|E|] - Uniformité orbitale - i = %.1f deg',inc_deg));

%% Erreurs
figure;
plot(alpha_vals,abs(E_edges_sim-E_edges_theo),'LineWidth',2); hold on;
plot(alpha_vals,abs(E_edges_sim-E_edges_sphere),'--','LineWidth',1.5);
grid on;

xlabel('\alpha_{max} en radians');
ylabel('Erreur absolue');
legend('Erreur simulation / theorie Delta', ...
       'Erreur simulation / theorie sphere', ...
       'Location','best');
title('Erreur sur le nombre moyen d''aretes');

save('prob_edges_walker_delta_results.mat', ...
    'N','alpha_vals','inc_deg','nQuad', ...
    'E_edges_sim','E_edges_theo','E_edges_sphere');

%% ============================================================
%  FONCTIONS LOCALES
%% ============================================================

function positions = sample_walker_delta(N, inc)
Omega = 2*pi*rand(N,1);
u = 2*pi*rand(N,1);

cO = cos(Omega);
sO = sin(Omega);
cu = cos(u);
su = sin(u);
ci = cos(inc);
si = sin(inc);

x = cO.*cu - sO.*su.*ci;
y = sO.*cu + cO.*su.*ci;
z = su.*si;

positions = [x y z];
end

function p = plink_delta_quadrature(alpha_max, inc, nQuad)
% Probabilite de lien Walker Delta avec uniformite orbitale.
% L'integrale sur DeltaOmega est resolue analytiquement.
% Il reste une quadrature 2D sur (u1,u2).

cmax = cos(alpha_max);
ci = cos(inc);
si = sin(inc);

[u,w] = gauss_legendre_interval(nQuad,0,2*pi);

u1 = u(:);
u2 = u(:).';

c1 = cos(u1);
s1 = sin(u1);
c2 = cos(u2);
s2 = sin(u2);

A = c1.*c2 + ci^2.*s1.*s2;
B = ci.*(s1.*c2 - c1.*s2);
C = si^2.*s1.*s2;

rho = sqrt(A.^2 + B.^2);

g = zeros(size(rho));
mask = rho > 1e-14;

q = zeros(size(rho));
q(mask) = (cmax - C(mask))./rho(mask);

g(mask & q <= -1) = 1;

middle = mask & q > -1 & q < 1;
g(middle) = acos(q(middle))/pi;

g(~mask) = double(C(~mask) >= cmax);

W = w(:)*w(:).';
p = sum(g.*W,'all')/(2*pi)^2;
end

function [x,w] = gauss_legendre_interval(n,a,b)
k = (1:n-1)';
beta = k./sqrt(4*k.^2-1);

J = diag(beta,1) + diag(beta,-1);

[V,D] = eig(J);
x0 = diag(D);

[x0,idx] = sort(x0);
V = V(:,idx);

w0 = 2*(V(1,:).^2)';

x = (b-a)/2*x0 + (a+b)/2;
w = (b-a)/2*w0;
end
