clear; clc; close all;

%% ============================================================
%  PROBABILITE DE LIEN - WALKER DELTA A UNIFORMITE ORBITALE
%
%  Comparaison entre :
%  - simulation Monte-Carlo de deux satellites Walker Delta ;
%  - formule semi-analytique obtenue par quadrature.
%% ============================================================

%% Parametres
numTests = 1e5;
alpha_vals = linspace(0,pi,50);

inc_deg = 90;
inc = deg2rad(inc_deg);

nQuad = 180;

P_sim = zeros(size(alpha_vals));
P_theo = zeros(size(alpha_vals));
P_sphere = zeros(size(alpha_vals));

%% Boucle sur alpha_max
for k = 1:length(alpha_vals)

    alpha_max = alpha_vals(k);

    %% Simulation Walker Delta
    pos1 = sample_walker_delta(numTests,inc);
    pos2 = sample_walker_delta(numTests,inc);

    dotProduct = sum(pos1.*pos2,2);
    dotProduct = max(min(dotProduct,1),-1);

    alpha = acos(dotProduct);

    P_sim(k) = mean(alpha <= alpha_max);

    %% Nouvelle formule semi-analytique Walker Delta
    P_theo(k) = plink_delta_quadrature(alpha_max,inc,nQuad);

    %% Ancienne reference uniforme sur la sphere
    P_sphere(k) = (1-cos(alpha_max))/2;
end

%% Affichage
figure;
plot(alpha_vals,P_theo,'LineWidth',2); hold on;
plot(alpha_vals,P_sim,'o','MarkerSize',4);
plot(alpha_vals,P_sphere,'--','LineWidth',1.5);
grid on;

xlabel('\alpha_{max} en radians');
ylabel('Probabilite de lien');
legend('Theorie Walker Delta','Simulation Monte-Carlo', ...
       'Uniforme sphere','Location','northwest');
title(sprintf('Probabilite de lien - Walker Delta - i = %.1f deg',inc_deg));

%% Erreurs
figure;
plot(alpha_vals,abs(P_sim-P_theo),'LineWidth',2); hold on;
plot(alpha_vals,abs(P_sim-P_sphere),'--','LineWidth',1.5);
grid on;

xlabel('\alpha_{max} en radians');
ylabel('Erreur absolue');
legend('Erreur simulation / theorie Delta', ...
       'Erreur simulation / theorie sphere', ...
       'Location','best');
title('Erreur sur la probabilite de lien');

save('prob_lien_walker_delta_results.mat', ...
    'alpha_vals','inc_deg','nQuad','P_sim','P_theo','P_sphere');

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
