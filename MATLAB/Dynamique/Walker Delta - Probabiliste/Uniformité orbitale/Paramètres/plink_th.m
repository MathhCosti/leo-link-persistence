clear; clc; close all;

%% ============================================================
%  p_link WALKER DELTA PAR REDUCTION SEMI-ANALYTIQUE
%
%  Modele :
%      Omega ~ U(0,2*pi)
%      u     ~ U(0,2*pi)
%      inclinaison i fixe
%
%  L'integrale sur DeltaOmega est resolue analytiquement.
%  Il reste une quadrature deterministe 2D sur (u1,u2).
%% ============================================================

%% Parametres
R_earth = 6371;      % km
h = 550;             % km
R = R_earth + h;

inc_deg = 58;
inc = deg2rad(inc_deg);

dmax = 1500;         % km

% Ordres de quadrature a tester pour verifier la convergence
quad_orders = [50 100 150 200 300 400];

%% Seuil de lien
alpha_max = 2*asin(min(dmax/(2*R),1));
cmax = cos(alpha_max);

ci = cos(inc);
si = sin(inc);

%% Reference uniforme sur la sphere
p_link_sphere = (1-cos(alpha_max))/2;

%% Calcul par quadrature de Gauss-Legendre
p_quad = zeros(size(quad_orders));

for q = 1:numel(quad_orders)
    n = quad_orders(q);

    % Noeuds et poids sur [0,2*pi]
    [u,w] = gauss_legendre_interval(n,0,2*pi);

    u1 = u(:);
    u2 = u(:).';

    c1 = cos(u1);
    s1 = sin(u1);
    c2 = cos(u2);
    s2 = sin(u2);

    % Produit scalaire :
    % Cdot = A*cos(DeltaOmega) + B*sin(DeltaOmega) + C
    A = c1.*c2 + ci^2.*s1.*s2;
    B = ci.*(s1.*c2 - c1.*s2);
    C = si^2.*s1.*s2;

    rho = sqrt(A.^2 + B.^2);

    % Condition :
    % rho*cos(DeltaOmega-delta) + C >= cmax
    qarg = zeros(size(rho));
    mask = rho > 1e-14;
    qarg(mask) = (cmax - C(mask))./rho(mask);

    % Fraction exacte des DeltaOmega satisfaisant la condition
    g = zeros(size(rho));

    g(mask & qarg <= -1) = 1;

    middle = mask & qarg > -1 & qarg < 1;
    g(middle) = acos(qarg(middle))/pi;

    % Cas degenere rho = 0
    g(~mask) = double(C(~mask) >= cmax);

    % Quadrature 2D :
    % p_link = 1/(2*pi)^2 * int int g(u1,u2) du1 du2
    W = w(:)*w(:).';
    p_quad(q) = sum(g.*W,'all')/(2*pi)^2;

    fprintf('Ordre %4d : p_link = %.10f\n',n,p_quad(q));
end

p_link_delta = p_quad(end);

%% Resultats
fprintf('\n=== Resultat semi-analytique ===\n');
fprintf('Inclinaison                         : %.2f deg\n',inc_deg);
fprintf('Rayon orbital                       : %.2f km\n',R);
fprintf('Distance maximale                   : %.2f km\n',dmax);
fprintf('Angle central maximal               : %.8f rad\n',alpha_max);
fprintf('p_link Walker Delta (quadrature)    : %.10f\n',p_link_delta);
fprintf('p_link uniforme sphere              : %.10f\n',p_link_sphere);
fprintf('Ecart relatif sphere -> Delta       : %.2f %%\n', ...
    100*(p_link_delta-p_link_sphere)/p_link_delta);

%% Convergence
figure;
plot(quad_orders,p_quad,'o-','LineWidth',1.6);
grid on;
xlabel('Ordre de quadrature par dimension');
ylabel('p_{link}^{Delta}');
title('Convergence de la quadrature semi-analytique');

%% Exemple : consequences sur les aretes et le degre
N = 204;

E_delta = nchoosek(N,2)*p_link_delta;
k_delta = (N-1)*p_link_delta;

E_sphere = nchoosek(N,2)*p_link_sphere;
k_sphere = (N-1)*p_link_sphere;

fprintf('\nPour N = %d :\n',N);
fprintf('Aretes moyennes Delta               : %.3f\n',E_delta);
fprintf('Degre moyen Delta                   : %.6f\n',k_delta);
fprintf('Aretes moyennes sphere              : %.3f\n',E_sphere);
fprintf('Degre moyen sphere                  : %.6f\n',k_sphere);

save('plink_results.mat', ...
    'R','h','inc_deg','inc','dmax','alpha_max','cmax', ...
    'quad_orders','p_quad','p_link_delta','p_link_sphere', ...
    'N','E_delta','k_delta','E_sphere','k_sphere');

%% ============================================================
%  Fonction locale : quadrature de Gauss-Legendre sur [a,b]
%% ============================================================
function [x,w] = gauss_legendre_interval(n,a,b)

    % Matrice de Jacobi de Legendre
    k = (1:n-1)';
    beta = k./sqrt(4*k.^2-1);

    J = diag(beta,1) + diag(beta,-1);

    % Valeurs propres = noeuds, premiers coefficients = poids
    [V,D] = eig(J);
    x0 = diag(D);

    [x0,idx] = sort(x0);
    V = V(:,idx);

    w0 = 2*(V(1,:).^2)';

    % Passage de [-1,1] a [a,b]
    x = (b-a)/2*x0 + (a+b)/2;
    w = (b-a)/2*w0;
end
