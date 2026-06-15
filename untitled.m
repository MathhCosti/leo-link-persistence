clear; clc; close all;

% Paramètres
R_earth = 6371;     % km
h = 550;            % altitude LEO en km (exemple de Starlink)
R = R_earth + h;    % rayon orbital

lambda = 1e-6;      % intensité en satellites / km^2
surface_sphere = 4*pi*R^2;

% Nombre de satellites : processus de Poisson sur la sphère
N = poissrnd(lambda * surface_sphere);

% Tirage uniforme sur la sphère
u = rand(N,1);                  % uniforme sur [0,1]
phi = 2*pi*rand(N,1);            % longitude uniforme
theta = acos(1 - 2*u);           % colatitude corrigée

% Coordonnées cartésiennes
x = R * sin(theta) .* cos(phi);
y = R * sin(theta) .* sin(phi);
z = R * cos(theta);

% Affichage
figure;
scatter3(x, y, z, 20, 'filled');
axis equal;
grid on;
xlabel('x (km)');
ylabel('y (km)');
zlabel('z (km)');
title(sprintf('Processus de Poisson sur sphère LEO, N = %d', N));