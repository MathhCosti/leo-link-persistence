clear; clc; close all;

%% Paramètres

R_E = 6371;               % rayon terrestre [km]

N = 200;                  % nombre de satellites
dmax = 1500;              % portée de lien [km]

numTests = 1000;          % simulations Monte-Carlo par point
eta = 0.1;               % percolation si Cmax/N >= eta
delta = 0.05;             % Hoeffding : confiance 1-delta

%% Valeurs de h balayées

h_vals = linspace(100, 2000, 60);      % [km]

%% Conversion h -> alpha_max avec dmax fixé

r = R_E + h_vals;

ratio = dmax ./ (2*r);
ratio = min(max(ratio, 0), 1);   % sécurité numérique

alpha_vals = 2 * asin(ratio);

%% Calcul de la percolation

res = percolation_alpha_sweep(N, alpha_vals, numTests, eta, delta);

%% Affichage

plot_percolation_curve(h_vals, res, ...
    'h [km]', ...
    sprintf('Percolation finie en fonction de h : N = %d, d_{max} = %.0f km', N, dmax), ...
    eta);

%% Affichage info Hoeffding

eps_H = sqrt(log(2/delta) / (2*numTests));

fprintf('Nombre de tests Monte-Carlo par point : %d\n', numTests);
fprintf('Niveau de confiance : %.2f %%\n', 100*(1-delta));
fprintf('Demi-largeur Hoeffding : %.4f\n', eps_H);