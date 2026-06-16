clear; clc; close all;

%% Paramètres

R_E = 6371;               % rayon terrestre [km]

h = 550;                  % altitude [km]
dmax = 1500;              % portée de lien [km]

numTests = 1000;          % simulations Monte-Carlo par point
eta = 0.5;                % percolation si Cmax/N >= eta
delta = 0.05;             % Hoeffding : confiance 1-delta

useEarthOccultation = false;
% false : seulement contrainte de distance inter-satellite
% true  : ajoute une contrainte de visibilité ligne de vue, non bloquée par la Terre

%% Valeurs de N balayées

N_vals = unique(round(linspace(50, 800, 40)));

%% Conversion dmax, h -> alpha_max

r = R_E + h;

ratio = dmax ./ (2*r);
ratio = min(max(ratio, 0), 1);   % sécurité numérique

alpha_max = 2 * asin(ratio);

%% Calcul de la percolation

res = percolation_N_sweep(N_vals, alpha_max, numTests, eta, delta);

%% Affichage

plot_percolation_curve(N_vals, res, ...
    'N', ...
    sprintf('Percolation finie en fonction de N : h = %.0f km, d_{max} = %.0f km', h, dmax), ...
    eta);

%% Affichage info Hoeffding

eps_H = sqrt(log(2/delta) / (2*numTests));

fprintf('Nombre de tests Monte-Carlo par point : %d\n', numTests);
fprintf('Niveau de confiance : %.2f %%\n', 100*(1-delta));
fprintf('Demi-largeur Hoeffding : %.4f\n', eps_H);