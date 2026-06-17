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

%% Ajout des seuils théoriques sur l'axe h
hold on;

% Seuils en degré moyen
deg_perc = 4.512;
deg_conn = (N-1) * (1 - N^(-1/(N-1)));

% Pour dmax fixé : deg(h) = (N-1) dmax^2 / (4 (R_E+h)^2)
% donc h_seuil = dmax/(2 sqrt(deg/(N-1))) - R_E
h_perc = dmax / (2 * sqrt(deg_perc/(N-1))) - R_E;
h_conn = dmax / (2 * sqrt(deg_conn/(N-1))) - R_E;

if h_perc >= min(h_vals) && h_perc <= max(h_vals)
    xline(h_perc, ':', ...
        sprintf('Seuil percolation: h = %.0f km', h_perc), ...
        'LineWidth', 1.5, ...
        'LabelOrientation', 'horizontal', ...
        'LabelVerticalAlignment', 'bottom');
end

if h_conn >= min(h_vals) && h_conn <= max(h_vals)
    xline(h_conn, ':', ...
        sprintf('Seuil connexité: h = %.0f km', h_conn), ...
        'LineWidth', 1.5, ...
        'LabelOrientation', 'horizontal', ...
        'LabelVerticalAlignment', 'top');
end

%% Affichage info Hoeffding

eps_H = sqrt(log(2/delta) / (2*numTests));

fprintf('Nombre de tests Monte-Carlo par point : %d\n', numTests);
fprintf('Niveau de confiance : %.2f %%\n', 100*(1-delta));
fprintf('Demi-largeur Hoeffding : %.4f\n', eps_H);