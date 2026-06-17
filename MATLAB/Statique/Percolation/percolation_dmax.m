clear; clc; close all;

%% Paramètres

R_E = 6371;               % rayon terrestre [km]

N = 200;                  % nombre de satellites
h = 550;                  % altitude [km]

numTests = 1000;          % simulations Monte-Carlo par point
eta = 0.5;                % percolation si Cmax/N >= eta
delta = 0.05;             % Hoeffding : confiance 1-delta

%% Valeurs de d_max balayées

dmax_vals = linspace(200, 4000, 60);      % [km]

%% Conversion d_max -> alpha_max

r = R_E + h;

ratio = dmax_vals ./ (2*r);
ratio = min(max(ratio, 0), 1);   % sécurité numérique

alpha_vals = 2 * asin(ratio);

%% Calcul de la percolation

res = percolation_alpha_sweep(N, alpha_vals, numTests, eta, delta);

%% Affichage

plot_percolation_curve(dmax_vals, res, ...
    'd_{max} [km]', ...
    sprintf('Percolation finie en fonction de d_{max} : N = %d, h = %.0f km', N, h), ...
    eta);

%% Ajout des seuils théoriques

% Dans ce modèle sphérique, pour deux satellites tirés uniformément
% sur une sphère de rayon orbital r, on a exactement :
%
%   p_link = d_max^2 / (4 r^2)
%
% Le degré moyen vaut donc :
%
%   d_bar = (N-1) p_link
%
% On convertit les seuils exprimés en degré moyen vers un d_max critique.

hold on;

% Seuil de percolation heuristique du graphe géométrique aléatoire 2D
deg_perc = 4.512;

% Seuil de connexité obtenu avec le critère exact E[# isolés] = 1
% Approximation courante : deg_conn ~= log(N)
deg_conn_exact = (N-1) * (1 - N^(-1/(N-1)));
deg_conn_log = log(N);

% Conversion degré moyen -> d_max critique [km]
dmax_perc = 2 * r * sqrt(deg_perc / (N-1));
dmax_conn_exact = 2 * r * sqrt(deg_conn_exact / (N-1));
dmax_conn_log = 2 * r * sqrt(deg_conn_log / (N-1));

% Lignes verticales seulement si elles tombent dans la fenêtre affichée
yl = ylim;

if dmax_perc >= min(dmax_vals) && dmax_perc <= max(dmax_vals)
    xline(dmax_perc, '--', ...
        sprintf('Percolation : %.0f km', dmax_perc), ...
        'LineWidth', 1.5, ...
        'LabelOrientation', 'horizontal', ...
        'LabelVerticalAlignment', 'bottom', ...
        'HandleVisibility', 'off');
end

if dmax_conn_exact >= min(dmax_vals) && dmax_conn_exact <= max(dmax_vals)
    xline(dmax_conn_exact, '--', ...
        sprintf('Connexité : %.0f km', dmax_conn_exact), ...
        'LineWidth', 1.5, ...
        'LabelOrientation', 'horizontal', ...
        'LabelVerticalAlignment', 'top', ...
        'HandleVisibility', 'off');
end

ylim(yl);

%% Affichage info Hoeffding

eps_H = sqrt(log(2/delta) / (2*numTests));

fprintf('Nombre de tests Monte-Carlo par point : %d\n', numTests);
fprintf('Niveau de confiance : %.2f %%\n', 100*(1-delta));
fprintf('Demi-largeur Hoeffding : %.4f\n', eps_H);
