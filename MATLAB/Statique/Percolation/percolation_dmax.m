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

%% Seuils théoriques sur l'axe dmax

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

deg_perc = 4.512;
deg_conn = (N-1) * (1 - N^(-1/(N-1)));

dmax_perc = 2 * r * sqrt(deg_perc / (N-1));
dmax_conn = 2 * r * sqrt(deg_conn / (N-1));

%% Approximation par seuils

epsilon_threshold = 0.01;

if isfinite(dmax_perc) && isfinite(dmax_conn) && abs(dmax_conn-dmax_perc) > eps
    k_sig = log((1-epsilon_threshold)/epsilon_threshold) / abs(dmax_conn-dmax_perc);
    s_sig = sign(dmax_conn-dmax_perc);
    res.P_threshold_approx = 1 ./ (1 + exp(-s_sig*k_sig*(dmax_vals-dmax_perc)));
else
    res.P_threshold_approx = NaN(size(dmax_vals));
end

%% Affichage

plot_percolation_curve(dmax_vals, res, ...
    'd_{max} [km]', ...
    sprintf('Percolation finie en fonction de d_{max} : N = %d, h = %.0f km', N, h), ...
    eta);

%% Ajout des seuils théoriques

hold on;

if dmax_perc >= min(dmax_vals) && dmax_perc <= max(dmax_vals)
    xline(dmax_perc, ':', ...
        sprintf('Seuil percolation: d_{max} = %.0f km', dmax_perc), ...
        'LineWidth', 1.5, ...
        'LabelOrientation', 'horizontal', ...
        'LabelVerticalAlignment', 'bottom', ...
        'DisplayName', 'Seuil percolation');
end

if dmax_conn >= min(dmax_vals) && dmax_conn <= max(dmax_vals)
    xline(dmax_conn, ':', ...
        sprintf('Seuil connexité: d_{max} = %.0f km', dmax_conn), ...
        'LineWidth', 1.5, ...
        'LabelOrientation', 'horizontal', ...
        'LabelVerticalAlignment', 'top', ...
        'DisplayName', 'Seuil connexité');
end

legend('show', 'Location', 'southeast');
ylim([0 1]);

%% Affichage info Hoeffding

eps_H = sqrt(log(2/delta) / (2*numTests));

fprintf('Nombre de tests Monte-Carlo par point : %d\n', numTests);
fprintf('Niveau de confiance : %.2f %%\n', 100*(1-delta));
fprintf('Demi-largeur Hoeffding : %.4f\n', eps_H);
fprintf('Seuil percolation théorique : dmax = %.2f km\n', dmax_perc);
fprintf('Seuil connexité théorique : dmax = %.2f km\n', dmax_conn);
