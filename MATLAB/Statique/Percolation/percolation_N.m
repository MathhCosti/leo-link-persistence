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

%% Ajout des seuils théoriques sur l'axe N
hold on;

% Probabilité de lien fixe pour h et dmax fixés
p_link = (1 - cos(alpha_max)) / 2;

% Seuil de percolation : (N-1) p_link = 4.512
deg_perc = 4.512;
N_perc = 1 + deg_perc / p_link;

if N_perc >= min(N_vals) && N_perc <= max(N_vals)
    xline(N_perc, ':', ...
        sprintf('Seuil percolation: N = %.0f', N_perc), ...
        'LineWidth', 1.5, ...
        'LabelOrientation', 'horizontal', ...
        'LabelVerticalAlignment', 'bottom');
end

% Seuil de connexité exact : p_link = 1 - N^(-1/(N-1))
% On le résout numériquement car N apparaît des deux côtés.
f_conn = @(x) 1 - x.^(-1./(x-1)) - p_link;
N_conn = NaN;

if f_conn(min(N_vals)) * f_conn(max(N_vals)) <= 0
    N_conn = fzero(f_conn, [min(N_vals), max(N_vals)]);
elseif abs(f_conn(min(N_vals))) < abs(f_conn(max(N_vals)))
    % Le seuil est probablement en-dessous de l'intervalle affiché
    N_conn = NaN;
else
    % Le seuil est probablement au-dessus de l'intervalle affiché
    N_conn = NaN;
end

if ~isnan(N_conn) && N_conn >= min(N_vals) && N_conn <= max(N_vals)
    xline(N_conn, ':', ...
        sprintf('Seuil connexité: N = %.0f', N_conn), ...
        'LineWidth', 1.5, ...
        'LabelOrientation', 'horizontal', ...
        'LabelVerticalAlignment', 'top');
end

%% Affichage info Hoeffding

eps_H = sqrt(log(2/delta) / (2*numTests));

fprintf('Nombre de tests Monte-Carlo par point : %d\n', numTests);
fprintf('Niveau de confiance : %.2f %%\n', 100*(1-delta));
fprintf('Demi-largeur Hoeffding : %.4f\n', eps_H);