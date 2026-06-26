clear; clc; close all;

%% Paramètres

N = 200;                  % nombre de satellites
numTests = 1000;          % nombre de simulations Monte-Carlo par alpha
eta = 0.5;                % seuil de percolation : Cmax/N >= eta
delta = 0.05;             % niveau d'erreur pour Hoeffding : confiance 1-delta

alpha_vals = linspace(0.01, pi/2, 60); % seuils angulaires testés

%% Calcul de la percolation et des bornes

res = percolation_alpha_sweep(N, alpha_vals, numTests, eta, delta);

%% Seuils théoriques en degré moyen

% Percolation 2D : degré moyen critique du modèle de Gilbert
deg_perc = 4.512;

% Connexité : seuil exact associé au critère E[# satellites isolés] = 1
deg_conn = (N-1) * (1 - N^(-1/(N-1)));

% Conversion degré moyen -> probabilité de lien
p_perc = deg_perc / (N-1);
p_conn = deg_conn / (N-1);

% Conversion probabilité de lien -> seuil angulaire
% p_link = (1 - cos(alpha_max))/2
alpha_perc = NaN;
alpha_conn = NaN;

if p_perc <= 1
    alpha_perc = acos(1 - 2*p_perc);
end

if p_conn <= 1
    alpha_conn = acos(1 - 2*p_conn);
end

%% Approximation par seuils

% On modélise la transition par une sigmoïde centrée sur le seuil de
% percolation. Le seuil de percolation est donc le point où P ~= 1/2.
% On règle la raideur pour que P ~= 1-epsilon au seuil de connexité.
epsilon_threshold = 0.01;

if isfinite(alpha_perc) && isfinite(alpha_conn) && abs(alpha_conn-alpha_perc) > eps
    k_sig = log((1-epsilon_threshold)/epsilon_threshold) / abs(alpha_conn-alpha_perc);
    s_sig = sign(alpha_conn-alpha_perc);
    res.P_threshold_approx = 1 ./ (1 + exp(-s_sig*k_sig*(alpha_vals-alpha_perc)));
else
    res.P_threshold_approx = NaN(size(alpha_vals));
end

%% Affichage principal avec toutes les bornes

plot_percolation_curve(alpha_vals, res, ...
    '\alpha_{max} en radians', ...
    sprintf('Percolation finie'), ...
    eta);

%% Ajout des seuils théoriques

hold on;

if isfinite(alpha_perc) && alpha_perc >= min(alpha_vals) && alpha_perc <= max(alpha_vals)
    xline(alpha_perc, ':', ...
        sprintf('Seuil percolation: \\alpha = %.3f rad', alpha_perc), ...
        'LineWidth', 1.5, ...
        'LabelOrientation', 'horizontal', ...
        'LabelVerticalAlignment', 'bottom', ...
        'DisplayName', 'Seuil percolation');
end

if isfinite(alpha_conn) && alpha_conn >= min(alpha_vals) && alpha_conn <= max(alpha_vals)
    xline(alpha_conn, ':', ...
        sprintf('Seuil connexité: \\alpha = %.3f rad', alpha_conn), ...
        'LineWidth', 1.5, ...
        'LabelOrientation', 'horizontal', ...
        'LabelVerticalAlignment', 'top', ...
        'DisplayName', 'Seuil connexité');
end

legend('show', 'Location', 'southeast');
ylim([0 1]);

%% Affichage de l'erreur Hoeffding

eps_H = sqrt(log(2/delta) / (2*numTests));

fprintf('Nombre de tests Monte-Carlo : %d\n', numTests);
fprintf('Niveau de confiance : %.2f %%\n', 100*(1-delta));
fprintf('Demi-largeur Hoeffding : %.4f\n', eps_H);
fprintf('Seuil percolation théorique : alpha = %.4f rad\n', alpha_perc);
fprintf('Seuil connexité théorique : alpha = %.4f rad\n', alpha_conn);
