clear; clc; close all;

%% Paramètres

N = 200;                  % nombre de satellites
numTests = 1000;          % nombre de simulations Monte-Carlo par alpha
eta = 0.5;                % seuil de percolation : Cmax/N >= eta
delta = 0.05;             % niveau d'erreur pour Hoeffding : confiance 1-delta

alpha_vals = linspace(0.01, pi/2, 60); % seuils angulaires testés

m = ceil(eta * N);        % taille minimale de la composante géante

%% Stockage des résultats

P_perc_hat = zeros(size(alpha_vals));

Hoeffding_low = zeros(size(alpha_vals));
Hoeffding_up  = zeros(size(alpha_vals));

Bound_upper_nonisolated = zeros(size(alpha_vals));
Bound_upper_edges = zeros(size(alpha_vals));
Bound_upper_math = zeros(size(alpha_vals));
Bound_upper_tree = zeros(size(alpha_vals));

Bound_lower_star = zeros(size(alpha_vals));

%% Boucle principale

for k = 1:length(alpha_vals)

    alpha_max = alpha_vals(k);

    perc_events = zeros(numTests, 1);

    for t = 1:numTests

        %% Génération uniforme de N satellites sur la sphère unité

        positions = randn(N, 3);
        positions = positions ./ vecnorm(positions, 2, 2);

        %% Angles entre satellites

        cosAlpha = positions * positions';

        % sécurité numérique
        cosAlpha = max(min(cosAlpha, 1), -1);

        alpha = acos(cosAlpha);

        %% Graphe géométrique aléatoire

        A = (alpha <= alpha_max);
        A = A & ~eye(N);

        G = graph(A);

        %% Taille de la plus grande composante

        comp = conncomp(G);
        comp_sizes = accumarray(comp', 1);

        Cmax = max(comp_sizes);

        %% Événement de percolation finie

        perc_events(t) = (Cmax >= m);

    end

    %% Estimation Monte-Carlo

    P_perc_hat(k) = mean(perc_events);

    %% Borne probabiliste de Hoeffding

    eps_H = sqrt(log(2/delta) / (2*numTests));

    Hoeffding_low(k) = max(0, P_perc_hat(k) - eps_H);
    Hoeffding_up(k)  = min(1, P_perc_hat(k) + eps_H);

    %% Probabilité théorique de lien entre deux satellites

    p_link = (1 - cos(alpha_max)) / 2;

    %% Borne mathématique supérieure 1 :
    % Si une composante de taille >= m existe,
    % alors il existe au moins m satellites non isolés.
    %
    % P(Cmax >= m) <= E[# satellites non isolés] / m

    p_isolated = (1 - p_link)^(N-1);
    E_nonisolated = N * (1 - p_isolated);

    Bound_upper_nonisolated(k) = min(1, E_nonisolated / m);

    %% Borne mathématique supérieure 2 :
    % Si une composante de taille >= m existe,
    % alors le graphe contient au moins m-1 arêtes.
    %
    % P(Cmax >= m) <= E[# arêtes] / (m-1)

    E_edges = nchoosek(N, 2) * p_link;

    Bound_upper_edges(k) = min(1, E_edges / (m - 1));

    %% Borne mathématique supérieure 3 :
    if p_link == 0
        Bound_upper_tree(k) = 0;
    else
        logBound = gammaln(N+1) ...
             - gammaln(m+1) ...
             - gammaln(N-m+1) ...
             + (m-2)*log(m) ...
             + (m-1)*log(p_link);

        Bound_upper_tree(k) = min(1, exp(logBound));
    end

    %% Meilleure des deux bornes supérieures

    Bound_upper_math(k) = min(Bound_upper_nonisolated(k), Bound_upper_edges(k));
    Bound_upper_math(k) = min(Bound_upper_tree(k), Bound_upper_math(k));

    %% Borne mathématique inférieure conservative :
    % Borne inférieure par étoile :
    % Si le satellite 1 a au moins m-1 voisins,
    % alors il existe une composante de taille au moins m.
    
    p_link = (1 - cos(alpha_max)) / 2;
    
    Bound_lower_star(k) = 1 - binocdf(m-2, N-1, p_link);

end

%% Affichage principal

figure;
hold on; grid on;

% Intervalle de Hoeffding
fill([alpha_vals fliplr(alpha_vals)], ...
     [Hoeffding_low fliplr(Hoeffding_up)], ...
     [0.8 0.8 0.8], ...
     'EdgeColor', 'none', ...
     'FaceAlpha', 0.5);

plot(alpha_vals, P_perc_hat, 'LineWidth', 2);

plot(alpha_vals, Bound_upper_math, '--', 'LineWidth', 2);
plot(alpha_vals, Bound_lower_star, '--', 'LineWidth', 2);

xlabel('\alpha_{max} en radians');
ylabel('Probabilité');
title(['Percolation finie : P(C_{max}/N \geq ', num2str(eta), ')']);

legend('Intervalle Hoeffding', ...
       'Estimation Monte-Carlo', ...
       'Borne supérieure mathématique', ...
       'Borne inférieure conservative', ...
       'Location', 'southeast');

ylim([0 1]);

%% Affichage de l'erreur Hoeffding

eps_H = sqrt(log(2/delta) / (2*numTests));

fprintf('Nombre de tests Monte-Carlo : %d\n', numTests);
fprintf('Niveau de confiance : %.2f %%\n', 100*(1-delta));
fprintf('Demi-largeur Hoeffding : %.4f\n', eps_H);