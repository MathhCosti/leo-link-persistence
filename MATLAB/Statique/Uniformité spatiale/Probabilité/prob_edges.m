clear; clc; close all;

%% Paramètres

N = 200;                         % nombre de satellites
numTests = 1000;                 % nombre de simulations par valeur de alpha
alpha_vals = linspace(0.01, pi/2, 60); % valeurs de alpha_max testées

E_edges_sim = zeros(size(alpha_vals));
E_edges_theo = zeros(size(alpha_vals));

%% Boucle sur les valeurs de alpha_max

for k = 1:length(alpha_vals)

    alpha_max = alpha_vals(k);

    edges_tests = zeros(numTests, 1);

    for t = 1:numTests

        %% Génération uniforme de N satellites sur la sphère unité

        % Méthode : vecteurs gaussiens normalisés
        positions = randn(N, 3);
        positions = positions ./ vecnorm(positions, 2, 2);

        %% Calcul des angles entre tous les satellites

        cosAlpha = positions * positions';

        % Sécurité numérique pour éviter acos(1.00000001)
        cosAlpha = max(min(cosAlpha, 1), -1);

        alpha = acos(cosAlpha);

        %% Matrice d'adjacence

        A = (alpha <= alpha_max);

        % Suppression des boucles i -> i
        A = A & ~eye(N);

        %% Nombre d'arêtes dans le graphe

        % Comme A est symétrique, chaque arête est comptée deux fois :
        % une fois en A(i,j), une fois en A(j,i).
        nb_edges = sum(A(:)) / 2;

        edges_tests(t) = nb_edges;

    end

    %% Moyenne simulée du nombre d'arêtes

    E_edges_sim(k) = mean(edges_tests);

    %% Formule théorique

    p_link = (1 - cos(alpha_max)) / 2;

    E_edges_theo(k) = nchoosek(N, 2) * p_link;

end

%% Affichage théorie / simulation

figure;
plot(alpha_vals, E_edges_theo, 'LineWidth', 2);
hold on;
plot(alpha_vals, E_edges_sim, 'o', 'MarkerSize', 4);
grid on;

xlabel('\alpha_{max} en radians');
ylabel('Nombre moyen d''arêtes');
legend('Théorie', 'Simulation Monte-Carlo', 'Location', 'northwest');
title('Vérification de E[|E|] = C(N,2) p_{link}');

%% Affichage de l'erreur absolue

figure;
plot(alpha_vals, abs(E_edges_sim - E_edges_theo), 'LineWidth', 2);
grid on;

xlabel('\alpha_{max} en radians');
ylabel('Erreur absolue');
title('Erreur entre simulation et théorie');