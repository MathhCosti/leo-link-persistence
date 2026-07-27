clear; clc; close all;

%% Paramètres
N = 100;                         % nombre de satellites
numTests = 1000;                 % nombre de simulations par valeur de alpha
alpha_vals = linspace(0, pi, 50); % valeurs de alpha_max en radians

E_deg_sim = zeros(size(alpha_vals));
E_deg_theo = zeros(size(alpha_vals));

%% Boucle sur les valeurs de alpha_max
for k = 1:length(alpha_vals)

    alpha_max = alpha_vals(k);

    deg_mean_tests = zeros(numTests, 1);

    for t = 1:numTests

        %% Génération uniforme de N satellites sur la sphère unité

        % Méthode : vecteurs gaussiens normalisés
        positions = randn(N, 3);
        positions = positions ./ vecnorm(positions, 2, 2);

        %% Calcul des angles entre tous les satellites

        cosAlpha = positions * positions';

        % Sécurité numérique
        cosAlpha = max(min(cosAlpha, 1), -1);

        alpha = acos(cosAlpha);

        %% Matrice d'adjacence

        A = (alpha <= alpha_max);

        % Suppression des liens d'un satellite vers lui-même
        A = A & ~eye(N);

        %% Degré de chaque satellite

        deg = sum(A, 2);

        % Moyenne des degrés dans cette réalisation
        deg_mean_tests(t) = mean(deg);

    end

    %% Espérance simulée du degré
    E_deg_sim(k) = mean(deg_mean_tests);

    %% Formule théorique
    p_link = (1 - cos(alpha_max)) / 2;
    E_deg_theo(k) = (N - 1) * p_link;

end

%% Affichage comparaison théorie / simulation

figure;
plot(alpha_vals, E_deg_theo, 'LineWidth', 2);
hold on;
plot(alpha_vals, E_deg_sim, 'o', 'MarkerSize', 4);
grid on;

xlabel('\alpha_{max} en radians');
ylabel('Espérance du degré');
legend('Théorie', 'Simulation Monte-Carlo', 'Location', 'northwest');
title('Espérance du degré moyen en fonction de \alpha_{max}');

%% Erreur absolue

figure;
plot(alpha_vals, abs(E_deg_sim - E_deg_theo), 'LineWidth', 2);
grid on;

xlabel('\alpha_{max} en radians');
ylabel('Erreur absolue');
title('Erreur entre simulation et théorie');