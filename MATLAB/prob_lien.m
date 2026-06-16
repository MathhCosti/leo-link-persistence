clear; clc; close all;

%% Paramètres
numTests = 1e5;                  % nombre de tests Monte-Carlo par valeur de alpha
alpha_vals = linspace(0, pi, 50); % valeurs de alpha_max entre 0 et pi radians

P_sim = zeros(size(alpha_vals));
P_theo = zeros(size(alpha_vals));

%% Boucle sur les valeurs de alpha_max
for k = 1:length(alpha_vals)

    alpha_max = alpha_vals(k);

    % --- Génération uniforme de deux satellites sur la sphère unité ---

    % Satellite 1
    phi1 = 2*pi*rand(numTests,1);
    u1 = 2*rand(numTests,1) - 1;      % cos(theta1) uniforme dans [-1,1]
    theta1 = acos(u1);

    x1 = sin(theta1).*cos(phi1);
    y1 = sin(theta1).*sin(phi1);
    z1 = cos(theta1);

    % Satellite 2
    phi2 = 2*pi*rand(numTests,1);
    u2 = 2*rand(numTests,1) - 1;      % cos(theta2) uniforme dans [-1,1]
    theta2 = acos(u2);

    x2 = sin(theta2).*cos(phi2);
    y2 = sin(theta2).*sin(phi2);
    z2 = cos(theta2);

    % --- Calcul de l'angle entre les deux satellites ---

    dotProduct = x1.*x2 + y1.*y2 + z1.*z2;

    % Sécurité numérique pour éviter acos(1.00000001)
    dotProduct = max(min(dotProduct, 1), -1);

    alpha = acos(dotProduct);

    % --- Estimation Monte-Carlo ---
    P_sim(k) = mean(alpha <= alpha_max);

    % --- Formule théorique ---
    P_theo(k) = (1 - cos(alpha_max)) / 2;
end

%% Affichage des résultats
figure;
plot(alpha_vals, P_theo, 'LineWidth', 2);
hold on;
plot(alpha_vals, P_sim, 'o', 'MarkerSize', 4);
grid on;

xlabel('\alpha_{max} en radians');
ylabel('Probabilité de lien');
legend('Théorie', 'Simulation Monte-Carlo', 'Location', 'northwest');
title('Probabilité d''avoir un lien entre deux satellites');

%% Affichage de l'erreur
figure;
plot(alpha_vals, abs(P_sim - P_theo), 'LineWidth', 2);
grid on;

xlabel('\alpha_{max} en radians');
ylabel('Erreur absolue');
title('Erreur entre simulation et théorie');