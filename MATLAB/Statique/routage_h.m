% Paramètres
R_earth = 6371;     % km
lambda = 9e-7;      % intensité en satellites / km^2

dmax = 1500;        % km

%% Paramètres de simulation
h_values = linspace(500, 2000, 50);         % altitude LEO en km
nSim = 50;                                  % nombre de simulations par h

P_routing_mean = zeros(size(h_values));
N_mean = zeros(size(h_values));

%% Boucle sur h
for idx = 1:length(h_values)

    h = h_values(idx);
    R = R_earth + h;

    surface_sphere = 4*pi*R^2;

    if idx == 1

        % Nombre de satellites : processus de Poisson sur la sphère
        N = poissrnd(lambda * surface_sphere);

        % Tirage uniforme sur la sphère
        u = rand(N,1);                   % uniforme sur [0,1]
        phi = 2*pi*rand(N,1);            % longitude uniforme
        theta = acos(1 - 2*u);           % colatitude corrigée

    end

    % Coordonnées cartésiennes
    x = R * sin(theta) .* cos(phi);
    y = R * sin(theta) .* sin(phi);
    z = R * cos(theta);

    positions = [x y z];

    P_routing_sim = zeros(nSim,1);
    N_sim = zeros(nSim,1);

    for sim = 1:nSim

        %% 1. Construction du graphe géométrique
        D = squareform(pdist(positions));

        A = (D <= dmax) & (D > 0);
        G = graph(A);

        %% 2. Probabilité de routage multi-sauts
        comp = conncomp(G);
        component_sizes = histcounts(comp, 1:(max(comp)+1));

        P_routing = sum(component_sizes .* (component_sizes - 1)) / (N * (N - 1));

        %% Stockage
        P_routing_sim(sim) = P_routing;
    end

    %% Moyenne sur les simulations
    P_routing_mean(idx) = mean(P_routing_sim);
    N_mean(idx) = mean(N_sim);

    fprintf("h = %.2e | N moyen = %.1f | P_routing = %.3f\n", ...
        h, N_mean(idx), P_routing_mean(idx));
end

%% Affichage des résultats
figure;
plot(h_values, P_routing_mean, 's-', 'LineWidth', 1.5);
grid on;

xlabel('altitude (km)');
ylabel('Probabilité');
title('Probabilité de lien de routage en fonction de altitude');
legend('P(routage multi-sauts)', 'Location', 'best');

hold off;