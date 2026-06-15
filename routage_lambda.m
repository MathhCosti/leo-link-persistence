% Paramètres
R_earth = 6371;     % km
h = 550;            % altitude LEO en km (exemple de Starlink)
R = R_earth + h;    % rayon orbital

lambda = 9e-7;      % intensité en satellites / km^2
surface_sphere = 4*pi*R^2;

dmax = 1500; % km

%% Paramètres de simulation
lambda_values = linspace(1e-8, 1.2e-6, 25);  % satellites / km^2
nSim = 50;                                  % nombre de simulations par lambda

P_routing_mean = zeros(size(lambda_values));
N_mean = zeros(size(lambda_values));

%% Boucle sur lambda
for idx = 1:length(lambda_values)

    lambda = lambda_values(idx);

    P_routing_sim = zeros(nSim,1);
    N_sim = zeros(nSim,1);

    for sim = 1:nSim

        %% 1. Tirage du nombre de satellites
        N = poissrnd(lambda * surface_sphere);
        N_sim(sim) = N;

        % Cas où il y a trop peu de satellites
        if N < 2
            P_routing_sim(sim) = 0;
            continue;
        end

        %% 2. Tirage uniforme sur la sphère
        u = rand(N,1);
        phi = 2*pi*rand(N,1);
        theta = acos(1 - 2*u);

        x = R * sin(theta) .* cos(phi);
        y = R * sin(theta) .* sin(phi);
        z = R * cos(theta);

        positions = [x y z];

        %% 3. Construction du graphe géométrique
        D = squareform(pdist(positions));

        A = (D <= dmax) & (D > 0);
        G = graph(A);

        %% 4. Probabilité de routage multi-sauts
        comp = conncomp(G);
        component_sizes = histcounts(comp, 1:(max(comp)+1));

        P_routing = sum(component_sizes .* (component_sizes - 1)) / (N * (N - 1));

        %% Stockage
        P_routing_sim(sim) = P_routing;
    end

    %% Moyenne sur les simulations
    P_routing_mean(idx) = mean(P_routing_sim);
    N_mean(idx) = mean(N_sim);

    fprintf("lambda = %.2e | N moyen = %.1f | P_routing = %.3f\n", ...
        lambda, N_mean(idx), P_routing_mean(idx));
end

%% Affichage des résultats
figure;
plot(lambda_values, P_routing_mean, 's-', 'LineWidth', 1.5);
grid on;

xlabel('\lambda (satellites / km^2)');
ylabel('Probabilité');
title('Probabilité de lien de routage en fonction de \lambda');
legend('P(routage multi-sauts)', 'Location', 'best');

hold off;