%% Paramètres de simulation
dmax_values = linspace(600, 2000, 50);      % altitude LEO en km
nSim = 50;                                  % nombre de simulations par dmax

P_routing_mean = zeros(size(dmax_values));
N_mean = zeros(size(dmax_values));

R_earth = 6371;     % km
h = 550;            % altitude LEO en km (exemple de Starlink)
R = R_earth + h;    % rayon orbital

lambda = 9e-7;      % intensité en satellites / km^2
surface_sphere = 4*pi*R^2;

% Nombre de satellites : processus de Poisson sur la sphère
N = poissrnd(lambda * surface_sphere);

% Tirage uniforme sur la sphère
u = rand(N,1);                  % uniforme sur [0,1]
phi = 2*pi*rand(N,1);            % longitude uniforme
theta = acos(1 - 2*u);           % colatitude corrigée

% Coordonnées cartésiennes
x = R * sin(theta) .* cos(phi);
y = R * sin(theta) .* sin(phi);
z = R * cos(theta);

positions = [x y z];

%% Boucle sur h
for idx = 1:length(dmax_values)

    dmax = dmax_values(idx);

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

    fprintf("dmax = %.2e | N moyen = %.1f | P_routing = %.3f\n", ...
        dmax, N_mean(idx), P_routing_mean(idx));
end

%% Affichage des résultats
figure;
plot(dmax_values, P_routing_mean, 's-', 'LineWidth', 1.5);
grid on;

xlabel('d_{max} (km)');
ylabel('Probabilité');
title('Probabilité de lien de routage en fonction de d_{max}');
legend('P(routage multi-sauts)', 'Location', 'best');

hold off;