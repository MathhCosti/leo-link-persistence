% Positions unitaires initiales
r0 = positions ./ vecnorm(positions, 2, 2);  % N x 3

% Création d'une direction tangentielle aléatoire pour chaque satellite
random_vec = randn(N, 3);

% Projection dans le plan tangent à la sphère
tangent = random_vec - dot(random_vec, r0, 2) .* r0;
tangent = tangent ./ vecnorm(tangent, 2, 2);

% Pôle orbital : normal au plan orbital
orbital_poles = cross(r0, tangent, 2);
orbital_poles = orbital_poles ./ vecnorm(orbital_poles, 2, 2);

mu = 398600;          % km^3/s^2, paramètre gravitationnel terrestre
omega = sqrt(mu / R^3);  % rad/s

% Temps de simulation
dt = 10;              % pas de temps en secondes
Tmax = 6000;          % durée totale en secondes
time_values = 0:dt:Tmax;

% Stockage optionnel de la probabilité de routage
P_routing_time = zeros(length(time_values), 1);

for k = 1:length(time_values)

    t = time_values(k);

    % Mouvement orbital
    coswt = cos(omega * t);
    sinwt = sin(omega * t);

    r_t_unit = r0 * coswt + cross(orbital_poles, r0, 2) * sinwt;
    positions_t = R * r_t_unit;

    x_t = positions_t(:,1);
    y_t = positions_t(:,2);
    z_t = positions_t(:,3);

    % Construction du graphe à l'instant t
    D = squareform(pdist(positions_t));
    A = (D <= dmax) & (D > 0);
    G = graph(A);

    % Probabilité de routage à l'instant t
    comp = conncomp(G);
    component_sizes = histcounts(comp, 1:(max(comp)+1));

    P_routing_time(k) = sum(component_sizes .* (component_sizes - 1)) / (N * (N - 1));
end

% Affichage
figure;
plot(time_values, P_routing_time, 'LineWidth', 1.5);
grid on;
xlabel('Temps (s)');
ylabel('P(routage multi-sauts)');
title('Évolution temporelle de la probabilité de routage');