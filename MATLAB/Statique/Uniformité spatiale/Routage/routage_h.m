%% Probabilité de routage multi-sauts en fonction de h
% Simulation Monte-Carlo + modèles théoriques :
% 1) lien direct
% 2) composante géante type Erdős-Rényi
% 3) approximation géométrique corrigée par un seuil de percolation spatial

clear; clc; close all;

%% Paramètres
R_earth = 6371;     % km
lambda = 5e-7;      % intensité en satellites / km^2
dmax = 1500;        % km

%% Paramètres de simulation
h_values = linspace(500, 2000, 50);          % altitude LEO en km
nSim = 50;                                   % nombre de simulations par h

P_routing_mean = zeros(size(h_values));
P_routing_std  = zeros(size(h_values));
N_mean         = zeros(size(h_values));

%% Résultats théoriques / approximatifs
P_direct_theory       = zeros(size(h_values));
P_routing_ER_theory   = zeros(size(h_values));
P_routing_geo_theory  = zeros(size(h_values));
k_mean_theory         = zeros(size(h_values));
Nbar_theory           = zeros(size(h_values));

% Seuil géométrique effectif.
% Dans un graphe géométrique 2D, le seuil de percolation apparaît pour un
% degré moyen critique plus grand que 1, typiquement autour de 4 à 5.
k_crit_geo = 4.512;

%% Boucle sur h
for idx = 1:length(h_values)

    h = h_values(idx);
    R = R_earth + h;
    surface_sphere = 4*pi*R^2;
    Nbar = lambda * surface_sphere;

    P_routing_sim = zeros(nSim,1);
    N_sim = zeros(nSim,1);

    for sim = 1:nSim

        %% 1. Nombre de satellites : processus de Poisson sur la sphère
        N = poissrnd(Nbar);
        N_sim(sim) = N;

        if N < 2
            P_routing_sim(sim) = 0;
            continue;
        end

        %% 2. Tirage uniforme des satellites sur la sphère
        u = rand(N,1);
        phi = 2*pi*rand(N,1);
        theta = acos(1 - 2*u);

        %% 3. Coordonnées cartésiennes
        x = R * sin(theta) .* cos(phi);
        y = R * sin(theta) .* sin(phi);
        z = R * cos(theta);

        positions = [x y z];

        %% 4. Construction du graphe géométrique
        D = squareform(pdist(positions));
        A = (D <= dmax) & (D > 0);
        G = graph(A);

        %% 5. Probabilité de routage multi-sauts
        % Deux satellites sont routables s'ils appartiennent à la même
        % composante connexe. Si les composantes ont des tailles s_k :
        % P_routing = sum_k s_k(s_k-1)/(N(N-1)).
        comp = conncomp(G);
        component_sizes = histcounts(comp, 1:(max(comp)+1));
        P_routing = sum(component_sizes .* (component_sizes - 1)) / (N * (N - 1));

        P_routing_sim(sim) = P_routing;
    end

    %% Moyenne sur les simulations
    P_routing_mean(idx) = mean(P_routing_sim);
    P_routing_std(idx)  = std(P_routing_sim);
    N_mean(idx) = mean(N_sim);

    %% Théorie du lien direct
    % Pour deux points uniformes sur une sphère de rayon R, la probabilité
    % que leur distance corde soit <= d vaut d^2/(4R^2), tant que d <= 2R.
    % On tronque aussi à la distance maximale compatible avec la LOS.
    d_LOS = 2*sqrt(R^2 - R_earth^2);
    r_eff = min([dmax, d_LOS, 2*R]);
    p_link = r_eff^2/(4*R^2);

    %% Degré moyen théorique
    % On utilise Nbar plutôt qu'un N fixé, car la simulation tire N selon
    % une loi de Poisson à chaque réalisation.
    kbar = max(Nbar - 1, 0) * p_link;

    %% Modèle 1 : approximation Erdős-Rényi classique
    S_ER = giant_component_fraction_ER(kbar);

    %% Modèle 2 : approximation géométrique corrigée
    % Dans un graphe géométrique, le seuil de percolation est repoussé
    % vers k_crit_geo > 1. On remplace donc le degré moyen kbar par
    % un degré moyen effectif kbar/k_crit_geo dans l'équation de la
    % composante géante.
    
    k_eff = kbar - k_crit_geo + 1;
    
    if k_eff <= 1
        S_geo = 0;
    else
        S_geo = 0.5;  % initialisation non nulle pour converger vers la solution positive
    
        for iter = 1:1000
            S_new = 1 - exp(-k_eff * S_geo);
    
            if abs(S_new - S_geo) < 1e-10
                break;
            end
    
            S_geo = S_new;
        end
    end

    P_direct_theory(idx)      = p_link;
    k_mean_theory(idx)        = kbar;
    Nbar_theory(idx)          = Nbar;
    P_routing_ER_theory(idx)  = S_ER^2;
    P_routing_geo_theory(idx) = S_geo^2;

    fprintf("h = %.1f | Nbar = %.1f | k_th = %.2f | P_sim = %.3f | P_ER = %.3f | P_geo = %.3f\n", ...
        h, Nbar, kbar, P_routing_mean(idx), P_routing_ER_theory(idx), P_routing_geo_theory(idx));
end

%% Affichage des résultats
figure;
plot(h_values, P_routing_mean, 's-', 'LineWidth', 1.5); hold on;
plot(h_values, P_routing_ER_theory, '--', 'LineWidth', 1.5);
plot(h_values, P_routing_geo_theory, '-.', 'LineWidth', 1.8);
plot(h_values, P_direct_theory, ':', 'LineWidth', 1.5);
grid on;

xlabel('altitude h (km)');
ylabel('Probabilité');
title('Probabilité de routage en fonction de h');
legend('Simulation : routage multi-sauts', ...
       'Théorie ER : composante géante', ...
       'Approx. géométrique corrigée', ...
       'Théorie : lien direct', ...
       'Location', 'best');

%% Affichage du degré moyen théorique
figure;
plot(h_values, k_mean_theory, '-', 'LineWidth', 1.5); hold on;
yline(1, ':', 'LineWidth', 1.5);
yline(k_crit_geo, '-.', 'LineWidth', 1.5);
plot(h_values, log(max(Nbar_theory, 1)), '--', 'LineWidth', 1.5);
grid on;

xlabel('altitude h (km)');
ylabel('Degré moyen théorique');
title('Seuils approximatifs : percolation et connexité');
legend('k = (Nbar-1)p_{link}', ...
       'seuil ER : k=1', ...
       'seuil géométrique approx. : k_c', ...
       'seuil connexité approx. : log(Nbar)', ...
       'Location', 'best');

%% Fonction locale : composante géante Erdős-Rényi
function S = giant_component_fraction_ER(k)
    if k <= 1
        S = 0;
        return;
    end

    S = 1e-6;
    for it = 1:1000
        S_new = 1 - exp(-k*S);
        if abs(S_new - S) < 1e-12
            break;
        end
        S = S_new;
    end
end
