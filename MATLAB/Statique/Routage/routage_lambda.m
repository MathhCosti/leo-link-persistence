%% Probabilité de routage multi-sauts en fonction de lambda
% Simulation Monte-Carlo + modèles théoriques :
% 1) lien direct
% 2) composante géante type Erdős-Rényi
% 3) approximation géométrique corrigée par un seuil de percolation spatial

clear; clc; close all;

%% Paramètres
R_earth = 6371;     % km
h = 550;            % altitude LEO en km
R = R_earth + h;    % rayon orbital

surface_sphere = 4*pi*R^2;
dmax = 1500;        % km

%% Paramètres de simulation
lambda_values = linspace(1e-8, 1.2e-6, 25);  % satellites / km^2
nSim = 50;                                    % nombre de simulations par lambda

P_routing_mean = zeros(size(lambda_values));
P_routing_std  = zeros(size(lambda_values));
N_mean         = zeros(size(lambda_values));

%% Résultats théoriques / approximatifs
P_direct_theory       = zeros(size(lambda_values));   % probabilité de lien direct
P_routing_ER_theory   = zeros(size(lambda_values));   % approximation Erdős-Rényi
P_routing_geo_theory  = zeros(size(lambda_values));   % approximation géométrique corrigée
k_mean_theory         = zeros(size(lambda_values));   % degré moyen

% Seuil géométrique effectif.
% Dans un graphe géométrique 2D, le seuil de percolation apparaît pour un
% degré moyen critique plus grand que 1, typiquement autour de 4 à 5.
k_crit_geo = 4.512;

%% Géométrie du lien sur la sphère
% Deux satellites sont liés si leur distance de corde est <= dmax.
% Pour deux points uniformes sur la sphère de rayon R :
% p_link = surface calotte / surface sphère = dmax^2/(4R^2), si dmax <= 2R.
% Si on ajoute une contrainte LOS Terre, on borne aussi dmax par d_LOS.
d_LOS = 2*sqrt(R^2 - R_earth^2);   % distance max sans occultation terrestre
r_eff = min([dmax, d_LOS, 2*R]);
p_link = r_eff^2/(4*R^2);

%% Boucle sur lambda
for idx = 1:length(lambda_values)

    lambda = lambda_values(idx);

    P_routing_sim = zeros(nSim,1);
    N_sim = zeros(nSim,1);

    for sim = 1:nSim

        %% 1. Nombre de satellites : processus de Poisson sur la sphère
        N = poissrnd(lambda * surface_sphere);
        N_sim(sim) = N;

        if N < 2
            P_routing_sim(sim) = 0;
            continue;
        end

        %% 2. Tirage uniforme des satellites sur la sphère
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

    %% Théorie approchée
    Nbar = lambda * surface_sphere;

    % Degré moyen théorique.
    % On utilise Nbar plutôt qu'un N fixé, car la simulation tire N selon
    % une loi de Poisson à chaque réalisation.
    kbar = max(Nbar - 1, 0) * p_link;

    %% Modèle 1 : approximation Erdős-Rényi classique
    % Ce modèle suppose des arêtes indépendantes. Il déclenche souvent la
    % percolation trop tôt pour un graphe géométrique.
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
    P_routing_ER_theory(idx)  = S_ER^2;
    P_routing_geo_theory(idx) = S_geo^2;

    fprintf("lambda = %.2e | Nbar = %.1f | k_th = %.2f | P_sim = %.3f | P_ER = %.3f | P_geo = %.3f\n", ...
        lambda, Nbar, kbar, P_routing_mean(idx), P_routing_ER_theory(idx), P_routing_geo_theory(idx));
end

%% Affichage des résultats : probabilité en fonction de lambda
figure;
plot(lambda_values, P_routing_mean, 's-', 'LineWidth', 1.5); hold on;
plot(lambda_values, P_routing_ER_theory, '--', 'LineWidth', 1.5);
plot(lambda_values, P_routing_geo_theory, '-.', 'LineWidth', 1.8);
plot(lambda_values, P_direct_theory, ':', 'LineWidth', 1.5);
grid on;

xlabel('\lambda (satellites / km^2)');
ylabel('Probabilité');
title('Probabilité de routage en fonction de \lambda');
legend('Simulation : routage multi-sauts', ...
       'Théorie ER : composante géante', ...
       'Approx. géométrique corrigée', ...
       'Théorie : lien direct', ...
       'Location', 'best');

%% Affichage des résultats : probabilité en fonction du nombre moyen de satellites N
figure;
plot(N_mean, P_routing_mean, 'o-', 'LineWidth', 1.5); hold on;
plot(lambda_values*surface_sphere, P_routing_ER_theory, '--', 'LineWidth', 1.5);
plot(lambda_values*surface_sphere, P_routing_geo_theory, '-.', 'LineWidth', 1.8);
grid on;

xlabel('Nombre moyen de satellites N');
ylabel('Probabilité');
title('Probabilité de routage en fonction de N');
legend('Simulation : routage multi-sauts', ...
       'Théorie ER : composante géante', ...
       'Approx. géométrique corrigée', ...
       'Location', 'best');

%% Affichage du degré moyen théorique
figure;
plot(lambda_values, k_mean_theory, '-', 'LineWidth', 1.5); hold on;
yline(1, ':', 'LineWidth', 1.5);
yline(k_crit_geo, '-.', 'LineWidth', 1.5);
plot(lambda_values, log(max(lambda_values*surface_sphere, 1)), '--', 'LineWidth', 1.5);
grid on;

xlabel('\lambda (satellites / km^2)');
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
