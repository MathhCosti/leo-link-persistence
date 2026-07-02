clear; clc; close all;

%% ============================================================
%  NOMBRE DE LIENS INTER-SATELLITES TEMPOREL
%  Version : mean empirique + minimum théorique + maximum théorique
%
%  Sorties principales :
%  - num_edges_all(iter,k) : nombre de liens à l'instant t_k
%  - mean_edges(t)        : moyenne empirique sur les itérations
%  - L_min_theory_uniform : minimum théorique, configuration uniforme
%  - L_max_theory_polar   : maximum théorique, concentration polaire
%% ============================================================

%% Paramètres physiques
R_earth = 6371;      % km
h = 550;             % km
R = R_earth + h;     % rayon orbital

mu = 398600;              % km^3/s^2
omega = sqrt(mu / R^3);   % vitesse angulaire orbitale rad/s

%% Paramètres du processus de Poisson
lambda = 4e-7;       % satellites / km^2
surface_sphere = 4*pi*R^2;
N_mean_theory = lambda * surface_sphere;

%% Paramètres des liens et du temps
dmax = 1500;      % km
dt = 60;          % pas temporel en secondes
Tmax = 12000;     % durée totale de simulation

time_values = 0:dt:Tmax;
Nt = length(time_values);

%% Nombre d'itérations Monte-Carlo
n_iter = 100;

%% ============================================================
%  GRANDEURS THÉORIQUES
%% ============================================================

% Distance maximale réellement utilisable : seuil radio + contrainte LOS.
d_LOS = 2*sqrt(R^2 - R_earth^2);
d_eff = min(dmax, d_LOS);

% Angle central correspondant à la distance euclidienne d_eff :
% d_eff = 2 R sin(alpha_max/2)
alpha_max = 2*asin(d_eff/(2*R));

% Probabilité qu'une paire de satellites uniformes soit liée.
% La calotte accessible autour d'un satellite a une aire
% S_cap = 2*pi*R^2*(1 - cos(alpha_max)), donc :
p_link_uniform = (1 - cos(alpha_max))/2;
% Équivalent exact avec la distance de corde : p_link_uniform = d_eff^2/(4*R^2).

% Minimum théorique : répartition uniforme sur la sphère.
% On utilise N_mean_theory pour obtenir une valeur déterministe indépendante
% des fluctuations de N entre les itérations.
L_min_theory_uniform = N_mean_theory*(N_mean_theory - 1)/2 * p_link_uniform;

A_demi_bande_max = 2*pi*R^2*sin(0.5*alpha_max);
mu_calotte_max = lambda * A_demi_bande_max;

%% Analyse par strates polaires
% On découpe la région polaire en couronnes successives :
% [0, alpha_max/2], [alpha_max/2, alpha_max], etc.
% Pour chaque strate, on estime le nombre moyen de satellites arrivant dans
% la strate et son degré moyen local, mais sans filtrage par percolation.
% Pas des strates : ici alpha_max/2.
% Les strates sont construites jusqu'à pi/2 sans jamais dépasser.
% Aucune strate n'est rejetée par un critère de percolation.
beta_step_strates = alpha_max;
beta_max_strates = pi/2;

strates_polaires = strates_polaires(lambda, R, alpha_max, ...
    'beta_step', beta_step_strates, ...
    'beta_max', beta_max_strates, ...
    'verbose', true);

% Maximum théorique raffiné : somme des contributions des strates actives.
% Pour chaque strate active i, on utilise :
%   L_i ≈ (1/2) * mu_i * k_i
% pour un pôle, puis on multiplie par 2 pour les deux pôles.
% Cette version utilise uniquement la formule du degre moyen local,
% sans correction particuliere pour la calotte centrale.
[L_max_theory_strates, details_liens_strates] = ...
    liens_theoriques(strates_polaires);

%% Stockage Monte-Carlo
num_edges_all = zeros(n_iter, Nt);
N_all = zeros(n_iter, 1);

%% Stockage empirique par strate polaire active
n_strates_active = height(strates_polaires.active_table);
links_strate_internal_all = zeros(n_iter, Nt, n_strates_active);
links_strate_incident_all = zeros(n_iter, Nt, n_strates_active);
n_sat_strate_all = zeros(n_iter, Nt, n_strates_active);

%% ============================================================
%  BOUCLE MONTE-CARLO
%% ============================================================

for it = 1:n_iter

    %% Nombre de satellites
    N = poissrnd(lambda * surface_sphere);
    N_all(it) = N;

    fprintf('Itération %d / %d : N = %d satellites\n', it, n_iter, N);

    %% Génération uniforme des positions initiales sur la sphère
    u = rand(N,1);
    phi = 2*pi*rand(N,1);
    theta = acos(1 - 2*u);

    %% Sens de rotation défini par un plan séparateur passant par les pôles
    rotation_sign = ones(N, 1);
    rotation_sign(phi >= pi) = -1;

    %% Evolution temporelle
    for k = 1:Nt

        t = time_values(k);

        %% Mouvement orbital
        phi_t = phi;
        theta_t = theta + rotation_sign * omega * t;

        x_t = R * sin(theta_t) .* cos(phi_t);
        y_t = R * sin(theta_t) .* sin(phi_t);
        z_t = R * cos(theta_t);

        positions_t = [x_t y_t z_t];

        %% Matrice des distances et liens inter-satellites
        D = squareform(pdist(positions_t));
        A = (D <= d_eff) & (D > 0);

        %% Nombre de liens inter-satellites
        num_edges_all(it,k) = nnz(triu(A,1));

        %% Liens empiriques par strate polaire active
        if n_strates_active > 0
            [L_incident_strate, n_sat_strate] = ...
                liens_empiriques(A, z_t, R, strates_polaires);

            links_strate_incident_all(it,k,:) = reshape(L_incident_strate, 1, 1, []);
            n_sat_strate_all(it,k,:) = reshape(n_sat_strate, 1, 1, []);
        end
    end
end

%% ============================================================
%  STATISTIQUES EMPIRIQUES
%% ============================================================

mean_edges = mean(num_edges_all, 1);      % mean empirique temporel
mean_edges_global = mean(mean_edges);     % moyenne globale temps + itérations

%% Sinusoïde théorique approximée à partir du minimum et du maximum théoriques
% Le nombre de liens présente deux pics par orbite : la période des liens est
% donc approximativement T_links = T_orb/2 = pi/omega.
T_orb = 2*pi/omega;
T_links = T_orb/2;

L_sin_center = 0.5 * (L_max_theory_strates + L_min_theory_uniform);
L_sin_amp    = 0.5 * (L_max_theory_strates - L_min_theory_uniform);

% On choisit la phase pour que la sinusoïde soit minimale à t = 0,
% comme dans la dynamique orbitale simulée.
L_sin_theory = L_sin_center - L_sin_amp * cos(2*pi*time_values/T_links);

% Sécurité numérique : évite de petites valeurs hors bornes dues à l'arrondi.
L_sin_theory = max(min(L_sin_theory, L_max_theory_strates), L_min_theory_uniform);

% Instant où le nombre moyen total de liens est maximal.
[mean_edges_peak, idx_peak_edges] = max(mean_edges);
t_peak_edges = time_values(idx_peak_edges);

% Statistiques empiriques par strate, moyennées sur les itérations.
if n_strates_active > 0
    mean_links_internal_strate_t = reshape(mean(links_strate_internal_all, 1), Nt, n_strates_active);
    mean_links_incident_strate_t = reshape(mean(links_strate_incident_all, 1), Nt, n_strates_active);
    mean_n_sat_strate_t = reshape(mean(n_sat_strate_all, 1), Nt, n_strates_active);

    L_internal_peak = mean_links_internal_strate_t(idx_peak_edges, :)';
    L_incident_peak = mean_links_incident_strate_t(idx_peak_edges, :)';
    n_sat_peak = mean_n_sat_strate_t(idx_peak_edges, :)';

    L_internal_max_t = max(mean_links_internal_strate_t, [], 1)';
    L_incident_max_t = max(mean_links_incident_strate_t, [], 1)';

    empirical_strates_comparison = strates_polaires.active_table(:, ...
        {'index', 'beta_in_deg', 'beta_out_deg', 'mu', 'k_mean_local'});
    empirical_strates_comparison.N_emp_peak = n_sat_peak;
    empirical_strates_comparison.L_internal_emp_peak = L_internal_peak;
    empirical_strates_comparison.L_incident_emp_peak = L_incident_peak;
    empirical_strates_comparison.L_internal_emp_max_t = L_internal_max_t;
    empirical_strates_comparison.L_incident_emp_max_t = L_incident_max_t;

    if ~isempty(details_liens_strates) && ...
            ismember('L_two_poles', details_liens_strates.Properties.VariableNames)

        % On aligne les deux tables par indice de strate.
        [is_found, loc] = ismember(empirical_strates_comparison.index, details_liens_strates.index);

        % Initialisation avec NaN au cas où certaines strates empiriques
        % ne sont pas présentes dans la table théorique.
        empirical_strates_comparison.L_theory_degree = NaN(height(empirical_strates_comparison),1);

        empirical_strates_comparison.L_theory_degree(is_found) = ...
            details_liens_strates.L_two_poles(loc(is_found));

        empirical_strates_comparison.ratio_incident_peak_theory = ...
            empirical_strates_comparison.L_incident_emp_peak ./ ...
            max(empirical_strates_comparison.L_theory_degree, eps);

        empirical_strates_comparison.ratio_internal_peak_theory = ...
            empirical_strates_comparison.L_internal_emp_peak ./ ...
            max(empirical_strates_comparison.L_theory_degree, eps);
    end
else
    mean_links_internal_strate_t = [];
    mean_links_incident_strate_t = [];
    mean_n_sat_strate_t = [];
    empirical_strates_comparison = table();
end

mean_N = mean(N_all);
min_N = min(N_all);
max_N = max(N_all);

fprintf('\nNombre théorique moyen de satellites : %.2f\n', N_mean_theory);
fprintf('Nombre empirique moyen de satellites : %.2f\n', mean_N);
fprintf('Nombre minimal empirique de satellites : %d\n', min_N);
fprintf('Nombre maximal empirique de satellites : %d\n', max_N);

fprintf('\nalpha_max = %.4f rad = %.2f deg\n', alpha_max, rad2deg(alpha_max));
fprintf('p_link uniforme = %.6f\n', p_link_uniform);
fprintf('Minimum théorique uniforme : %.2f liens\n', L_min_theory_uniform);
fprintf('Maximum théorique par strates poissonien : %.2f liens\n', L_max_theory_strates);
fprintf('Mean empirique globale : %.2f liens\n', mean_edges_global);
fprintf('Surface demi-bande max : %.2e km^2\n', A_demi_bande_max);
fprintf('Nombre moyen de satellites par calotte au maximum : %.2f\n', mu_calotte_max);

fprintf('\n--- Analyse par strates polaires ---\n');
fprintf('Pas angulaire des strates : %.4f rad = %.2f deg\n', ...
    beta_step_strates, rad2deg(beta_step_strates));
fprintf('Strates considerees jusqu''a beta_max = %.4f rad = %.2f deg\n', ...
    beta_max_strates, rad2deg(beta_max_strates));
fprintf('Nombre de strates considerees : %d\n', height(strates_polaires.active_table));

fprintf('Maximum strates, version poissonienne : %.2f liens\n', L_max_theory_strates);
fprintf('Période orbitale : %.2f s\n', T_orb);
fprintf('Période théorique des liens : %.2f s\n', T_links);
fprintf('Sinusoïde théorique : centre = %.2f, amplitude = %.2f\n', L_sin_center, L_sin_amp);
fprintf('Contributions par strate active :\n');
disp(details_liens_strates);

%% ============================================================
%  GRAPHE : MEAN EMPIRIQUE + MIN/MAX THÉORIQUES
%% ============================================================

figure;
hold on;

plot(time_values, mean_edges, 'k', 'LineWidth', 2.2);

plot(time_values, L_sin_theory, 'm--', ...
    'LineWidth', 2.0);

yline(L_min_theory_uniform, 'g--', ...
    sprintf('Minimum th. uniforme = %.1f', L_min_theory_uniform), ...
    'LineWidth', 1.8, ...
    'LabelHorizontalAlignment', 'left');

yline(L_max_theory_strates, 'b--', ...
    sprintf('Maximum th. strates = %.1f', L_max_theory_strates), ...
    'LineWidth', 2.0, ...
    'LabelHorizontalAlignment', 'left');

yline(mean_edges_global, 'k:', ...
    sprintf('Mean empirique globale = %.1f', mean_edges_global), ...
    'LineWidth', 1.5, ...
    'LabelHorizontalAlignment', 'left');

grid on;
xlabel('Temps (s)');
ylabel('Nombre de liens inter-satellites');
title(sprintf('Liens inter-satellites : mean empirique et bornes théoriques — %d itérations', n_iter));
legend('Mean empirique temporel', ...
       'Sinusoïde théorique min/max', ...
       'Minimum théorique uniforme', ...
       'Maximum théorique par strates', ...
       'Mean empirique globale', ...
       'Location', 'best');

%% ============================================================
%  GRAPHE : DEGRÉ MOYEN PAR STRATE POLAIRE
%% ============================================================

figure;
hold on;

beta_mid_deg = rad2deg(strates_polaires.all_table.beta_mid);
bar(beta_mid_deg, strates_polaires.all_table.k_mean_local);
grid on;
xlabel('Milieu de la strate polaire (degrés depuis le pôle)');
ylabel('Degré moyen local estimé');
title('Degré moyen local par strate polaire');

%% ============================================================
%  GRAPHE : COMPARAISON EMPIRIQUE / THÉORIE PAR STRATE
%% ============================================================

if n_strates_active > 0 && ~isempty(empirical_strates_comparison)
    figure;
    hold on;

    x = empirical_strates_comparison.index;

    if ismember('L_theory_degree', empirical_strates_comparison.Properties.VariableNames)
        Y = [empirical_strates_comparison.L_incident_emp_peak, ...
             empirical_strates_comparison.L_theory_degree];
        bar(x, Y);
        legend('Empirique', ...
               'Théorie', ...
               'Location', 'best');
    else
        Y = empirical_strates_comparison.L_incident_emp_peak;
        bar(x, Y);
        legend('Empirique', ...
               'Location', 'best');
    end

    grid on;
    xlabel('Indice de strate active');
    ylabel('Nombre de liens');
    title(sprintf('Liens par strate au pic empirique moyen, t = %.0f s', t_peak_edges));
end


%% ============================================================
%  GRAPHE 3D : STRATES POLAIRES CONSERVÉES
%% ============================================================

% plot_strates_polaires_3D(R, alpha_max, strates_polaires, ...
%     'show_inactive', true, ...
%     'sphere_alpha', 0.04, ...
%     'strate_alpha', 0.30);

%% ============================================================
%  SAUVEGARDE
%% ============================================================

save('liens_inter.mat', ...
    'R_earth', 'h', 'R', 'mu', 'omega', ...
    'lambda', 'surface_sphere', 'N_mean_theory', ...
    'dmax', 'd_LOS', 'd_eff', 'alpha_max', 'p_link_uniform', ...
    'dt', 'Tmax', 'time_values', 'n_iter', ...
    'A_demi_bande_max', 'mu_calotte_max', ...
    'beta_step_strates', 'beta_max_strates', 'strates_polaires', ...
    'details_liens_strates', 'empirical_strates_comparison', ...
    'links_strate_incident_all', 'n_sat_strate_all', ...
    'mean_links_internal_strate_t', 'mean_links_incident_strate_t', 'mean_n_sat_strate_t', ...
    'idx_peak_edges', 't_peak_edges', 'mean_edges_peak', ...
    'T_orb', 'T_links', 'L_sin_center', 'L_sin_amp', 'L_sin_theory', ...
    'L_max_theory_strates', ...
    'N_all', 'num_edges_all', 'mean_edges', 'mean_edges_global', ...
    'L_min_theory_uniform');

fprintf('\nAnalyse terminée.\n');
fprintf('Résultats sauvegardés dans liens_inter_satellites_theorie_min_max_mean_empirique.mat\n');
