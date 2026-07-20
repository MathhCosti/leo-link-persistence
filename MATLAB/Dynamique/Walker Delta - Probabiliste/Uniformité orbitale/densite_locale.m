clear; clc; close all;

%% DENSITE LOCALE WALKER DELTA : THEORIE ET SIMULATION
%
% Hypotheses orbitales :
%   Omega ~ U(0,2*pi)
%   u     ~ U(0,2*pi)
%   inclinaison commune i
%
% La latitude verifie : sin(phi) = sin(i)*sin(u).

%% Parametres physiques
R_earth = 6371;      % km
h = 550;             % km
R = R_earth + h;     % rayon orbital en km

inc_deg = 58;
inc = deg2rad(inc_deg);

%% Parametres de constellation et de simulation
N = 204;
n_realizations = 10000;   % nombre de constellations independantes
n_bins = 80;              % nombre de bandes de latitude
rng(1);                   % reproductibilite

%% ============================================================
% 1. DENSITE THEORIQUE CONTINUE
%% ============================================================

n_phi = 2000;
eps_phi = 1e-6;
phi = linspace(-inc + eps_phi, inc - eps_phi, n_phi);

% Densite de probabilite de latitude
f_phi_theory = cos(phi) ./ ...
    (pi * sqrt(sin(inc)^2 - sin(phi).^2));

% Densite surfacique locale theorique
lambda_theory = N ./ ...
    (2*pi^2*R^2 * sqrt(sin(inc)^2 - sin(phi).^2));

%% ============================================================
% 2. DENSITE EMPIRIQUE PAR BANDES DE LATITUDE
%% ============================================================

% Generation des phases orbitales uniformes
u_samples = 2*pi*rand(N*n_realizations,1);

% Latitude induite par le mouvement orbital
phi_samples = asin(sin(inc).*sin(u_samples));

% Bornes et centres des bandes
phi_edges = linspace(-inc,inc,n_bins+1);
phi_centers = 0.5*(phi_edges(1:end-1) + phi_edges(2:end));

% Nombre total d'occurrences dans chaque bande
counts_total = histcounts(phi_samples,phi_edges);

% Nombre moyen de satellites par constellation et par bande
mean_count_per_bin = counts_total / n_realizations;

% Aire exacte de chaque bande de latitude
% A([phi_a,phi_b]) = 2*pi*R^2*(sin(phi_b)-sin(phi_a))
area_bins = 2*pi*R^2 .* ...
    (sin(phi_edges(2:end)) - sin(phi_edges(1:end-1)));

% Densite empirique moyenne dans chaque bande
lambda_empirical = mean_count_per_bin ./ area_bins;

%% ============================================================
% 3. VALEUR THEORIQUE MOYENNE DANS CHAQUE BANDE
%% ============================================================

% Fonction de repartition analytique de la latitude :
% F_phi(x) = 1/2 + asin(sin(x)/sin(i))/pi, pour |x| <= i.
F_edges = 0.5 + asin(sin(phi_edges)./sin(inc))/pi;
F_edges(1) = 0;
F_edges(end) = 1;

% Probabilite theorique d'appartenir a chaque bande
probability_bins_theory = diff(F_edges);

% Densite theorique moyenne exacte dans chaque bande
lambda_theory_bins = N*probability_bins_theory ./ area_bins;

%% ============================================================
% 4. DENSITES MOYENNES DE REFERENCE
%% ============================================================

surface_band = 4*pi*R^2*sin(inc);
lambda_band_mean = N/surface_band;
lambda_sphere_mean = N/(4*pi*R^2);

density_ratio_theory = lambda_theory/lambda_band_mean;
density_ratio_empirical = lambda_empirical/lambda_band_mean;

%% ============================================================
% 5. VERIFICATIONS ET ERREURS
%% ============================================================

% Reconstruction de N a partir des densites par bandes
N_reconstructed_theory = sum(lambda_theory_bins .* area_bins);
N_reconstructed_empirical = sum(lambda_empirical .* area_bins);

% Erreurs entre theorie et simulation sur les valeurs moyennees par bande
absolute_error = lambda_empirical - lambda_theory_bins;
rmse = sqrt(mean(absolute_error.^2));
relative_l1_error = ...
    sum(abs(absolute_error).*area_bins) / ...
    sum(lambda_theory_bins.*area_bins);

%% ============================================================
% 6. COMPARAISON SUR LE MEME GRAPHE
%% ============================================================

figure;
plot(rad2deg(phi),lambda_theory,'LineWidth',2); hold on;
plot(rad2deg(phi_centers),lambda_theory_bins,'--','LineWidth',1.8);
plot(rad2deg(phi_centers),lambda_empirical,'o','MarkerSize',4, ...
    'LineWidth',1.1);

yline(lambda_band_mean,'-.', ...
    sprintf('Moyenne bande = %.3e',lambda_band_mean), ...
    'LineWidth',1.3);
yline(lambda_sphere_mean,':', ...
    sprintf('Moyenne sphere = %.3e',lambda_sphere_mean), ...
    'LineWidth',1.3);

grid on;
xlabel('Latitude \phi (deg)');
ylabel('\lambda_{\Delta}(\phi) (satellites/km^2)');
title(sprintf(['Densite locale Walker Delta : theorie et simulation ' ...
    '(N = %d, i = %.1f deg)'],N,inc_deg));
legend('Theorie continue', ...
       'Theorie moyennee par bande', ...
       'Simulation empirique', ...
       'Densite moyenne sur la bande', ...
       'Densite moyenne sur la sphere', ...
       'Location','best');
hold off;

%% Facteur de concentration
figure;
plot(rad2deg(phi),density_ratio_theory,'LineWidth',2); hold on;
plot(rad2deg(phi_centers),density_ratio_empirical,'o', ...
    'MarkerSize',4,'LineWidth',1.1);
yline(1,'--','Moyenne de bande','LineWidth',1.3);
grid on;
xlabel('Latitude \phi (deg)');
ylabel('\lambda_{\Delta}(\phi)/\bar{\lambda}_{bande}');
title('Facteur de concentration locale : theorie et simulation');
legend('Theorie continue','Simulation empirique', ...
    'Moyenne de bande','Location','best');
hold off;

%% Densite de probabilite de latitude
figure;
plot(rad2deg(phi),f_phi_theory,'LineWidth',2); hold on;

% Histogramme normalise en densite de probabilite
histogram(rad2deg(phi_samples),rad2deg(phi_edges), ...
    'Normalization','pdf','DisplayStyle','stairs','LineWidth',1.2);

grid on;
xlabel('Latitude \phi (deg)');
ylabel('Densite de probabilite');
title('Distribution de latitude : theorie et simulation');

% Conversion de la densite theorique rad^{-1} vers deg^{-1}
% pour etre compatible avec l'histogramme exprime en degres.
cla;
plot(rad2deg(phi),f_phi_theory*pi/180,'LineWidth',2); hold on;
histogram(rad2deg(phi_samples),rad2deg(phi_edges), ...
    'Normalization','pdf','DisplayStyle','stairs','LineWidth',1.2);
grid on;
xlabel('Latitude \phi (deg)');
ylabel('f_{\phi}^{\Delta}(\phi) (deg^{-1})');
title('Distribution de latitude : theorie et simulation');
legend('Theorie','Simulation empirique','Location','best');
hold off;

%% ============================================================
% 7. AFFICHAGE CONSOLE
%% ============================================================

fprintf('\n=== Densite locale Walker Delta ===\n');
fprintf('Nombre de satellites N                    : %d\n',N);
fprintf('Nombre de realisations                    : %d\n',n_realizations);
fprintf('Inclinaison                               : %.2f deg\n',inc_deg);
fprintf('Rayon orbital                             : %.2f km\n',R);
fprintf('Surface de la bande accessible            : %.6e km^2\n',surface_band);
fprintf('Densite moyenne sur la bande              : %.6e sat/km^2\n', ...
    lambda_band_mean);
fprintf('Densite moyenne sur toute la sphere       : %.6e sat/km^2\n', ...
    lambda_sphere_mean);
fprintf('N reconstruit par la theorie par bandes   : %.6f\n', ...
    N_reconstructed_theory);
fprintf('N reconstruit par la simulation           : %.6f\n', ...
    N_reconstructed_empirical);
fprintf('RMSE theorie/empirique                     : %.6e sat/km^2\n',rmse);
fprintf('Erreur L1 relative ponderee par les aires : %.4f %%\n', ...
    100*relative_l1_error);

%% ============================================================
% 8. SAUVEGARDE
%% ============================================================

save('densite_locale_walker_delta_comparaison.mat', ...
    'R_earth','h','R','inc_deg','inc','N', ...
    'n_realizations','n_bins', ...
    'phi','f_phi_theory','lambda_theory', ...
    'phi_edges','phi_centers','area_bins', ...
    'lambda_theory_bins','lambda_empirical', ...
    'lambda_band_mean','lambda_sphere_mean', ...
    'N_reconstructed_theory','N_reconstructed_empirical', ...
    'rmse','relative_l1_error');

fprintf(['Resultats sauvegardes dans ' ...
    'densite_locale_walker_delta_comparaison.mat\n']);
