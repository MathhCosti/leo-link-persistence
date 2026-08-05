%% verif_loi_distances_liens.m
% Verification empirique de la loi des distances conditionnellement
% a l'existence d'un lien, puis comparaison au modele local uniforme :
%
%   f_D(d) = 2d/dmax^2,       0 <= d <= dmax
%   F_D(d) = (d/dmax)^2.
%
% Le script utilise les graphes sauvegardes dans
% analysis_temp_results.mat.
%
% Il produit :
%   1) PDF empirique globale et PDF theorique ;
%   2) CDF empirique globale et CDF theorique ;
%   3) PDF empirique par latitude ;
%   4) densite empirique au bord f_D(dmax) par latitude ;
%   5) probabilite empirique d'appartenir a une couronne proche de dmax.
%
% Sortie :
%   verif_loi_distances_liens_results.mat

clear; clc; close all;

%% ============================================================
%  1. Parametres numeriques
%% ============================================================

n_time_samples = 300;     % nombre maximal d'instants analyses
n_distance_bins = 50;     % classes pour les PDF
n_latitude_bins = 24;     % classes de latitude du milieu du lien

% Epaisseur de la couronne utilisee pour estimer f_D(dmax).
% Elle est exprimee comme une fraction de dmax.
boundary_fraction = 0.05;

% Nombre minimal de liens dans une tranche de latitude pour afficher
% une estimation locale fiable.
min_links_per_lat_bin = 100;

%% ============================================================
%  2. Chargement
%% ============================================================

script_dir = fileparts(mfilename('fullpath'));

analysis_candidates = {
    fullfile(script_dir,'analysis_temp_results.mat')
    fullfile(script_dir,'..','analysis_temp_results.mat')
};

analysis_file = '';

for k = 1:numel(analysis_candidates)
    if isfile(analysis_candidates{k})
        analysis_file = analysis_candidates{k};
        break;
    end
end

if isempty(analysis_file)
    error('Fichier analysis_temp_results.mat introuvable.');
end

A = load(analysis_file);

required_fields = {'Positions','Adjacency','dmax'};

for k = 1:numel(required_fields)
    if ~isfield(A,required_fields{k})
        error('analysis_temp_results.mat doit contenir %s.', ...
            required_fields{k});
    end
end

Positions = A.Positions;
Adjacency = A.Adjacency;
dmax = double(A.dmax);

Nt = min(numel(Positions),numel(Adjacency));

if Nt < 1
    error('Aucun instant temporel disponible.');
end

R = mean(vecnorm(double(Positions{1}),2,2));

%% ============================================================
%  3. Choix des instants et grilles
%% ============================================================

n_time_samples = min(n_time_samples,Nt);
time_indices = unique(round(linspace(1,Nt,n_time_samples)));

% Inclinaison / latitude maximale observee.
all_lat_first = asin(max(min( ...
    double(Positions{1}(:,3))/R,1),-1));

lat_max = max(abs(all_lat_first));

if isfield(A,'inc')
    lat_max = max(lat_max,double(A.inc));
end

latitude_edges = linspace(-lat_max,lat_max,n_latitude_bins+1);
latitude_centers = ...
    0.5*(latitude_edges(1:end-1)+latitude_edges(2:end));

distance_edges = linspace(0,dmax,n_distance_bins+1);
distance_centers = ...
    0.5*(distance_edges(1:end-1)+distance_edges(2:end));

distance_bin_widths = diff(distance_edges);

boundary_width = boundary_fraction*dmax;
boundary_lower = dmax-boundary_width;

%% ============================================================
%  4. Extraction des distances de tous les liens
%% ============================================================

all_link_distances = [];
all_link_latitudes = [];

link_count_per_lat_bin = zeros(1,n_latitude_bins);
boundary_count_per_lat_bin = zeros(1,n_latitude_bins);

distance_hist_per_lat_bin = ...
    zeros(n_latitude_bins,n_distance_bins);

for it = 1:numel(time_indices)

    t = time_indices(it);

    X = double(Positions{t});
    U = X./vecnorm(X,2,2);

    adjacency_t = logical(spones(Adjacency{t}));
    adjacency_t = adjacency_t | adjacency_t.';
    adjacency_t(1:size(adjacency_t,1)+1:end) = false;

    [ii,jj] = find(triu(adjacency_t,1));

    if isempty(ii)
        continue;
    end

    % Distance euclidienne (corde) entre satellites.
    link_vector = X(ii,:)-X(jj,:);
    link_distance = vecnorm(link_vector,2,2);

    % Latitude du milieu geometrique normalise du lien.
    midpoint = U(ii,:)+U(jj,:);
    midpoint_norm = vecnorm(midpoint,2,2);

    valid_midpoint = midpoint_norm > 1e-12;
    midpoint(valid_midpoint,:) = ...
        midpoint(valid_midpoint,:) ./ midpoint_norm(valid_midpoint);

    link_latitude = nan(numel(ii),1);
    link_latitude(valid_midpoint) = asin(max(min( ...
        midpoint(valid_midpoint,3),1),-1));

    % Protection contre de tres faibles erreurs numeriques.
    valid_link = ...
        isfinite(link_distance) ...
        & isfinite(link_latitude) ...
        & link_distance <= dmax*(1+1e-10);

    link_distance = min(link_distance(valid_link),dmax);
    link_latitude = link_latitude(valid_link);

    all_link_distances = [all_link_distances;link_distance]; %#ok<AGROW>
    all_link_latitudes = [all_link_latitudes;link_latitude]; %#ok<AGROW>

    lat_bin = discretize(link_latitude,latitude_edges);

    for b = 1:n_latitude_bins
        mask_b = lat_bin == b;

        if ~any(mask_b)
            continue;
        end

        d_b = link_distance(mask_b);

        link_count_per_lat_bin(b) = ...
            link_count_per_lat_bin(b)+numel(d_b);

        boundary_count_per_lat_bin(b) = ...
            boundary_count_per_lat_bin(b) ...
            + nnz(d_b >= boundary_lower);

        distance_hist_per_lat_bin(b,:) = ...
            distance_hist_per_lat_bin(b,:) ...
            + histcounts(d_b,distance_edges);
    end
end

n_links_total = numel(all_link_distances);

if n_links_total == 0
    error('Aucun lien extrait des instants analyses.');
end

%% ============================================================
%  5. Loi empirique globale et loi uniforme locale
%% ============================================================

% PDF empirique globale.
global_counts = histcounts(all_link_distances,distance_edges);
pdf_emp_global = ...
    global_counts ...
    ./ (n_links_total*distance_bin_widths);

% PDF theorique utilisee dans le modele.
pdf_uniform_local = ...
    2*distance_centers/dmax^2;

% CDF empirique aux centres de classes.
cdf_emp_global = cumsum(global_counts)/n_links_total;

% CDF theorique.
cdf_uniform_local = ...
    (distance_edges(2:end)/dmax).^2;

% Densite empirique au bord estimee par la derniere couronne.
boundary_probability_emp_global = ...
    nnz(all_link_distances >= boundary_lower)/n_links_total;

boundary_density_emp_global = ...
    boundary_probability_emp_global/boundary_width;

boundary_density_theory = 2/dmax;

% Probabilite theorique exacte de la couronne.
boundary_probability_theory_exact = ...
    1-(boundary_lower/dmax)^2;

% Approximation lineaire utilisee lorsque la couronne est mince.
boundary_probability_theory_linear = ...
    2*boundary_width/dmax;

%% ============================================================
%  6. Statistiques locales en latitude
%% ============================================================

pdf_emp_by_latitude = nan(n_latitude_bins,n_distance_bins);

for b = 1:n_latitude_bins
    if link_count_per_lat_bin(b) <= 0
        continue;
    end

    pdf_emp_by_latitude(b,:) = ...
        distance_hist_per_lat_bin(b,:) ...
        ./ (link_count_per_lat_bin(b)*distance_bin_widths);
end

boundary_probability_emp_by_latitude = safe_divide( ...
    boundary_count_per_lat_bin,link_count_per_lat_bin);

boundary_density_emp_by_latitude = ...
    boundary_probability_emp_by_latitude/boundary_width;

boundary_density_ratio_emp_theory = ...
    boundary_density_emp_by_latitude/boundary_density_theory;

%% ============================================================
%  7. Diagnostics statistiques
%% ============================================================

% Erreur de CDF de type Kolmogorov-Smirnov sur les distances normalisees.
x_sorted = sort(all_link_distances/dmax);
n = numel(x_sorted);
cdf_emp_sorted = (1:n)'/n;
cdf_th_sorted = x_sorted.^2;

ks_distance = max(abs(cdf_emp_sorted-cdf_th_sorted));

% Moments.
mean_distance_emp = mean(all_link_distances);
mean_distance_theory = 2*dmax/3;

second_moment_emp = mean(all_link_distances.^2);
second_moment_theory = dmax^2/2;

% RMSE de la PDF globale sur les centres.
rmse_pdf_global = sqrt(mean( ...
    (pdf_emp_global-pdf_uniform_local).^2));

%% ============================================================
%  8. Figures
%% ============================================================

figure;
bar(distance_centers/dmax, ...
    pdf_emp_global*dmax,1, ...
    'FaceAlpha',0.35, ...
    'DisplayName','PDF empirique');
hold on;
plot(distance_centers/dmax, ...
    pdf_uniform_local*dmax, ...
    'LineWidth',2.2, ...
    'DisplayName','PDF uniforme locale : 2x');
grid on;
xlabel('Distance normalisee x=d/d_{max}');
ylabel('Densite normalisee');
title('Loi des distances conditionnellement a l''existence d''un lien');
legend('Location','best');
hold off;

figure;
stairs(distance_edges(2:end)/dmax, ...
    cdf_emp_global, ...
    'LineWidth',2, ...
    'DisplayName','CDF empirique');
hold on;
plot(distance_edges(2:end)/dmax, ...
    cdf_uniform_local, ...
    '--','LineWidth',2, ...
    'DisplayName','CDF uniforme locale : x^2');
grid on;
xlabel('Distance normalisee x=d/d_{max}');
ylabel('F_D(d)');
title('Fonction de repartition des distances de lien');
legend('Location','best');
hold off;

figure;
imagesc( ...
    distance_centers/dmax, ...
    rad2deg(latitude_centers), ...
    pdf_emp_by_latitude*dmax);
axis xy;
colorbar;
xlabel('Distance normalisee x=d/d_{max}');
ylabel('Latitude du milieu du lien (deg)');
title('PDF empirique locale des distances');

figure;
hold on;
valid_lat = ...
    link_count_per_lat_bin >= min_links_per_lat_bin ...
    & isfinite(boundary_density_emp_by_latitude);

plot(rad2deg(latitude_centers(valid_lat)), ...
    boundary_density_emp_by_latitude(valid_lat)*dmax, ...
    'LineWidth',2, ...
    'DisplayName','d_{max} f_{D|lien,\phi}^{emp}(d_{max})');

yline(boundary_density_theory*dmax, ...
    '--','LineWidth',2, ...
    'DisplayName','Modele uniforme local : 2');

grid on;
xlabel('Latitude du milieu du lien (deg)');
ylabel('Densite au bord normalisee');
title('Densite des distances au voisinage de d_{max}');
legend('Location','best');
hold off;

figure;
hold on;
plot(rad2deg(latitude_centers(valid_lat)), ...
    boundary_probability_emp_by_latitude(valid_lat), ...
    'LineWidth',2, ...
    'DisplayName','Probabilite empirique');

yline(boundary_probability_theory_exact, ...
    '--','LineWidth',2, ...
    'DisplayName',sprintf( ...
        'Theorie exacte = %.4f', ...
        boundary_probability_theory_exact));

yline(boundary_probability_theory_linear, ...
    ':','LineWidth',2, ...
    'DisplayName',sprintf( ...
        'Approximation lineaire = %.4f', ...
        boundary_probability_theory_linear));

grid on;
xlabel('Latitude du milieu du lien (deg)');
ylabel(sprintf('P(D \\geq %.2f d_{max} | lien,\\phi)', ...
    1-boundary_fraction));
title('Probabilite d''etre proche du seuil de rupture');
legend('Location','best');
hold off;

figure;
plot(rad2deg(latitude_centers), ...
    link_count_per_lat_bin, ...
    'LineWidth',1.8);
grid on;
xlabel('Latitude du milieu du lien (deg)');
ylabel('Nombre de liens echantillonnes');
title('Echantillonnage des liens selon la latitude');

%% ============================================================
%  9. Affichage console
%% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' VERIFICATION DE LA LOI DES DISTANCES DES LIENS\n');
fprintf('============================================================\n');
fprintf('Fichier analyse                    : %s\n',analysis_file);
fprintf('Instants analyses                  : %d\n',numel(time_indices));
fprintf('Nombre total de liens              : %d\n',n_links_total);
fprintf('dmax                               : %.8f km\n',dmax);
fprintf('Largeur de couronne                : %.8f km\n',boundary_width);
fprintf('------------------------------------------------------------\n');
fprintf('E[D] empirique                     : %.8f km\n',mean_distance_emp);
fprintf('E[D] theorie uniforme locale       : %.8f km\n',mean_distance_theory);
fprintf('Rapport moyenne emp/theorie        : %.8f\n', ...
    mean_distance_emp/mean_distance_theory);
fprintf('E[D^2] empirique                   : %.8f km^2\n',second_moment_emp);
fprintf('E[D^2] theorie                     : %.8f km^2\n',second_moment_theory);
fprintf('------------------------------------------------------------\n');
fprintf('P(couronne) empirique globale      : %.10f\n', ...
    boundary_probability_emp_global);
fprintf('P(couronne) theorie exacte         : %.10f\n', ...
    boundary_probability_theory_exact);
fprintf('P(couronne) theorie lineaire       : %.10f\n', ...
    boundary_probability_theory_linear);
fprintf('f_D(dmax) empirique approx         : %.10e km^-1\n', ...
    boundary_density_emp_global);
fprintf('f_D(dmax) theorie                  : %.10e km^-1\n', ...
    boundary_density_theory);
fprintf('Rapport bord emp/theorie           : %.10f\n', ...
    boundary_density_emp_global/boundary_density_theory);
fprintf('------------------------------------------------------------\n');
fprintf('Distance KS sur CDF                : %.6e\n',ks_distance);
fprintf('RMSE PDF globale                   : %.6e km^-1\n', ...
    rmse_pdf_global);
fprintf('============================================================\n');

%% ============================================================
%  10. Sauvegarde
%% ============================================================

output_file = fullfile(script_dir, ...
    'verif_loi_distances_liens_results.mat');

save(output_file, ...
    'analysis_file','R','dmax', ...
    'n_time_samples','time_indices', ...
    'n_distance_bins','n_latitude_bins', ...
    'boundary_fraction','boundary_width','boundary_lower', ...
    'min_links_per_lat_bin', ...
    'distance_edges','distance_centers','distance_bin_widths', ...
    'latitude_edges','latitude_centers', ...
    'all_link_distances','all_link_latitudes', ...
    'global_counts','pdf_emp_global','pdf_uniform_local', ...
    'cdf_emp_global','cdf_uniform_local', ...
    'distance_hist_per_lat_bin','pdf_emp_by_latitude', ...
    'link_count_per_lat_bin','boundary_count_per_lat_bin', ...
    'boundary_probability_emp_global', ...
    'boundary_probability_theory_exact', ...
    'boundary_probability_theory_linear', ...
    'boundary_density_emp_global','boundary_density_theory', ...
    'boundary_probability_emp_by_latitude', ...
    'boundary_density_emp_by_latitude', ...
    'boundary_density_ratio_emp_theory', ...
    'mean_distance_emp','mean_distance_theory', ...
    'second_moment_emp','second_moment_theory', ...
    'ks_distance','rmse_pdf_global');

fprintf('Resultats sauvegardes dans %s\n',output_file);

%% ============================================================
%  Fonction locale
%% ============================================================

function ratio = safe_divide(num,den)
    ratio = nan(size(num));
    valid = den > 0;
    ratio(valid) = num(valid)./den(valid);
end
