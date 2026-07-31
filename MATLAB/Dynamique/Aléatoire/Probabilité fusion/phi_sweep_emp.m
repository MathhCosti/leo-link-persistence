%% phi_sweep_empi.m
% Estimation empirique de phi_sweep :
%
%   phi_sweep_empi = somme(A_new,sat) / somme(A_geom)
%
% avec :
%   A_geom = 2*dmax*ell
%   A_new,sat = |B_i(t+dt) \ B_i(t)|
%
% Le code utilise analysis_temp_results.mat contenant Positions, R et dmax.
%
% La zone de liaison est définie de manière cohérente par la distance
% corde ||x-center|| <= dmax, aussi bien pour l'échantillonnage que
% pour le test d'appartenance.

clear; clc; close all;

%% Paramètres Monte-Carlo
n_mc_per_satellite = 30000;
rng(1);
max_transitions = inf;

%% Chargement
script_dir = fileparts(mfilename('fullpath'));
candidate_files = {
    fullfile(script_dir, 'analysis_temp_results.mat')
    fullfile(script_dir, 'analysis_temp_results(4).mat')
    fullfile(script_dir, '..', 'analysis_temp_results.mat')
    fullfile(script_dir, '..', 'analysis_temp_results(4).mat')
};

input_file = '';
for i_file = 1:numel(candidate_files)
    if isfile(candidate_files{i_file})
        input_file = candidate_files{i_file};
        break;
    end
end

if isempty(input_file)
    error(['Fichier introuvable. Place analysis_temp_results.mat ', ...
           'dans le dossier du script ou dans son dossier parent.']);
end

S = load(input_file);
required_fields = {'Positions','R','dmax'};
for i_field = 1:numel(required_fields)
    if ~isfield(S,required_fields{i_field})
        error('Le fichier doit contenir la variable %s.',required_fields{i_field});
    end
end

Positions = S.Positions;
R = S.R;
dmax = S.dmax;

if ~iscell(Positions)
    error('Positions doit être une cellule temporelle.');
end

Nt = numel(Positions);
if Nt < 2
    error('Il faut au moins deux instants temporels.');
end

n_transitions = min(Nt-1,max_transitions);

%% Stockage
phi_sweep_t = nan(n_transitions,1);
area_geom_t = zeros(n_transitions,1);
area_new_sat_t = zeros(n_transitions,1);

phi_sweep_satellite_all = [];
ell_satellite_all = [];
area_geom_satellite_all = [];
area_new_satellite_all = [];

n_new_samples_t = zeros(n_transitions,1);
n_total_samples_t = zeros(n_transitions,1);

% Pour une calotte définie par une distance corde <= dmax,
% l'aire sphérique vaut exactement pi*dmax^2.
disk_area_cap = pi*dmax^2;

%% Boucle temporelle
for k = 1:n_transitions
    P_t = Positions{k};
    P_next = Positions{k+1};

    if size(P_t,2) ~= 3 || size(P_next,2) ~= 3
        error('Positions{%d} et Positions{%d} doivent être N x 3.',k,k+1);
    end
    if size(P_t,1) ~= size(P_next,1)
        error('Le nombre de satellites change à la transition %d.',k);
    end

    N = size(P_t,1);

    for i = 1:N
        old_center = P_t(i,:);
        new_center = P_next(i,:);

        % Déplacement géodésique du centre entre t et t+dt
        % Les centres restent sur la sphère de rayon R.
        cos_alpha_move = dot(old_center,new_center) / ...
            (norm(old_center)*norm(new_center));
        cos_alpha_move = min(max(cos_alpha_move,-1),1);
        alpha_move = acos(cos_alpha_move);
        ell_i = R*alpha_move;

        if ell_i <= 1e-12
            continue;
        end

        % Aire géométrique théorique
        A_geom_i = 2*dmax*ell_i;

        % Aire réellement nouvelle du satellite
        sample_points = sample_spherical_disk(new_center,R,dmax,n_mc_per_satellite);
        covered_by_old_disk = points_covered_by_center(sample_points,old_center,dmax);
        is_new_for_satellite = ~covered_by_old_disk;

        n_new_i = nnz(is_new_for_satellite);
        A_new_i = disk_area_cap * n_new_i / n_mc_per_satellite;

        area_geom_t(k) = area_geom_t(k) + A_geom_i;
        area_new_sat_t(k) = area_new_sat_t(k) + A_new_i;
        n_new_samples_t(k) = n_new_samples_t(k) + n_new_i;
        n_total_samples_t(k) = n_total_samples_t(k) + n_mc_per_satellite;

        phi_i = A_new_i / A_geom_i;
        phi_sweep_satellite_all(end+1,1) = phi_i; %#ok<SAGROW>
        ell_satellite_all(end+1,1) = ell_i; %#ok<SAGROW>
        area_geom_satellite_all(end+1,1) = A_geom_i; %#ok<SAGROW>
        area_new_satellite_all(end+1,1) = A_new_i; %#ok<SAGROW>
    end

    if area_geom_t(k) > 0
        phi_sweep_t(k) = area_new_sat_t(k)/area_geom_t(k);
    end

    if mod(k,max(1,floor(n_transitions/10))) == 0 || k == n_transitions
        fprintf('Transition %d / %d terminée\n',k,n_transitions);
    end
end

%% Agrégation
total_area_geom = sum(area_geom_t);
total_area_new_sat = sum(area_new_sat_t);

if total_area_geom > 0
    phi_sweep_empi = total_area_new_sat / total_area_geom;
else
    phi_sweep_empi = NaN;
    warning('Aucune aire géométrique balayée non nulle.');
end

valid_t = ~isnan(phi_sweep_t);
phi_sweep_mean_t = mean(phi_sweep_t(valid_t));
phi_sweep_std_t = std(phi_sweep_t(valid_t));
phi_sweep_mean_satellite = mean(phi_sweep_satellite_all);
phi_sweep_std_satellite = std(phi_sweep_satellite_all);

%% Comparaison à la formule plane exacte
ell_clamped = min(max(ell_satellite_all,0),2*dmax);
A_inter_exact = 2*dmax^2 .* acos(ell_clamped./(2*dmax)) ...
    - 0.5 .* ell_clamped .* sqrt(max(4*dmax^2-ell_clamped.^2,0));
A_new_exact_all = pi*dmax^2 - A_inter_exact;
phi_exact_all = A_new_exact_all ./ max(2*dmax.*ell_clamped,eps);
phi_exact_all(ell_clamped <= 1e-12) = 1;
phi_sweep_exact_ratio = sum(A_new_exact_all) / sum(2*dmax.*ell_clamped);

%% Affichage
fprintf('\n');
fprintf('=================================================================\n');
fprintf(' ESTIMATION EMPIRIQUE DE phi_sweep\n');
fprintf('=================================================================\n');
fprintf('Fichier chargé                          : %s\n',input_file);
fprintf('Transitions utilisées                  : %d\n',n_transitions);
fprintf('Points MC par satellite/transition     : %d\n',n_mc_per_satellite);
fprintf('-----------------------------------------------------------------\n');
fprintf('Aire géométrique totale                : %.8e km^2\n',total_area_geom);
fprintf('Aire nouvelle satellite totale         : %.8e km^2\n',total_area_new_sat);
fprintf('-----------------------------------------------------------------\n');
fprintf('phi_sweep_empi                           : %.6f\n',phi_sweep_empi);
fprintf('Moyenne temporelle                     : %.6f\n',phi_sweep_mean_t);
fprintf('Écart-type temporel                    : %.6f\n',phi_sweep_std_t);
fprintf('Moyenne simple par satellite           : %.6f\n',phi_sweep_mean_satellite);
fprintf('Écart-type par satellite               : %.6f\n',phi_sweep_std_satellite);
fprintf('phi_sweep plan exact                   : %.6f\n',phi_sweep_exact_ratio);
fprintf('Écart MC - plan exact                  : %.6e\n', ...
    phi_sweep_empi - phi_sweep_exact_ratio);
fprintf('=================================================================\n\n');

%% Axe temporel
if isfield(S,'time_values') && numel(S.time_values) >= n_transitions+1
    x = S.time_values(1:n_transitions);
    x = x(:);
    x_label = 'Temps avant transition (s)';
else
    x = (1:n_transitions).';
    x_label = 'Indice de transition';
end

%% Figures
figure;
plot(x,phi_sweep_t,'LineWidth',1.4);
grid on;
xlabel(x_label);
ylabel('\phi_{sweep}(t)');
title('Facteur géométrique temporel \phi_{sweep}');
yline(phi_sweep_empi,'--',sprintf('Agrégée = %.4f',phi_sweep_empi), ...
    'LabelHorizontalAlignment','left');

figure;
histogram(phi_sweep_satellite_all);
grid on;
xlabel('\phi_{sweep,i}');
ylabel('Nombre de satellites-transition');
title('Distribution de \phi_{sweep} par satellite');

figure;
plot(x,area_geom_t,'LineWidth',1.2);
hold on;
plot(x,area_new_sat_t,'LineWidth',1.2);
grid on;
xlabel(x_label);
ylabel('Aire [km^2]');
title('Aire géométrique et aire réellement nouvelle');
legend('Aire géométrique','Aire nouvelle satellite','Location','best');

%% Sauvegarde
output_file = fullfile(script_dir,'phi_sweep_empi_results.mat');
save(output_file, ...
    'phi_sweep_empi','phi_sweep_mean_t','phi_sweep_std_t', ...
    'phi_sweep_mean_satellite','phi_sweep_std_satellite', ...
    'phi_sweep_exact_ratio','phi_sweep_t','phi_sweep_satellite_all', ...
    'phi_exact_all','area_geom_t','area_new_sat_t', ...
    'area_geom_satellite_all','area_new_satellite_all', ...
    'ell_satellite_all','total_area_geom','total_area_new_sat', ...
    'n_new_samples_t','n_total_samples_t','n_mc_per_satellite', ...
    'n_transitions','R','dmax','x');

fprintf('Résultats sauvegardés dans %s\n',output_file);

%% Fonctions locales
function points = sample_spherical_disk(center,R,dmax,n_points)
    % Échantillonnage uniforme en surface dans la calotte sphérique
    % définie par la condition de distance corde
    %
    %     ||x-center|| <= dmax.
    %
    % L'angle géocentrique maximal vérifie
    %
    %     dmax = 2R sin(alpha_max/2).

    e_r = center / norm(center);

    ref = [0 0 1];
    if abs(dot(e_r,ref)) > 0.95
        ref = [1 0 0];
    end

    e1 = cross(e_r,ref);
    e1 = e1 / norm(e1);

    e2 = cross(e_r,e1);
    e2 = e2 / norm(e2);

    alpha_max = 2*asin(min(dmax/(2*R),1));

    % Uniformité surfacique sur la calotte :
    % cos(alpha) uniforme sur [cos(alpha_max),1].
    u = rand(n_points,1);
    cos_alpha = 1 - u*(1-cos(alpha_max));
    sin_alpha = sqrt(max(1-cos_alpha.^2,0));

    theta = 2*pi*rand(n_points,1);

    tangent_direction = ...
        cos(theta).*e1 + sin(theta).*e2;

    points = ...
        R*cos_alpha.*e_r + ...
        R*sin_alpha.*tangent_direction;
end

function covered = points_covered_by_center(sample_points,center,dmax)
    differences = sample_points-center;
    D2 = sum(differences.^2,2);
    covered = D2 <= dmax^2;
end
