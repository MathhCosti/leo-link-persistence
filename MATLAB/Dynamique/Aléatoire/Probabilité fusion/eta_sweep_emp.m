%% eta_sweep_empi_corrige.m
% Estimation empirique de la fraction d'aire nouvellement explorée
% par un satellite appartenant à une composante :
%
% eta_sweep_empi =
%   aire du nouveau disque non couverte par la composante à t
%   ---------------------------------------------------------
%   aire du nouveau disque non couverte par le satellite i à t
%
% Pour chaque satellite i et chaque transition t -> t+dt, on échantillonne
% uniformément dans sa zone de liaison à l'instant t+dt.
%
% - "aire brute nouvellement balayée" :
%       points qui ne sont plus dans la zone du même satellite à t ;
%
% - "aire réellement nouvelle pour la composante" :
%       points qui ne sont dans la zone d'aucun satellite de la
%       composante de i à l'instant t.
%
% La valeur principale est un rapport agrégé :
%
% eta_sweep_empi = somme(A_nouvelle_composante) / somme(A_nouvelle_satellite)

clear; clc; close all;

%% ============================================================
%  1. Paramètres Monte-Carlo
%% ============================================================

n_mc_per_satellite = 2000;
rng(1);

% Mettre inf pour utiliser toutes les transitions.
max_transitions = inf;

%% ============================================================
%  2. Chargement
%% ============================================================

script_dir = fileparts(mfilename('fullpath'));
input_file = fullfile(script_dir, '..', 'analysis_temp_results.mat');

if isempty(input_file)
    error(['Fichier introuvable. Place analysis_temp_results.mat ', ...
           'dans le dossier du script ou son dossier parent.']);
end

S = load(input_file);

required_fields = {'Adjacency', 'Positions', 'R', 'dmax'};
for i_field = 1:numel(required_fields)
    if ~isfield(S, required_fields{i_field})
        error('Le fichier doit contenir la variable %s.', ...
              required_fields{i_field});
    end
end

Adjacency = S.Adjacency;
Positions = S.Positions;
R = S.R;
dmax = S.dmax;

if ~iscell(Adjacency) || ~iscell(Positions)
    error('Adjacency et Positions doivent être des cellules temporelles.');
end

Nt = min(numel(Adjacency), numel(Positions));
if Nt < 2
    error('Il faut au moins deux instants temporels.');
end

n_transitions = min(Nt - 1, max_transitions);

%% ============================================================
%  3. Stockage
%% ============================================================

eta_sweep_t = nan(n_transitions,1);

area_raw_new_t       = zeros(n_transitions,1);
area_component_new_t = zeros(n_transitions,1);

n_raw_new_samples_t       = zeros(n_transitions,1);
n_component_new_samples_t = zeros(n_transitions,1);
n_total_samples_t         = zeros(n_transitions,1);

eta_sweep_satellite_all = [];

disk_area_planar = pi * dmax^2;

%% ============================================================
%  4. Boucle temporelle
%% ============================================================

for k = 1:n_transitions

    A_t = normalize_adjacency(Adjacency{k});

    P_t = Positions{k};
    P_next = Positions{k+1};

    if size(P_t,2) ~= 3 || size(P_next,2) ~= 3
        error('Positions{%d} et Positions{%d} doivent être N x 3.', ...
              k, k+1);
    end

    if size(P_t,1) ~= size(P_next,1) || size(P_t,1) ~= size(A_t,1)
        error('Tailles incompatibles à la transition %d.', k);
    end

    N = size(P_t,1);
    labels_t = conncomp(graph(A_t));

    for i = 1:N

        old_center = P_t(i,:);
        new_center = P_next(i,:);

        component_nodes = find(labels_t == labels_t(i));
        component_positions = P_t(component_nodes,:);

        %% Points uniformes dans le disque tangent au nouveau centre
        sample_points = sample_spherical_disk( ...
            new_center, R, dmax, n_mc_per_satellite);

        %% Aire brute nouvellement balayée par le satellite i
        covered_by_old_self = points_covered_by_component( ...
            sample_points, old_center, dmax);

        is_raw_new = ~covered_by_old_self;

        %% Aire réellement nouvelle pour toute la composante
        covered_by_old_component = points_covered_by_component( ...
            sample_points, component_positions, dmax);

        is_component_new = ~covered_by_old_component;

        n_raw_i = nnz(is_raw_new);
        n_component_i = nnz(is_component_new);

        % Par construction, toute aire nouvelle pour la composante
        % doit aussi être nouvelle pour le satellite lui-même.
        if any(is_component_new & ~is_raw_new)
            warning('Incohérence numérique détectée à t=%d, satellite=%d.', ...
                    k, i);
        end

        area_raw_i = disk_area_planar * n_raw_i / n_mc_per_satellite;
        area_component_i = ...
            disk_area_planar * n_component_i / n_mc_per_satellite;

        area_raw_new_t(k) = area_raw_new_t(k) + area_raw_i;
        area_component_new_t(k) = ...
            area_component_new_t(k) + area_component_i;

        n_raw_new_samples_t(k) = ...
            n_raw_new_samples_t(k) + n_raw_i;
        n_component_new_samples_t(k) = ...
            n_component_new_samples_t(k) + n_component_i;
        n_total_samples_t(k) = ...
            n_total_samples_t(k) + n_mc_per_satellite;

        if n_raw_i > 0
            eta_sweep_satellite_all(end+1,1) = ...
                n_component_i / n_raw_i; %#ok<SAGROW>
        end
    end

    if area_raw_new_t(k) > 0
        eta_sweep_t(k) = ...
            area_component_new_t(k) / area_raw_new_t(k);
    end

    if mod(k, max(1,floor(n_transitions/10))) == 0 ...
            || k == n_transitions
        fprintf('Transition %d / %d terminée\n', k, n_transitions);
    end
end

%% ============================================================
%  5. Agrégation
%% ============================================================

total_area_raw_new = sum(area_raw_new_t);
total_area_component_new = sum(area_component_new_t);

if total_area_raw_new > 0
    eta_sweep_empi = ...
        total_area_component_new / total_area_raw_new;
else
    eta_sweep_empi = NaN;
    warning(['Aucune aire brute nouvellement balayée n''a été détectée. ', ...
             'Augmente n_mc_per_satellite.']);
end

valid_t = ~isnan(eta_sweep_t);

if any(valid_t)
    eta_sweep_mean_t = mean(eta_sweep_t(valid_t));
    eta_sweep_std_t = std(eta_sweep_t(valid_t));
else
    eta_sweep_mean_t = NaN;
    eta_sweep_std_t = NaN;
end

if ~isempty(eta_sweep_satellite_all)
    eta_sweep_mean_satellite = mean(eta_sweep_satellite_all);
    eta_sweep_std_satellite = std(eta_sweep_satellite_all);
else
    eta_sweep_mean_satellite = NaN;
    eta_sweep_std_satellite = NaN;
end

%% ============================================================
%  6. Affichage
%% ============================================================

fprintf('\n');
fprintf('=================================================================\n');
fprintf(' ESTIMATION CORRIGEE DE eta_sweep_empi\n');
fprintf('=================================================================\n');
fprintf('Fichier chargé                              : %s\n', input_file);
fprintf('Transitions utilisées                      : %d\n', n_transitions);
fprintf('Points MC par satellite et transition      : %d\n', ...
        n_mc_per_satellite);
fprintf('-----------------------------------------------------------------\n');
fprintf('Aire brute nouvellement balayée             : %.8e km^2\n', ...
        total_area_raw_new);
fprintf('Aire nouvelle pour les composantes          : %.8e km^2\n', ...
        total_area_component_new);
fprintf('-----------------------------------------------------------------\n');
fprintf('eta_sweep_empi = Anew_comp / Anew_sat        : %.6f\n', ...
        eta_sweep_empi);
fprintf('Moyenne temporelle                          : %.6f\n', ...
        eta_sweep_mean_t);
fprintf('Écart-type temporel                         : %.6f\n', ...
        eta_sweep_std_t);
fprintf('Moyenne simple par satellite                : %.6f\n', ...
        eta_sweep_mean_satellite);
fprintf('Écart-type par satellite                    : %.6f\n', ...
        eta_sweep_std_satellite);
fprintf('-----------------------------------------------------------------\n');
fprintf('Échantillons dans aire brute nouvelle       : %d\n', ...
        sum(n_raw_new_samples_t));
fprintf('Échantillons dans aire nouvelle composante  : %d\n', ...
        sum(n_component_new_samples_t));
fprintf('=================================================================\n\n');

%% ============================================================
%  7. Axe temporel et figures
%% ============================================================

if isfield(S,'time_values') && numel(S.time_values) >= n_transitions+1
    x = S.time_values(1:n_transitions);
    x = x(:);
    x_label = 'Temps avant transition (s)';
else
    x = (1:n_transitions).';
    x_label = 'Indice de transition';
end

figure;
plot(x, eta_sweep_t, 'LineWidth', 1.4);
grid on;
xlabel(x_label);
ylabel('\eta_{sweep}(t)');
ylim([0 1]);
title('Fraction d''aire brute réellement nouvelle pour la composante');

if ~isnan(eta_sweep_empi)
    yline(eta_sweep_empi, '--', ...
        sprintf('Agrégée = %.4f', eta_sweep_empi), ...
        'LabelHorizontalAlignment', 'left');
end

figure;
plot(x, area_raw_new_t, 'LineWidth', 1.2);
hold on;
plot(x, area_component_new_t, 'LineWidth', 1.2);
grid on;
xlabel(x_label);
ylabel('Aire [km^2]');
title('Aires nouvelles par transition');
legend('Nouvelle pour le satellite', ...
       'Nouvelle pour la composante', ...
       'Location','best');

figure;
histogram(eta_sweep_satellite_all);
grid on;
xlabel('\eta_{sweep,i}');
ylabel('Nombre de satellites-transition');
title('Distribution de \eta_{sweep} par satellite');

%% ============================================================
%  8. Sauvegarde
%% ============================================================

output_file = fullfile(script_dir, 'eta_sweep_emp_results.mat');

save(output_file, ...
    'eta_sweep_empi', ...
    'eta_sweep_mean_t', ...
    'eta_sweep_std_t', ...
    'eta_sweep_mean_satellite', ...
    'eta_sweep_std_satellite', ...
    'eta_sweep_t', ...
    'eta_sweep_satellite_all', ...
    'area_raw_new_t', ...
    'area_component_new_t', ...
    'total_area_raw_new', ...
    'total_area_component_new', ...
    'n_raw_new_samples_t', ...
    'n_component_new_samples_t', ...
    'n_total_samples_t', ...
    'n_mc_per_satellite', ...
    'n_transitions', ...
    'R', 'dmax', 'x');

fprintf('Résultats sauvegardés dans %s\n', output_file);

%% ========================================================================
% Fonctions locales
%% ========================================================================

function A = normalize_adjacency(A_in)
    A = logical(spones(A_in));
    A = A | A.';
    A(1:size(A,1)+1:end) = false;
    A = sparse(A);
end

function points = sample_spherical_disk(center, R, radius, n_points)
    % Échantillonnage uniforme dans le disque du plan tangent, puis
    % projection exponentielle sur la sphère.

    e_r = center / norm(center);

    % Construction d'une base orthonormée tangentielle
    ref = [0 0 1];
    if abs(dot(e_r,ref)) > 0.95
        ref = [1 0 0];
    end

    e1 = cross(e_r,ref);
    e1 = e1 / norm(e1);
    e2 = cross(e_r,e1);
    e2 = e2 / norm(e2);

    rho = radius * sqrt(rand(n_points,1));
    theta = 2*pi*rand(n_points,1);

    tangent_vectors = ...
        (rho.*cos(theta)).*e1 + ...
        (rho.*sin(theta)).*e2;

    points = expmap_sphere(center,tangent_vectors,R);
end

function points = expmap_sphere(origin,tangent_vectors,R)

    n_points = size(tangent_vectors,1);
    points = zeros(n_points,3);

    e_r = origin / norm(origin);
    rho = vecnorm(tangent_vectors,2,2);

    zero_mask = rho <= 1e-14;
    nonzero_mask = ~zero_mask;

    if any(zero_mask)
        points(zero_mask,:) = repmat(origin,nnz(zero_mask),1);
    end

    if any(nonzero_mask)
        directions = tangent_vectors(nonzero_mask,:) ...
            ./ rho(nonzero_mask);

        angles = rho(nonzero_mask) / R;

        points(nonzero_mask,:) = ...
            R.*cos(angles).*e_r + ...
            R.*sin(angles).*directions;
    end
end

function covered = points_covered_by_component( ...
        sample_points,component_positions,dmax)

    if isvector(component_positions) && numel(component_positions)==3
        component_positions = reshape(component_positions,1,3);
    end

    n_points = size(sample_points,1);
    covered = false(n_points,1);

    block_size = 2000;

    for i_start = 1:block_size:n_points
        i_end = min(i_start+block_size-1,n_points);
        block = sample_points(i_start:i_end,:);

        dx = block(:,1) - component_positions(:,1).';
        dy = block(:,2) - component_positions(:,2).';
        dz = block(:,3) - component_positions(:,3).';

        D2 = dx.^2 + dy.^2 + dz.^2;

        covered(i_start:i_end) = any(D2 <= dmax^2,2);
    end
end
