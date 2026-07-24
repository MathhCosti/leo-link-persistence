%% vitesse_radiale_depuis_analyse_delta.m
% Calcule la vitesse radiale inter-satellites à partir du fichier
% leo_zigzag_analysis_results_delta.mat.
%
% Pour un lien (i,j), la vitesse radiale exacte est
%
%   v_rad,ij = ((r_j-r_i) . (v_j-v_i)) / ||r_j-r_i||.
%
% v_rad > 0 : les satellites s'éloignent ;
% v_rad < 0 : les satellites se rapprochent.
%
% Le script calcule :
% - la vitesse radiale signée de tous les liens existants ;
% - la partie sortante max(v_rad,0) pertinente pour les ruptures ;
% - la vitesse radiale des liens proches de dmax ;
% - des moyennes temporelles et globales ;
% - une estimation empirique de q_break par flux de frontière.
%
% Les vitesses individuelles sont reconstruites par différences finies
% à partir des positions sauvegardées.

clearvars;
clc;
close all;

%% ========================= FICHIER ==================================

script_dir = fileparts(mfilename('fullpath'));
data_file = fullfile(script_dir, 'leo_zigzag_analysis_results_delta.mat');

if ~isfile(data_file)
    error('Fichier introuvable : %s', data_file);
end

S = load(data_file, ...
    'Positions', 'Adjacency', 'time_values', ...
    'N', 'R', 'dmax', 'dt', 'inc_deg');

required = {'Positions','Adjacency','time_values','N','R','dmax','dt'};
for k = 1:numel(required)
    if ~isfield(S, required{k})
        error('Variable manquante dans le .mat : %s', required{k});
    end
end

Positions   = S.Positions;
Adjacency   = S.Adjacency;
time_values = S.time_values(:);
N           = double(S.N);
R           = double(S.R);
dmax        = double(S.dmax);
dt          = double(S.dt);

Nt = numel(time_values);

if numel(Positions) ~= Nt || numel(Adjacency) ~= Nt
    error('Positions, Adjacency et time_values ont des tailles incompatibles.');
end

%% ====================== PARAMETRES D'ANALYSE ========================

% Largeur de la couche proche de la frontière de communication.
% Les liens de distance >= dmax - boundary_width sont considérés
% comme des liens proches de la rupture.
boundary_width = 100; % km

% Tolérance numérique sur les distances.
distance_tolerance = 1e-9;

%% ================= RECONSTRUCTION DES VITESSES ======================

Velocities = cell(Nt,1);

for k = 1:Nt
    Xk = double(Positions{k});

    if size(Xk,1) ~= N || size(Xk,2) ~= 3
        error('Positions{%d} doit être une matrice N x 3.', k);
    end

    if k == 1
        % Différence avant
        dt_local = time_values(2) - time_values(1);
        Velocities{k} = ...
            (double(Positions{2}) - double(Positions{1})) / dt_local;

    elseif k == Nt
        % Différence arrière
        dt_local = time_values(Nt) - time_values(Nt-1);
        Velocities{k} = ...
            (double(Positions{Nt}) - double(Positions{Nt-1})) / dt_local;

    else
        % Différence centrée
        dt_local = time_values(k+1) - time_values(k-1);
        Velocities{k} = ...
            (double(Positions{k+1}) - double(Positions{k-1})) / dt_local;
    end
end

%% ===================== CALCUL PAR LIEN ==============================

mean_vrad_signed_t   = NaN(Nt,1);
mean_vrad_abs_t      = NaN(Nt,1);
mean_vrad_out_t      = NaN(Nt,1);
mean_vrad_in_t       = NaN(Nt,1);

mean_vrad_out_boundary_t = NaN(Nt,1);
mean_vrad_all_boundary_t = NaN(Nt,1);

n_links_t            = zeros(Nt,1);
n_outgoing_links_t   = zeros(Nt,1);
n_boundary_links_t   = zeros(Nt,1);
n_out_boundary_t     = zeros(Nt,1);

all_vrad_signed      = [];
all_vrad_out         = [];
all_vrad_boundary    = [];
all_vrad_out_boundary = [];
all_distances        = [];

for k = 1:Nt

    X = double(Positions{k});
    V = double(Velocities{k});
    A = logical(Adjacency{k});

    [ii,jj] = find(triu(A,1));

    if isempty(ii)
        continue;
    end

    dr = X(jj,:) - X(ii,:);
    dv = V(jj,:) - V(ii,:);

    distances = sqrt(sum(dr.^2,2));

    valid = distances > distance_tolerance;
    dr = dr(valid,:);
    dv = dv(valid,:);
    distances = distances(valid);

    vrad = sum(dr .* dv,2) ./ distances;
    vrad_out = max(vrad,0);
    vrad_in = max(-vrad,0);

    boundary_mask = distances >= (dmax - boundary_width);
    outgoing_mask = vrad > 0;
    outgoing_boundary_mask = boundary_mask & outgoing_mask;

    n_links_t(k)          = numel(vrad);
    n_outgoing_links_t(k) = sum(outgoing_mask);
    n_boundary_links_t(k) = sum(boundary_mask);
    n_out_boundary_t(k)   = sum(outgoing_boundary_mask);

    mean_vrad_signed_t(k) = mean(vrad);
    mean_vrad_abs_t(k)    = mean(abs(vrad));
    mean_vrad_out_t(k)    = mean(vrad_out);
    mean_vrad_in_t(k)     = mean(vrad_in);

    if any(boundary_mask)
        mean_vrad_all_boundary_t(k) = mean(vrad(boundary_mask));
    end

    if any(outgoing_boundary_mask)
        mean_vrad_out_boundary_t(k) = ...
            mean(vrad(outgoing_boundary_mask));
    end

    all_vrad_signed = [all_vrad_signed; vrad]; %#ok<AGROW>
    all_vrad_out    = [all_vrad_out; vrad_out]; %#ok<AGROW>
    all_distances   = [all_distances; distances]; %#ok<AGROW>

    if any(boundary_mask)
        all_vrad_boundary = ...
            [all_vrad_boundary; vrad(boundary_mask)]; %#ok<AGROW>
    end

    if any(outgoing_boundary_mask)
        all_vrad_out_boundary = ...
            [all_vrad_out_boundary; ...
             vrad(outgoing_boundary_mask)]; %#ok<AGROW>
    end
end

%% ====================== MOYENNES GLOBALES ===========================

mean_vrad_signed = mean(all_vrad_signed,'omitnan');
mean_vrad_abs    = mean(abs(all_vrad_signed),'omitnan');

% Moyenne de la partie positive sur tous les liens :
% E[(dot D)_+ | lien]
mean_vrad_out = mean(all_vrad_out,'omitnan');

% Moyenne conditionnelle sachant que le lien s'éloigne :
% E[dot D | dot D > 0, lien]
positive_values = all_vrad_signed(all_vrad_signed > 0);
mean_vrad_given_outgoing = mean(positive_values,'omitnan');

% Moyenne sortante conditionnée par proximité de dmax :
% E[dot D | dot D > 0, D >= dmax-boundary_width]
mean_vrad_out_boundary = ...
    mean(all_vrad_out_boundary,'omitnan');

fraction_outgoing = ...
    numel(positive_values) / max(numel(all_vrad_signed),1);

fraction_boundary = ...
    numel(all_vrad_boundary) / max(numel(all_vrad_signed),1);

%% ========== ESTIMATION DE q_break PAR FLUX DE FRONTIERE =============

% Approximation plane classique :
%
% q_break ~= 2*Delta_t/dmax * E[(dot D)_+ | lien].
%
% Cette expression utilise la partie positive moyenne sur tous les liens.
q_break_flux_all_links = ...
    2 * dt / dmax * mean_vrad_out;

% Variante utilisant uniquement les liens proches de la frontière.
% La densité au bord n'est alors plus remplacée par 2/dmax :
% on estime directement la fraction de liens qui franchiraient dmax
% pendant un pas sous évolution linéaire.
q_break_crossing_estimate = NaN;

if ~isempty(all_vrad_signed)
    predicted_distances = all_distances + all_vrad_signed * dt;
    q_break_crossing_estimate = ...
        mean((all_distances <= dmax) & ...
             (predicted_distances > dmax));
end

%% ============================ AFFICHAGE =============================

fprintf('\n');
fprintf('====================================================================\n');
fprintf(' Vitesse radiale depuis analyse temporelle Walker Delta\n');
fprintf('====================================================================\n');
fprintf('N                                      : %d\n', N);
fprintf('Nombre de temps                        : %d\n', Nt);
fprintf('dt                                     : %.6f s\n', dt);
fprintf('R                                      : %.6f km\n', R);
fprintf('dmax                                   : %.6f km\n', dmax);

if isfield(S,'inc_deg')
    fprintf('Inclinaison                            : %.6f deg\n', ...
        double(S.inc_deg));
end

fprintf('Liens observés au total                : %d\n', ...
    numel(all_vrad_signed));
fprintf('--------------------------------------------------------------------\n');
fprintf('E[dot D | lien]                        : %.6f km/s\n', ...
    mean_vrad_signed);
fprintf('E[|dot D| | lien]                      : %.6f km/s\n', ...
    mean_vrad_abs);
fprintf('E[(dot D)_+ | lien]                    : %.6f km/s\n', ...
    mean_vrad_out);
fprintf('E[dot D | dot D>0, lien]               : %.6f km/s\n', ...
    mean_vrad_given_outgoing);
fprintf('Fraction de liens sortants             : %.6f\n', ...
    fraction_outgoing);
fprintf('--------------------------------------------------------------------\n');
fprintf('Largeur couche frontière               : %.3f km\n', ...
    boundary_width);
fprintf('Fraction de liens dans cette couche    : %.6f\n', ...
    fraction_boundary);
fprintf('E[dot D | sortant, proche de dmax]      : %.6f km/s\n', ...
    mean_vrad_out_boundary);
fprintf('--------------------------------------------------------------------\n');
fprintf('q_break flux, tous liens               : %.8f\n', ...
    q_break_flux_all_links);
fprintf('q_break franchissement direct          : %.8f\n', ...
    q_break_crossing_estimate);
fprintf('====================================================================\n\n');

%% ============================ GRAPHES ===============================

figure;
plot(time_values, mean_vrad_out_t, 'LineWidth', 1.5);
grid on;
xlabel('Temps (s)');
ylabel('E[(\dot D)_+ | lien] (km/s)');
title('Vitesse radiale sortante moyenne des liens');

figure;
plot(time_values, mean_vrad_out_boundary_t, 'LineWidth', 1.5);
grid on;
xlabel('Temps (s)');
ylabel('Vitesse radiale sortante (km/s)');
title(sprintf('Liens proches de d_{max} (couche %.0f km)', ...
    boundary_width));

figure;
histogram(all_vrad_signed,50,'Normalization','probability');
grid on;
xlabel('\dot D (km/s)');
ylabel('Fréquence');
title('Distribution des vitesses radiales signées');

figure;
scatter(all_distances, all_vrad_signed, 6, 'filled');
grid on;
xlabel('Distance du lien (km)');
ylabel('\dot D (km/s)');
title('Vitesse radiale en fonction de la longueur du lien');

%% ============================ SAUVEGARDE =============================

results = struct();

results.N = N;
results.R = R;
results.dmax = dmax;
results.dt = dt;
results.time_values = time_values;
results.boundary_width = boundary_width;

results.mean_vrad_signed_t = mean_vrad_signed_t;
results.mean_vrad_abs_t = mean_vrad_abs_t;
results.mean_vrad_out_t = mean_vrad_out_t;
results.mean_vrad_in_t = mean_vrad_in_t;
results.mean_vrad_out_boundary_t = mean_vrad_out_boundary_t;

results.mean_vrad_signed = mean_vrad_signed;
results.mean_vrad_abs = mean_vrad_abs;
results.mean_vrad_out = mean_vrad_out;
results.mean_vrad_given_outgoing = mean_vrad_given_outgoing;
results.mean_vrad_out_boundary = mean_vrad_out_boundary;

results.fraction_outgoing = fraction_outgoing;
results.fraction_boundary = fraction_boundary;

results.q_break_flux_all_links = q_break_flux_all_links;
results.q_break_crossing_estimate = q_break_crossing_estimate;

results.all_vrad_signed = all_vrad_signed;
results.all_distances = all_distances;

save('vitesse_radiale_delta_results.mat','results','-v7.3');

fprintf('Résultats sauvegardés dans vitesse_radiale_delta_results.mat\n');
