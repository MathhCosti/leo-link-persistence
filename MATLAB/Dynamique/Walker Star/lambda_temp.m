%% lambda_effective_temp_liens_inter.m
% Calcule une densite effective locale de satellites lambda_eff(t)
% a partir des donnees contenues dans liens_inter.mat.
%
% Idee : on garde la densite de satellites par strate
%     lambda_m(t) = N_m(t) / A_m
% puis on la moyenne avec les poids donnes par le nombre theorique de liens
% par strate L_m^th(t) issu de liens_inter.mat :
%
%     lambda_eff(t) = sum_m L_m^th(t) lambda_m(t) / sum_m L_m^th(t).
%
% Cette quantite peut ensuite etre injectee dans
%     p_merge(t) = 1 - exp(-2 lambda_eff(t) dmax v_rel(t) dt).

clear; clc; close all;

%% Chargement des donnees
script_dir = fileparts(mfilename('fullpath'));
mat_file = fullfile(script_dir, 'Nombre liens', 'liens_inter.mat');

if ~isfile(mat_file)
    error('Fichier %s introuvable. Place-le dans le dossier courant MATLAB.', mat_file);
end

S = load(mat_file);

required_vars = {'time_values','R','dmax','mean_links_incident_strate_t','mean_n_sat_strate_t'};
for q = 1:numel(required_vars)
    if ~isfield(S, required_vars{q})
        error('Variable manquante dans %s : %s', mat_file, required_vars{q});
    end
end

time_values = double(S.time_values(:));
R = double(S.R);
dmax = double(S.dmax);

% Nombre theorique de liens par strate en fonction du temps
L_th_t = double(S.mean_links_incident_strate_t);

% Nombre moyen de satellites par strate en fonction du temps
N_strate_t = double(S.mean_n_sat_strate_t);

% Remise au format [Nt x M]
Nt = numel(time_values);
if size(L_th_t,1) ~= Nt && size(L_th_t,2) == Nt
    L_th_t = L_th_t.';
end
if size(N_strate_t,1) ~= Nt && size(N_strate_t,2) == Nt
    N_strate_t = N_strate_t.';
end

if size(L_th_t,1) ~= Nt
    error('Dimension incompatible : mean_links_incident_strate_t doit avoir Nt lignes.');
end
if size(N_strate_t,1) ~= Nt
    error('Dimension incompatible : mean_n_sat_strate_t doit avoir Nt lignes.');
end
if size(L_th_t,2) ~= size(N_strate_t,2)
    error('Les nombres de strates de L_th_t et N_strate_t ne coincident pas.');
end

M = size(L_th_t,2);

%% Reconstruction des strates en latitude absolue
% Les strates couvrent la demi-sphere en latitude absolue ell in [0, pi/2].
% Comme elles representent les deux hemispheres repliees, l'aire de la strate
% [ell_a, ell_b] est :
%     A_m = 4*pi*R^2*(sin(ell_b)-sin(ell_a)).

if isfield(S, 'beta_step_strates')
    beta_step = double(S.beta_step_strates);
else
    beta_step = (pi/2)/M;
end

ell_edges = (0:M) * beta_step;
ell_edges = min(ell_edges, pi/2);
ell_edges(end) = pi/2;

% Si le dernier pas a cree deux bornes identiques, on repasse en decoupage uniforme.
if any(diff(ell_edges) <= 1e-12)
    warning('Strates degeneratees detectees. Utilisation d''un decoupage uniforme de [0, pi/2].');
    ell_edges = linspace(0, pi/2, M+1);
end

ell_centers = 0.5 * (ell_edges(1:end-1) + ell_edges(2:end));
ell_centers_deg = rad2deg(ell_centers);

A_strates = 4*pi*R^2 * (sin(ell_edges(2:end)) - sin(ell_edges(1:end-1)));
A_strates = reshape(A_strates, 1, []);

%% Correction d'orientation de N_strate_t par rapport aux latitudes
% Convention choisie dans ce script : colonnes de 1 a M = equateur -> pole.
% On verifie l'orientation de mean_n_sat_strate_t en regardant la densite
% moyenne aux extremites. Dans ce modele, la densite moyenne est maximale
% dans les strates polaires. Si la premiere strate est beaucoup plus dense
% que la derniere, on retourne N_strate_t.
lambda_tmp_no_flip = N_strate_t ./ A_strates;
lambda_mean_no_flip = mean(lambda_tmp_no_flip, 1, 'omitnan');

flip_N_strate = lambda_mean_no_flip(1) > lambda_mean_no_flip(end);
if flip_N_strate
    N_strate_t = fliplr(N_strate_t);
    fprintf('Orientation des strates : N_strate_t retourne pour suivre equateur -> pole.\n');
else
    fprintf('Orientation des strates : N_strate_t conserve tel quel equateur -> pole.\n');
end

%% Densite de satellites par strate corrigee
lambda_strate_t = N_strate_t ./ A_strates;  % [satellites / km^2]

%% Correction d'orientation des liens theoriques
% Convention choisie : colonnes de 1 a M = equateur -> pole.
% mean_links_incident_strate_t peut etre indexe dans le sens inverse selon le
% script generateur. On choisit l'orientation qui donne la meilleure coherence
% entre la densite moyenne de satellites corrigee et les poids moyens de liens.
lambda_strate_mean_tmp = mean(lambda_strate_t, 1, 'omitnan');
L_mean_no_flip = mean(L_th_t, 1, 'omitnan');
L_mean_flip    = fliplr(L_mean_no_flip);

c_no = corr(lambda_strate_mean_tmp(:), L_mean_no_flip(:), 'Rows','complete');
c_fl = corr(lambda_strate_mean_tmp(:), L_mean_flip(:),    'Rows','complete');

if isnan(c_no), c_no = -Inf; end
if isnan(c_fl), c_fl = -Inf; end

flip_L_th = c_fl > c_no;
if flip_L_th
    L_th_t = fliplr(L_th_t);
    fprintf('Orientation des strates : L_th_t retourne pour correspondre a N_strate_t.\n');
else
    fprintf('Orientation des strates : L_th_t conserve tel quel.\n');
end

%% Densite effective ponderee par les liens theoriques
sum_L_th = sum(L_th_t, 2);
lambda_eff_t = NaN(Nt,1);

valid = sum_L_th > 0;
lambda_eff_t(valid) = sum(L_th_t(valid,:) .* lambda_strate_t(valid,:), 2) ./ sum_L_th(valid);

%% References globales
surface_sphere = 4*pi*R^2;

if isfield(S, 'N_mean_theory')
    N_ref = double(S.N_mean_theory);
elseif isfield(S, 'N_all')
    N_ref = mean(double(S.N_all(:)), 'omitnan');
else
    N_ref = mean(sum(N_strate_t,2), 'omitnan');
end

lambda_global = N_ref / surface_sphere;

if isfield(S, 'lambda')
    lambda_file = double(S.lambda);
else
    lambda_file = lambda_global;
end

%% Affichage console
fprintf('\n=== Densite effective lambda depuis liens_inter.mat ===\n');
fprintf('Nombre de strates M                         : %d\n', M);
fprintf('Rayon orbital R                             : %.2f km\n', R);
fprintf('dmax                                        : %.2f km\n', dmax);
fprintf('N reference                                 : %.2f\n', N_ref);
fprintf('lambda globale N/(4*pi*R^2)                 : %.6e sat/km^2\n', lambda_global);
fprintf('lambda stockee dans liens_inter.mat          : %.6e sat/km^2\n', lambda_file);
fprintf('lambda_eff moyen temporel                   : %.6e sat/km^2\n', mean(lambda_eff_t,'omitnan'));
fprintf('lambda_eff min / max                        : %.6e / %.6e sat/km^2\n', ...
    min(lambda_eff_t,[],'omitnan'), max(lambda_eff_t,[],'omitnan'));

%% Figure 1 : evolution temporelle de lambda_eff(t)
figure;
plot(time_values, lambda_eff_t, 'LineWidth', 1.8); hold on;
yline(lambda_global, '--', sprintf('lambda globale = %.2e', lambda_global), ...
    'LabelHorizontalAlignment','left');
yline(mean(lambda_eff_t,'omitnan'), ':', sprintf('moyenne lambda_{eff} = %.2e', mean(lambda_eff_t,'omitnan')), ...
    'LabelHorizontalAlignment','left');
grid on;
xlabel('Temps (s)');
ylabel('\lambda_{sat,eff}(t) (sat/km^2)');
title('\lambda_{sat,eff}(t) ponderee par les liens theoriques par strate');
legend('\lambda_{sat,eff}(t)', '\lambda globale', 'Moyenne temporelle', 'Location','best');

%% Figure 2 : densites par strate moyennees temporellement
lambda_strate_mean = mean(lambda_strate_t, 1, 'omitnan');
L_th_mean = mean(L_th_t, 1, 'omitnan');

figure;
yyaxis left;
bar(ell_centers_deg, lambda_strate_mean, 0.8);
ylabel('\lambda_m moyen (sat/km^2)');

yyaxis right;
plot(ell_centers_deg, L_th_mean, '-o', 'LineWidth', 1.5);
ylabel('L_m^{th} moyen');

grid on;
xlabel('Latitude absolue de la strate (degres)');
title('Densite satellitaire moyenne et poids de liens par strate');
legend('\lambda_m moyen', 'L_m^{th} moyen', 'Location','best');

%% Sauvegarde
save('lambda_effective_temp_liens_inter_results.mat', ...
    'time_values', 'lambda_eff_t', 'lambda_global', 'lambda_file', ...
    'lambda_strate_t', 'lambda_strate_mean', 'N_strate_t', ...
    'L_th_t', 'L_th_mean', 'A_strates', 'ell_edges', 'ell_centers', ...
    'ell_centers_deg', 'R', 'dmax', 'N_ref', 'flip_L_th', 'flip_N_strate');

fprintf('\nResultats sauvegardes dans lambda_effective_temp_liens_inter_results.mat\n');
