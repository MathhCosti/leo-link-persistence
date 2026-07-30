%% lambda_effective_temp_delta_spatial_pond_satellites.m
% Densité satellitaire effective temporelle pour un Walker-Delta
% à uniformité spatiale initiale.
%
% La densité locale de la strate m est :
%
%     lambda_m(t) = N_m(t) / A_m
%
% avec, pour une strate de latitude absolue
%
%     phi_m <= |phi| < phi_{m+1},
%
% l'aire Nord + Sud :
%
%     A_m = 4*pi*R^2 [sin(phi_{m+1}) - sin(phi_m)].
%
% La densité effective est pondérée par le nombre de satellites
% présents dans chaque strate :
%
%     lambda_eff_sat(t)
%       = sum_m N_m(t) lambda_m(t) / sum_m N_m(t).
%
% Cette pondération évite de surpondérer les strates déjà très
% connectées, notamment près des latitudes extrêmes.
%
% Le script charge par défaut :
%
%     liens_inter_delta_spatial_min_max_analytique_empirique.mat
%
% Il accepte également un fichier contenant des liens incidents si cette
% variable est disponible, mais privilégie les liens internes afin d'éviter
% le double comptage des liens entre deux strates.

clear; clc; close all;

%% ============================================================
%  FICHIER D'ENTREE
%% ============================================================

script_dir = fileparts(mfilename('fullpath'));
liens_file = fullfile(script_dir, 'Nombre liens', 'liens_quadrature_results.mat');

if isempty(liens_file)
    error(['Aucun fichier de résultats Walker-Delta spatial trouvé. ', ...
        'Place dans le même dossier l''un des fichiers suivants :\n', ...
        '  liens_inter_delta_spatial_min_max_analytique_empirique.mat\n', ...
        '  liens_inter_delta_spatial_corrige.mat\n', ...
        '  liens_inter_delta_spatial.mat']);
end

fprintf('Fichier chargé : %s\n', liens_file);
S = load(liens_file);

%% ============================================================
%  VARIABLES OBLIGATOIRES
%% ============================================================

required_vars = {'time_values','R','dmax','inc'};

for q = 1:numel(required_vars)
    if ~isfield(S, required_vars{q})
        error('Variable manquante dans %s : %s', ...
            liens_file, required_vars{q});
    end
end

time_values = double(S.time_values(:));
R = double(S.R);
dmax = double(S.dmax);
inc = double(S.inc);
inc_deg = rad2deg(inc);

Nt = numel(time_values);

%% ============================================================
%  NOMBRE DE SATELLITES PAR STRATE
%% ============================================================

if isfield(S, 'mean_n_sat_strate_t')
    N_strate_t = double(S.mean_n_sat_strate_t);

elseif isfield(S, 'N_strate_theory_time')
    N_strate_t = double(S.N_strate_theory_time);

else
    error(['Aucune variable de nombre de satellites par strate trouvée. ', ...
        'Variable attendue : mean_n_sat_strate_t ou N_strate_theory_time.']);
end

% Format [Nt x M]
if size(N_strate_t,1) ~= Nt && size(N_strate_t,2) == Nt
    N_strate_t = N_strate_t.';
end

if size(N_strate_t,1) ~= Nt
    error('N_strate_t doit avoir Nt = %d lignes.', Nt);
end

M = size(N_strate_t,2);

%% ============================================================
%  POIDS SATELLITAIRES PAR STRATE
%% ============================================================

% La pondération utilise directement le nombre de satellites présents
% dans chaque strate :
%
%     W_m(t) = N_m(t).
%
% Cela évite le double renforcement suivant :
%
%     forte densité -> beaucoup de liens -> poids encore plus élevé.
%
% Cette correction est particulièrement importante lorsque inc = 90 deg,
% car les dernières strates proches des pôles ont une aire très faible.

N_weight_t = N_strate_t;
weight_type = 'nombre de satellites par strate';

%% ============================================================
%  RECONSTRUCTION DES BORDS DE STRATES DELTA
%% ============================================================

% Cas privilégié : la table de comparaison sauvegardée contient les
% véritables frontières de latitude.
if isfield(S, 'strates_comparison') && ...
        all(ismember({'phi_in','phi_out'}, ...
        S.strates_comparison.Properties.VariableNames))

    phi_in = double(S.strates_comparison.phi_in(:));
    phi_out = double(S.strates_comparison.phi_out(:));

elseif isfield(S, 'comparison_strates') && ...
        all(ismember({'latitude_inner_deg','latitude_outer_deg'}, ...
        S.comparison_strates.Properties.VariableNames))

    phi_in = deg2rad(double( ...
        S.comparison_strates.latitude_inner_deg(:)));
    phi_out = deg2rad(double( ...
        S.comparison_strates.latitude_outer_deg(:)));

elseif isfield(S, 'strates_delta') && ...
        isfield(S.strates_delta, 'active_table')

    T = S.strates_delta.active_table;

    if all(ismember({'latitude_inner','latitude_outer'}, ...
            T.Properties.VariableNames))
        phi_in = double(T.latitude_inner(:));
        phi_out = double(T.latitude_outer(:));
    else
        error('La table strates_delta ne contient pas les latitudes attendues.');
    end

else
    % Repli simple : découpage uniforme de [0,inc].
    warning(['Frontières exactes de strates absentes. ', ...
        'Utilisation d''un découpage uniforme de [0,inc].']);
    phi_edges = linspace(0,inc,M+1);
    phi_in = phi_edges(1:end-1).';
    phi_out = phi_edges(2:end).';
end

if numel(phi_in) ~= M || numel(phi_out) ~= M
    error(['Le nombre de frontières de strates (%d) ne correspond pas ', ...
        'au nombre de colonnes de N_strate_t (%d).'], ...
        numel(phi_in), M);
end

phi_mid = 0.5*(phi_in+phi_out);
phi_mid_deg = rad2deg(phi_mid);

%% ============================================================
%  ORIENTATION DES STRATES
%% ============================================================

% Convention imposée : colonnes équateur -> latitude extrême.
[phi_in, order] = sort(phi_in, 'ascend');
phi_out = phi_out(order);
phi_mid = phi_mid(order);
phi_mid_deg = phi_mid_deg(order);

N_strate_t = N_strate_t(:,order);
N_weight_t = N_weight_t(:,order);

%% ============================================================
%  AIRES ET DENSITES LOCALES
%% ============================================================

A_strates = 4*pi*R^2 .* ...
    (sin(phi_out)-sin(phi_in));

A_strates = reshape(A_strates,1,[]);

if any(A_strates <= 0)
    error('Une ou plusieurs strates possèdent une aire non positive.');
end

lambda_strate_t = N_strate_t ./ A_strates;

%% ============================================================
%  DENSITE EFFECTIVE PONDEREE PAR LES SATELLITES
%% ============================================================

sum_weights = sum(N_weight_t,2);
lambda_eff_t = NaN(Nt,1);

valid = sum_weights > 0;

lambda_eff_t(valid) = ...
    sum(N_weight_t(valid,:).*lambda_strate_t(valid,:),2) ...
    ./ sum_weights(valid);

% Sous la forme développée :
%
%     lambda_eff(t)
%       = sum_m N_m(t)^2/A_m
%         ------------------
%            sum_m N_m(t)
%
% Aucun repli fondé sur les liens n'est nécessaire, puisque les poids
% sont précisément les nombres de satellites.

%% ============================================================
%  REFERENCES GLOBALES
%% ============================================================

surface_band = 4*pi*R^2*sin(inc);

if isfield(S, 'N_mean')
    N_ref = double(S.N_mean);

elseif isfield(S, 'N_mean_theory')
    N_ref = double(S.N_mean_theory);

elseif isfield(S, 'N_all')
    N_ref = mean(double(S.N_all(:)), 'omitnan');

else
    N_ref = mean(sum(N_strate_t,2), 'omitnan');
end

lambda_global_band = N_ref/surface_band;

if isfield(S, 'lambda')
    lambda_file = double(S.lambda);
else
    lambda_file = lambda_global_band;
end

%% ============================================================
%  STATISTIQUES
%% ============================================================

lambda_strate_mean = mean(lambda_strate_t,1,'omitnan');
N_weight_mean = mean(N_weight_t,1,'omitnan');

lambda_eff_mean = mean(lambda_eff_t,'omitnan');
lambda_eff_min = min(lambda_eff_t,[],'omitnan');
lambda_eff_max = max(lambda_eff_t,[],'omitnan');

fprintf('\n=== Densité effective Walker-Delta spatial ===\n');
fprintf('Inclinaison                                 : %.2f deg\n', ...
    inc_deg);
fprintf('Nombre de strates                          : %d\n', M);
fprintf('Type de poids                              : %s\n', ...
    weight_type);
fprintf('Rayon orbital                              : %.2f km\n', R);
fprintf('dmax                                       : %.2f km\n', dmax);
fprintf('N de référence                             : %.3f\n', N_ref);
fprintf('Aire de la bande                           : %.6e km^2\n', ...
    surface_band);
fprintf('lambda globale dans la bande               : %.6e sat/km^2\n', ...
    lambda_global_band);
fprintf('lambda stockée dans le fichier             : %.6e sat/km^2\n', ...
    lambda_file);
fprintf('lambda_eff moyenne                         : %.6e sat/km^2\n', ...
    lambda_eff_mean);
fprintf('lambda_eff min / max                       : %.6e / %.6e\n', ...
    lambda_eff_min,lambda_eff_max);

%% ============================================================
%  FIGURE 1 : LAMBDA EFFECTIVE TEMPORELLE
%% ============================================================

figure;
hold on;
grid on;

plot(time_values,lambda_eff_t,'LineWidth',1.8, ...
    'DisplayName','\lambda_{eff,sat}(t)');

yline(lambda_global_band,'--', ...
    sprintf('\\lambda bande = %.2e',lambda_global_band), ...
    'LineWidth',1.2, ...
    'LabelHorizontalAlignment','left', ...
    'DisplayName','Densité globale dans la bande');

yline(lambda_eff_mean,':', ...
    sprintf('moyenne = %.2e',lambda_eff_mean), ...
    'LineWidth',1.2, ...
    'LabelHorizontalAlignment','left', ...
    'DisplayName','Moyenne temporelle');

xlabel('Temps (s)');
ylabel('\lambda_{eff,sat}^{\Delta,sp}(t) (sat/km^2)');
title('Densité effective pondérée par les satellites — Walker-Delta spatial');
legend('Location','best');
hold off;

%% ============================================================
%  FIGURE 2 : DENSITE ET POIDS PAR STRATE
%% ============================================================

figure;

yyaxis left;
bar(phi_mid_deg,lambda_strate_mean,0.8);
ylabel('\lambda_m moyen (sat/km^2)');

yyaxis right;
plot(phi_mid_deg,N_weight_mean,'-o','LineWidth',1.5);
ylabel('Nombre moyen de satellites');

grid on;
xlabel('Latitude absolue de la strate (degrés)');
title('Densité locale moyenne et poids satellitaire par strate');
legend('\lambda_m moyen',weight_type,'Location','best');

%% ============================================================
%  FIGURE 3 : CARTE TEMPORELLE DES DENSITES PAR STRATE
%% ============================================================

figure;
imagesc(time_values,phi_mid_deg,lambda_strate_t.');
axis xy;
colorbar;
xlabel('Temps (s)');
ylabel('Latitude absolue de la strate (degrés)');
title('\lambda_m(t) par strate — Walker-Delta spatial');

%% ============================================================
%  SAUVEGARDE
%% ============================================================

save('lambda_eff_th_results.mat', ...
    'time_values','lambda_eff_t','lambda_eff_mean', ...
    'lambda_eff_min','lambda_eff_max', ...
    'lambda_global_band','lambda_file', ...
    'lambda_strate_t','lambda_strate_mean', ...
    'N_strate_t','N_weight_t','N_weight_mean','weight_type', ...
    'A_strates','phi_in','phi_out','phi_mid','phi_mid_deg', ...
    'R','dmax','inc','inc_deg','surface_band','N_ref', ...
    'liens_file');

fprintf('\nRésultats sauvegardés dans ');
fprintf('lambda_effective_temp_delta_spatial_results.mat\n');
