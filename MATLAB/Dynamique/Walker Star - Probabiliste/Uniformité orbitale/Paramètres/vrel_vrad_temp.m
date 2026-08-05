%% vrel_vrad_emp_delta_latitude.m
% Diagnostic des relations entre vitesse orbitale, vitesse relative
% et vitesse radiale dans le modele Walker Delta a uniformite orbitale.
%
% Les grandeurs sont calculees sur les liens existants, puis moyennees
% temporellement par tranche de latitude.
%
% Pour un lien (i,j), la latitude retenue est celle du milieu geometrique
% des deux satellites, reprojete sur la sphere :
%
%   m_ij = (r_i + r_j)/||r_i+r_j||.
%
% Grandeurs calculees :
%   v_orb             : norme de la vitesse individuelle des satellites ;
%   v_rel             : ||v_j-v_i|| ;
%   v_rad_signed      : (v_j-v_i).e_ij ;
%   |v_rad|           : valeur absolue de la composante radiale ;
%   v_rad_out         : max(v_rad_signed,0) ;
%   v_rad_out_pos     : moyenne conditionnelle aux liens qui s'eloignent ;
%   v_rel et v_rad_out pour les liens dans la couche de rupture.
%
% Deux moyennes sont distinguees :
%   - moyenne temporelle par latitude :
%       chaque instant contenant au moins une observation a le meme poids ;
%   - moyenne agregee par latitude :
%       chaque satellite ou lien observe a le meme poids.
%
% Entree :
%   analysis_temp_results.mat
%
% Sortie :
%   vrel_vrad_emp_delta_latitude_results.mat

clear; clc; close all;

%% ============================================================
%  1. Chargement
%% ============================================================

script_dir = fileparts(mfilename('fullpath'));

candidate_files = {
    fullfile(script_dir, 'analysis_temp_results.mat')
    fullfile(script_dir, '..', 'analysis_temp_results.mat')
};

input_file = '';

for q = 1:numel(candidate_files)
    if isfile(candidate_files{q})
        input_file = candidate_files{q};
        break;
    end
end

if isempty(input_file)
    error(['Fichier analysis_temp_results.mat introuvable dans le ', ...
           'dossier du script ou dans son dossier parent.']);
end

S = load(input_file);

required_fields = {'Positions', 'Adjacency', 'dt', 'dmax'};

for q = 1:numel(required_fields)
    if ~isfield(S, required_fields{q})
        error('Le fichier doit contenir la variable %s.', ...
              required_fields{q});
    end
end

Positions = S.Positions;
Adjacency = S.Adjacency;
dt = double(S.dt);
dmax = double(S.dmax);

if ~iscell(Positions) || ~iscell(Adjacency)
    error('Positions et Adjacency doivent etre des cellules temporelles.');
end

Nt = min(numel(Positions), numel(Adjacency));

if Nt < 2
    error('Il faut au moins deux instants temporels.');
end

N = size(Positions{1},1);

for k = 1:Nt
    if ~isequal(size(Positions{k}), [N,3])
        error('Positions{%d} doit etre une matrice %d x 3.', k, N);
    end
end

if isfield(S,'time_values') && numel(S.time_values) >= Nt
    time_values = reshape(S.time_values(1:Nt), [], 1);
else
    time_values = (0:Nt-1).' * dt;
end

%% ============================================================
%  2. Parametres orbitaux de reference
%% ============================================================

if isfield(S,'R')
    R = double(S.R);
else
    R = mean(vecnorm(Positions{1},2,2));
end

if isfield(S,'mu')
    mu = double(S.mu);
else
    mu = 398600; % km^3/s^2
end

v_orb_theory = sqrt(mu/R);

if isfield(S,'inc')
    inc = double(S.inc);
elseif isfield(S,'inc_deg')
    inc = deg2rad(double(S.inc_deg));
else
    % Estimation par la latitude maximale effectivement observee.
    max_abs_lat = 0;
    for k = 1:Nt
        Pk = Positions{k};
        rk = vecnorm(Pk,2,2);
        lat_k = asin(max(min(Pk(:,3)./rk,1),-1));
        max_abs_lat = max(max_abs_lat,max(abs(lat_k)));
    end
    inc = max_abs_lat;
end

if inc > pi
    inc = deg2rad(inc);
end

inc_deg = rad2deg(inc);

%% ============================================================
%  3. Definition des tranches de latitude
%% ============================================================

% Largeur modifiable des tranches.
latitude_bin_width_deg = 5;

lat_limit_deg = ceil(inc_deg/latitude_bin_width_deg) ...
    * latitude_bin_width_deg;

if lat_limit_deg <= 0
    lat_limit_deg = 90;
end

latitude_edges_deg = ...
    -lat_limit_deg:latitude_bin_width_deg:lat_limit_deg;

if latitude_edges_deg(end) < lat_limit_deg
    latitude_edges_deg(end+1) = lat_limit_deg;
end

latitude_centers_deg = ...
    (latitude_edges_deg(1:end-1)+latitude_edges_deg(2:end))/2;

Nb = numel(latitude_centers_deg);

%% ============================================================
%  4. Vitesses instantanees
%
%  Difference centree aux instants interieurs ;
%  difference avant/arriere aux extremites.
%% ============================================================

Velocities = cell(Nt,1);

for k = 1:Nt
    if k == 1
        Velocities{k} = (Positions{2}-Positions{1})/dt;
    elseif k == Nt
        Velocities{k} = (Positions{Nt}-Positions{Nt-1})/dt;
    else
        Velocities{k} = ...
            (Positions{k+1}-Positions{k-1})/(2*dt);
    end
end

%% ============================================================
%  5. Stockage des moyennes instantanees par latitude
%% ============================================================

% Satellites
vorb_bin_t = nan(Nt,Nb);
n_sat_bin_t = zeros(Nt,Nb);

% Tous les liens
vrel_bin_t = nan(Nt,Nb);
vrad_signed_bin_t = nan(Nt,Nb);
vrad_abs_bin_t = nan(Nt,Nb);
vrad_out_bin_t = nan(Nt,Nb);
vrad_out_pos_bin_t = nan(Nt,Nb);
n_links_bin_t = zeros(Nt,Nb);
n_outgoing_links_bin_t = zeros(Nt,Nb);

% Liens dans la couche de rupture
vrel_border_bin_t = nan(Nt,Nb);
vrad_out_border_bin_t = nan(Nt,Nb);
n_border_links_bin_t = zeros(Nt,Nb);

% Sommes et effectifs pour les moyennes agregees
sum_vorb_bin = zeros(1,Nb);
count_vorb_bin = zeros(1,Nb);

sum_vrel_bin = zeros(1,Nb);
sum_vrad_signed_bin = zeros(1,Nb);
sum_vrad_abs_bin = zeros(1,Nb);
sum_vrad_out_bin = zeros(1,Nb);
count_links_bin = zeros(1,Nb);

sum_vrad_out_pos_bin = zeros(1,Nb);
count_outgoing_bin = zeros(1,Nb);

sum_vrel_border_bin = zeros(1,Nb);
sum_vrad_out_border_bin = zeros(1,Nb);
count_border_bin = zeros(1,Nb);

%% ============================================================
%  6. Calcul temporel par tranche de latitude
%% ============================================================

for k = 1:Nt

    P = double(Positions{k});
    V = double(Velocities{k});

    radius_sat = vecnorm(P,2,2);
    latitude_sat_deg = ...
        rad2deg(asin(max(min(P(:,3)./radius_sat,1),-1)));

    speed_sat = vecnorm(V,2,2);
    sat_bin = discretize(latitude_sat_deg,latitude_edges_deg);

    % --------------------------------------------------------
    % Vitesse orbitale par latitude
    % --------------------------------------------------------
    for b = 1:Nb
        mask_sat = sat_bin == b;
        n_sat_bin_t(k,b) = nnz(mask_sat);

        if any(mask_sat)
            values = speed_sat(mask_sat);
            vorb_bin_t(k,b) = mean(values,'omitnan');

            sum_vorb_bin(b) = ...
                sum_vorb_bin(b)+sum(values,'omitnan');
            count_vorb_bin(b) = ...
                count_vorb_bin(b)+nnz(isfinite(values));
        end
    end

    % --------------------------------------------------------
    % Liens existants
    % --------------------------------------------------------
    A = logical(spones(Adjacency{k}));
    A = A | A.';
    A(1:size(A,1)+1:end) = false;

    [I,J] = find(triu(A,1));

    if isempty(I)
        continue;
    end

    dr = P(J,:)-P(I,:);
    distance = vecnorm(dr,2,2);

    valid = isfinite(distance) & distance > 0;
    I = I(valid);
    J = J(valid);
    dr = dr(valid,:);
    distance = distance(valid);

    if isempty(I)
        continue;
    end

    e_rad = dr./distance;
    dv = V(J,:)-V(I,:);

    vrel = vecnorm(dv,2,2);
    vrad_signed = sum(dv.*e_rad,2);
    vrad_abs = abs(vrad_signed);
    vrad_out = max(vrad_signed,0);

    % Moyenne conditionnelle aux seuls liens qui s'eloignent.
    is_outgoing = vrad_signed > 0;

    % Couche de rupture sur un pas dt.
    is_border = is_outgoing & ...
        (distance + vrad_out*dt >= dmax);

    % Latitude du milieu geometrique du lien.
    midpoint = P(I,:)+P(J,:);
    midpoint_norm = vecnorm(midpoint,2,2);

    valid_midpoint = midpoint_norm > 0;
    latitude_link_deg = nan(size(midpoint_norm));

    midpoint_unit = midpoint(valid_midpoint,:) ...
        ./ midpoint_norm(valid_midpoint);

    latitude_link_deg(valid_midpoint) = ...
        rad2deg(asin(max(min(midpoint_unit(:,3),1),-1)));

    link_bin = discretize(latitude_link_deg,latitude_edges_deg);

    for b = 1:Nb
        mask_link = link_bin == b;
        n_links_bin_t(k,b) = nnz(mask_link);

        if any(mask_link)
            x_vrel = vrel(mask_link);
            x_signed = vrad_signed(mask_link);
            x_abs = vrad_abs(mask_link);
            x_out = vrad_out(mask_link);

            vrel_bin_t(k,b) = mean(x_vrel,'omitnan');
            vrad_signed_bin_t(k,b) = mean(x_signed,'omitnan');
            vrad_abs_bin_t(k,b) = mean(x_abs,'omitnan');
            vrad_out_bin_t(k,b) = mean(x_out,'omitnan');

            sum_vrel_bin(b) = ...
                sum_vrel_bin(b)+sum(x_vrel,'omitnan');
            sum_vrad_signed_bin(b) = ...
                sum_vrad_signed_bin(b)+sum(x_signed,'omitnan');
            sum_vrad_abs_bin(b) = ...
                sum_vrad_abs_bin(b)+sum(x_abs,'omitnan');
            sum_vrad_out_bin(b) = ...
                sum_vrad_out_bin(b)+sum(x_out,'omitnan');
            count_links_bin(b) = ...
                count_links_bin(b)+nnz(isfinite(x_vrel));
        end

        mask_out = mask_link & is_outgoing;
        n_outgoing_links_bin_t(k,b) = nnz(mask_out);

        if any(mask_out)
            x_out_pos = vrad_out(mask_out);
            vrad_out_pos_bin_t(k,b) = ...
                mean(x_out_pos,'omitnan');

            sum_vrad_out_pos_bin(b) = ...
                sum_vrad_out_pos_bin(b)+sum(x_out_pos,'omitnan');
            count_outgoing_bin(b) = ...
                count_outgoing_bin(b)+nnz(isfinite(x_out_pos));
        end

        mask_border = mask_link & is_border;
        n_border_links_bin_t(k,b) = nnz(mask_border);

        if any(mask_border)
            x_vrel_border = vrel(mask_border);
            x_vrad_border = vrad_out(mask_border);

            vrel_border_bin_t(k,b) = ...
                mean(x_vrel_border,'omitnan');
            vrad_out_border_bin_t(k,b) = ...
                mean(x_vrad_border,'omitnan');

            sum_vrel_border_bin(b) = ...
                sum_vrel_border_bin(b) ...
                + sum(x_vrel_border,'omitnan');

            sum_vrad_out_border_bin(b) = ...
                sum_vrad_out_border_bin(b) ...
                + sum(x_vrad_border,'omitnan');

            count_border_bin(b) = ...
                count_border_bin(b) ...
                + nnz(isfinite(x_vrel_border));
        end
    end
end

%% ============================================================
%  7. Moyennes temporelles par tranche
%
%  Chaque instant non vide a ici le meme poids.
%% ============================================================

vorb_time_mean_lat = mean(vorb_bin_t,1,'omitnan');

vrel_time_mean_lat = mean(vrel_bin_t,1,'omitnan');
vrad_signed_time_mean_lat = ...
    mean(vrad_signed_bin_t,1,'omitnan');
vrad_abs_time_mean_lat = ...
    mean(vrad_abs_bin_t,1,'omitnan');
vrad_out_time_mean_lat = ...
    mean(vrad_out_bin_t,1,'omitnan');
vrad_out_pos_time_mean_lat = ...
    mean(vrad_out_pos_bin_t,1,'omitnan');

vrel_border_time_mean_lat = ...
    mean(vrel_border_bin_t,1,'omitnan');
vrad_out_border_time_mean_lat = ...
    mean(vrad_out_border_bin_t,1,'omitnan');

%% ============================================================
%  7.b Modeles theoriques pour les liens proches de la rupture
%
%  Pour une paire locale constituee d'un satellite ascendant et
%  d'un satellite descendant a la latitude phi :
%
%    v_rel^th(phi)
%      = 2*v_orb*sqrt(sin(i)^2-sin(phi)^2)/cos(phi)
%
%  Le conditionnement par la couche de rupture pondere les grandes
%  projections radiales. Sous l'approximation d'une orientation
%  uniforme du lien relativement a v_rel :
%
%    v_rad,out^th(phi) = (pi/4)*v_rel^th(phi).
%% ============================================================

phi_centers_rad = deg2rad(latitude_centers_deg);

vrel_border_theory_lat = nan(1,Nb);
vrad_out_border_theory_lat = nan(1,Nb);

valid_theory = ...
    abs(phi_centers_rad) < inc & ...
    abs(cos(phi_centers_rad)) > 1e-12;

radicand = sin(inc)^2 - sin(phi_centers_rad(valid_theory)).^2;
radicand = max(radicand,0);

vrel_border_theory_lat(valid_theory) = ...
    2*v_orb_theory .* sqrt(radicand) ...
    ./ cos(phi_centers_rad(valid_theory));

vrad_out_border_theory_lat(valid_theory) = ...
    (pi/4) * vrel_border_theory_lat(valid_theory);

%% ============================================================
%  8. Moyennes agregees par tranche
%
%  Chaque observation individuelle a ici le meme poids.
%% ============================================================

vorb_aggregated_mean_lat = safe_divide( ...
    sum_vorb_bin,count_vorb_bin);

vrel_aggregated_mean_lat = safe_divide( ...
    sum_vrel_bin,count_links_bin);
vrad_signed_aggregated_mean_lat = safe_divide( ...
    sum_vrad_signed_bin,count_links_bin);
vrad_abs_aggregated_mean_lat = safe_divide( ...
    sum_vrad_abs_bin,count_links_bin);
vrad_out_aggregated_mean_lat = safe_divide( ...
    sum_vrad_out_bin,count_links_bin);

vrad_out_pos_aggregated_mean_lat = safe_divide( ...
    sum_vrad_out_pos_bin,count_outgoing_bin);

vrel_border_aggregated_mean_lat = safe_divide( ...
    sum_vrel_border_bin,count_border_bin);
vrad_out_border_aggregated_mean_lat = safe_divide( ...
    sum_vrad_out_border_bin,count_border_bin);

%% ============================================================
%  9. Rapports diagnostiques dans la couche de rupture
%% ============================================================

% Rapport entre la vitesse relative des liens proches de la rupture
% et la vitesse orbitale locale moyenne.
ratio_vrelborder_vorb_time_lat = ...
    vrel_border_time_mean_lat./vorb_time_mean_lat;

% Projection radiale sortante conditionnellement aux liens proches
% de la rupture.
ratio_vradborder_vrelborder_time_lat = ...
    vrad_out_border_time_mean_lat./vrel_border_time_mean_lat;

ratio_vrel_emp_theory_lat = ...
    vrel_border_time_mean_lat./vrel_border_theory_lat;

ratio_vrad_emp_theory_lat = ...
    vrad_out_border_time_mean_lat./vrad_out_border_theory_lat;

relative_error_vrel_lat = ...
    (vrel_border_time_mean_lat-vrel_border_theory_lat) ...
    ./ vrel_border_theory_lat;

relative_error_vrad_lat = ...
    (vrad_out_border_time_mean_lat-vrad_out_border_theory_lat) ...
    ./ vrad_out_border_theory_lat;

%% ============================================================
%  10. Valeurs globales dans la couche de rupture
%% ============================================================

v_orb_emp_global = ...
    sum(sum_vorb_bin)/max(sum(count_vorb_bin),1);

v_rel_border_emp_global = ...
    sum(sum_vrel_border_bin)/max(sum(count_border_bin),1);

v_rad_out_border_emp_global = ...
    sum(sum_vrad_out_border_bin)/max(sum(count_border_bin),1);

ratio_vrelborder_vorb_global = ...
    v_rel_border_emp_global/v_orb_emp_global;

ratio_vradborder_vrelborder_global = ...
    v_rad_out_border_emp_global/v_rel_border_emp_global;

%% ============================================================
%  11. Tableau recapitulatif
%% ============================================================

latitude_min_deg = latitude_edges_deg(1:end-1).';
latitude_max_deg = latitude_edges_deg(2:end).';
latitude_center_deg = latitude_centers_deg.';

results_by_latitude = table( ...
    latitude_min_deg, ...
    latitude_max_deg, ...
    latitude_center_deg, ...
    count_vorb_bin.', ...
    count_border_bin.', ...
    vorb_time_mean_lat.', ...
    vrel_border_time_mean_lat.', ...
    vrad_out_border_time_mean_lat.', ...
    vrel_border_theory_lat.', ...
    vrad_out_border_theory_lat.', ...
    ratio_vrelborder_vorb_time_lat.', ...
    ratio_vradborder_vrelborder_time_lat.', ...
    ratio_vrel_emp_theory_lat.', ...
    ratio_vrad_emp_theory_lat.', ...
    relative_error_vrel_lat.', ...
    relative_error_vrad_lat.', ...
    'VariableNames', { ...
    'latitude_min_deg', ...
    'latitude_max_deg', ...
    'latitude_center_deg', ...
    'n_sat_observations', ...
    'n_border_link_observations', ...
    'v_orb_time_mean_km_s', ...
    'v_rel_border_time_mean_km_s', ...
    'v_rad_out_border_time_mean_km_s', ...
    'v_rel_border_theory_km_s', ...
    'v_rad_out_border_theory_km_s', ...
    'ratio_vrelborder_vorb', ...
    'ratio_vradborder_vrelborder', ...
    'ratio_vrel_emp_theory', ...
    'ratio_vrad_emp_theory', ...
    'relative_error_vrel', ...
    'relative_error_vrad'});

disp(results_by_latitude);

%% ============================================================
%  12. Affichage numerique global
%% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' DIAGNOSTIC CINEMATIQUE - LIENS PROCHES DE LA RUPTURE\n');
fprintf('============================================================\n');
fprintf('Fichier charge                         : %s\n',input_file);
fprintf('Nombre de satellites                  : %d\n',N);
fprintf('Nombre d''instants                    : %d\n',Nt);
fprintf('Inclinaison utilisee                  : %.6f deg\n',inc_deg);
fprintf('Largeur des tranches                  : %.6f deg\n', ...
    latitude_bin_width_deg);
fprintf('------------------------------------------------------------\n');
fprintf('v_orb theorique = sqrt(mu/R)          : %.8f km/s\n', ...
    v_orb_theory);
fprintf('v_orb empirique globale              : %.8f km/s\n', ...
    v_orb_emp_global);
fprintf('v_rel | couche de rupture            : %.8f km/s\n', ...
    v_rel_border_emp_global);
fprintf('v_rad,out | couche de rupture        : %.8f km/s\n', ...
    v_rad_out_border_emp_global);
fprintf('------------------------------------------------------------\n');
fprintf('v_rel,bord / v_orb                   : %.8f\n', ...
    ratio_vrelborder_vorb_global);
fprintf('v_rad,out,bord / v_rel,bord          : %.8f\n', ...
    ratio_vradborder_vrelborder_global);
fprintf('Reference 1/pi                       : %.8f\n',1/pi);
fprintf('Reference 2/pi                       : %.8f\n',2/pi);
fprintf('Reference 4/pi                       : %.8f\n',4/pi);
fprintf('Reference 1/sqrt(2)                  : %.8f\n',1/sqrt(2));
fprintf('============================================================\n');

%% ============================================================
%  13. Figures centrees sur les liens proches de la rupture
%% ============================================================

figure;
plot(latitude_centers_deg,vorb_time_mean_lat,'LineWidth',1.5);
hold on;
plot(latitude_centers_deg,vrel_border_time_mean_lat,'LineWidth',1.5);
plot(latitude_centers_deg,vrel_border_theory_lat,'--','LineWidth',1.7);
plot(latitude_centers_deg,vrad_out_border_time_mean_lat,'LineWidth',1.5);
plot(latitude_centers_deg,vrad_out_border_theory_lat,'--','LineWidth',1.7);
yline(v_orb_theory,':', ...
    sprintf('v_{orb}^{th}=%.3f km/s',v_orb_theory));
grid on;
xlabel('Latitude du milieu du lien (deg)');
ylabel('Vitesse moyenne temporelle (km/s)');
title('Vitesses proches de la rupture : empirique et theorie');
legend('v_{orb} empirique', ...
       'v_{rel} empirique | rupture', ...
       'v_{rel}^{th}(\varphi)', ...
       'v_{rad,out} empirique | rupture', ...
       'v_{rad,out}^{th}(\varphi)', ...
       'v_{orb}^{th}', ...
       'Location','best');
hold off;

figure;
plot(latitude_centers_deg,ratio_vrelborder_vorb_time_lat, ...
    'LineWidth',1.5);
hold on;
plot(latitude_centers_deg, ...
    vrel_border_theory_lat./v_orb_theory, ...
    '--','LineWidth',1.7);
grid on;
xlabel('Latitude du milieu du lien (deg)');
ylabel('<v_{rel}|rupture>/v_{orb}');
title('Modele de vitesse relative selon la latitude');
legend('Rapport empirique', ...
       'Modele theorique', ...
       'Location','best');
hold off;

figure;
plot(latitude_centers_deg,ratio_vradborder_vrelborder_time_lat, ...
    'LineWidth',1.5);
hold on;
yline(pi/4,'--','\pi/4','LineWidth',1.7);
grid on;
xlabel('Latitude du milieu du lien (deg)');
ylabel('<v_{rad,out}|rupture>/<v_{rel}|rupture>');
title('Modele de projection radiale dans la couche de rupture');
legend('Rapport empirique','Modele \pi/4', ...
    'Location','best');
hold off;

figure;
plot(latitude_centers_deg,ratio_vrel_emp_theory_lat, ...
    'LineWidth',1.5);
hold on;
plot(latitude_centers_deg,ratio_vrad_emp_theory_lat, ...
    'LineWidth',1.5);
yline(1,'--','Accord parfait','LineWidth',1.5);
grid on;
xlabel('Latitude du milieu du lien (deg)');
ylabel('Valeur empirique / valeur theorique');
title('Qualite des modeles cinematiques selon la latitude');
legend('v_{rel}^{emp}/v_{rel}^{th}', ...
       'v_{rad,out}^{emp}/v_{rad,out}^{th}', ...
       'Location','best');
hold off;

figure;
plot(latitude_centers_deg,count_border_bin,'LineWidth',1.5);
grid on;
xlabel('Latitude du milieu du lien (deg)');
ylabel('Nombre total d''observations');
title('Nombre de liens proches de la rupture par latitude');

%% ============================================================
%  14. Sauvegarde
%% ============================================================

output_file = fullfile( ...
    script_dir,'vrel_vrad_temp_results.mat');

save(output_file, ...
    'input_file', ...
    'time_values', ...
    'R','mu','inc','inc_deg','dmax','dt', ...
    'v_orb_theory', ...
    'latitude_bin_width_deg', ...
    'latitude_edges_deg', ...
    'latitude_centers_deg', ...
    'vorb_bin_t', ...
    'vrel_border_bin_t', ...
    'vrad_out_border_bin_t', ...
    'n_sat_bin_t', ...
    'n_border_links_bin_t', ...
    'vorb_time_mean_lat', ...
    'vrel_border_time_mean_lat', ...
    'vrad_out_border_time_mean_lat', ...
    'phi_centers_rad', ...
    'vrel_border_theory_lat', ...
    'vrad_out_border_theory_lat', ...
    'ratio_vrel_emp_theory_lat', ...
    'ratio_vrad_emp_theory_lat', ...
    'relative_error_vrel_lat', ...
    'relative_error_vrad_lat', ...
    'vorb_aggregated_mean_lat', ...
    'vrel_border_aggregated_mean_lat', ...
    'vrad_out_border_aggregated_mean_lat', ...
    'ratio_vrelborder_vorb_time_lat', ...
    'ratio_vradborder_vrelborder_time_lat', ...
    'v_orb_emp_global', ...
    'v_rel_border_emp_global', ...
    'v_rad_out_border_emp_global', ...
    'ratio_vrelborder_vorb_global', ...
    'ratio_vradborder_vrelborder_global', ...
    'results_by_latitude');

fprintf('\nResultats sauvegardes dans :\n%s\n',output_file);

%% ============================================================
%  Fonction locale
%% ============================================================

function ratio = safe_divide(numerator,denominator)
    ratio = nan(size(numerator));
    valid = denominator > 0;
    ratio(valid) = numerator(valid)./denominator(valid);
end
