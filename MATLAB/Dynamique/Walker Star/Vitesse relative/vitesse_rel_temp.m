clear; clc; close all;

%% ============================================================
%  ANGLE EFFECTIF ENTRE VECTEURS VITESSES CONDITIONNE AUX LIENS
%  VERSION CORRIGEE
%
%  Correction importante : on construit explicitement les couples (i,j)
%  une seule fois, puis on utilise EXACTEMENT le meme ordre de couples pour :
%  - les distances d_ij,
%  - les vitesses relatives ||v_i-v_j||,
%  - les angles gamma_ij.
%
%  Cela evite le piege d'ordre entre pdist(...) et l'extraction par
%  matrice triangulaire, qui peut decaler les couples et rendre incoherente
%  la verification.
%% ============================================================

%% Parametres physiques
R_earth = 6371;      % km
h = 550;             % km
R = R_earth + h;     % rayon orbital

mu = 398600;              % km^3/s^2
omega = sqrt(mu / R^3);   % vitesse angulaire orbitale rad/s
v_orb = sqrt(mu / R);     % km/s

%% Parametres du processus de Poisson
lambda = 4e-7;       % satellites / km^2
surface_sphere = 4*pi*R^2;

N = poissrnd(lambda * surface_sphere);

fprintf('Nombre de satellites generes : N = %d\n', N);
fprintf('v_orb = %.4f km/s\n', v_orb);
fprintf('v_rel isotrope 4/pi v_orb = %.4f km/s\n', (4/pi)*v_orb);

%% Generation uniforme des positions initiales sur la sphere
u = rand(N,1);
phi = 2*pi*rand(N,1);
theta = acos(1 - 2*u);

%% Sens de rotation defini par un plan separateur passant par les poles
rotation_sign = ones(N, 1);
rotation_sign(y_from_angles(theta,phi) >= 0) = 1;
rotation_sign(y_from_angles(theta,phi) < 0) = -1;

%% Parametres des liens et du temps
dmax = 1500;      % km
dt = 60;          % pas temporel en secondes
Tmax = 12000;     % duree totale de simulation

time_values = 0:dt:Tmax;
Nt = length(time_values);

%% Couples i < j
% On utilise ces indices pour TOUTES les grandeurs.
[pair_i, pair_j] = find(triu(true(N), 1));
Npairs = numel(pair_i);

%% Stockage
gamma_mean_raw = NaN(Nt,1);       % E[gamma | lien], rad, diagnostic seulement
gamma_std_raw  = NaN(Nt,1);       % std(gamma | lien), rad
gamma_med_raw  = NaN(Nt,1);       % mediane(gamma | lien), rad

E_sin_half_link = NaN(Nt,1);      % E[sin(gamma/2) | lien]
gamma_eff_link = NaN(Nt,1);       % 2 asin(E[sin(gamma/2) | lien]), rad

vrel_link_direct = NaN(Nt,1);     % mean(||v_i-v_j|| | lien), km/s
vrel_link_from_gamma = NaN(Nt,1); % 2 v_orb E[sin(gamma/2) | lien], km/s
consistency_error = NaN(Nt,1);    % verification numerique

nb_links = zeros(Nt,1);

%% Strates d'altitude absolue du milieu des liens
% 0 deg = equateur, 90 deg = pole.
%
% Deux discretisations sont utilisees :
%   1) alt_edges_emp_deg : pour le comptage empirique des liens par altitude ;
%   2) alt_edges_th_deg  : pour la theorie, chargee depuis liens_inter.mat.
%
% La theorie du nombre de liens par strate ne vient plus du comptage
% empirique fait ici. Elle est chargee depuis liens_inter.mat, via
% mean_links_incident_strate_t.

n_alt_bins_emp = 30;
alt_edges_emp_deg = linspace(0, 90, n_alt_bins_emp+1);
alt_centers_emp_deg = 0.5 * (alt_edges_emp_deg(1:end-1) + alt_edges_emp_deg(2:end));
link_counts_alt_time_emp = zeros(Nt, n_alt_bins_emp);

%% Chargement du nombre THEORIQUE de liens par strate depuis liens_inter.mat
if exist('liens_inter.mat', 'file') ~= 2
    error(['Fichier liens_inter.mat introuvable. ', ...
           'Place-le dans le dossier courant MATLAB.']);
end

S_liens = load('liens_inter.mat');

required_vars = {'mean_links_incident_strate_t', 'time_values', ...
                 'beta_step_strates', 'beta_max_strates'};
for vv = 1:numel(required_vars)
    if ~isfield(S_liens, required_vars{vv})
        error('La variable %s est absente de liens_inter.mat.', required_vars{vv});
    end
end

% Dans liens_inter.mat, les strates sont indexees depuis le pole vers
% l'equateur. Pour gamma_lien_theorique_strates, on veut des altitudes
% absolues croissantes : 0 deg = equateur, 90 deg = pole.
L_th_pole_to_equator = double(S_liens.mean_links_incident_strate_t);
time_values_liens = double(S_liens.time_values(:));

M_th = size(L_th_pole_to_equator, 2);

beta_step_th = double(S_liens.beta_step_strates);
beta_max_th  = double(S_liens.beta_max_strates);

beta_edges_th = 0:beta_step_th:beta_max_th;
if beta_edges_th(end) < beta_max_th - 1e-12
    beta_edges_th = [beta_edges_th, beta_max_th];
end
% Securite : on force exactement M_th+1 bords.
if numel(beta_edges_th) ~= M_th + 1
    beta_edges_th = linspace(0, beta_max_th, M_th+1);
end

% Conversion colatitude depuis le pole beta -> altitude absolue depuis
% l'equateur ell = pi/2 - beta, puis ordre croissant.
alt_edges_th_rad = fliplr(pi/2 - beta_edges_th);
alt_edges_th_rad(1) = max(0, alt_edges_th_rad(1));
alt_edges_th_rad(end) = min(pi/2, alt_edges_th_rad(end));
alt_edges_th_deg = rad2deg(alt_edges_th_rad);
alt_centers_th_deg = 0.5 * (alt_edges_th_deg(1:end-1) + alt_edges_th_deg(2:end));

% Meme inversion pour les liens : colonnes equateur -> pole.
L_alt_th_from_mat = fliplr(L_th_pole_to_equator);

% Interpolation temporelle si la grille de liens_inter.mat differe.
if numel(time_values_liens) ~= Nt || any(abs(time_values_liens(:) - time_values(:)) > 1e-9)
    L_alt_th_time = interp1(time_values_liens, L_alt_th_from_mat, ...
                            time_values(:), 'linear', 'extrap');
    L_alt_th_time = max(L_alt_th_time, 0);
else
    L_alt_th_time = L_alt_th_from_mat;
end


%% ============================================================
%  CALCUL TEMPOREL
%% ============================================================

for k = 1:Nt

    t = time_values(k);

    %% Mouvement orbital
    phi_t = phi;
    theta_t = theta + rotation_sign * omega * t;
    theta_dot = rotation_sign * omega;

    %% Positions
    x_t = R * sin(theta_t) .* cos(phi_t);
    y_t = R * sin(theta_t) .* sin(phi_t);
    z_t = R * cos(theta_t);
    positions_t = [x_t y_t z_t];

    %% Vitesses orbitales
    vx_t = R * theta_dot .* cos(theta_t) .* cos(phi_t);
    vy_t = R * theta_dot .* cos(theta_t) .* sin(phi_t);
    vz_t = -R * theta_dot .* sin(theta_t);
    velocities_t = [vx_t vy_t vz_t];

    %% Distances, vitesses relatives et angles, avec le MEME ordre de couples
    dpos = positions_t(pair_i,:) - positions_t(pair_j,:);
    D_pairs = sqrt(sum(dpos.^2, 2));

    dvel = velocities_t(pair_i,:) - velocities_t(pair_j,:);
    V_pairs = sqrt(sum(dvel.^2, 2));

    vi = velocities_t(pair_i,:);
    vj = velocities_t(pair_j,:);

    ni = sqrt(sum(vi.^2, 2));
    nj = sqrt(sum(vj.^2, 2));

    cos_gamma_pairs = sum(vi .* vj, 2) ./ (ni .* nj);
    cos_gamma_pairs = max(-1, min(1, cos_gamma_pairs));
    gamma_pairs = acos(cos_gamma_pairs);  % rad, dans [0, pi]

    %% Condition de lien
    link_mask = (D_pairs <= dmax);
    nb_links(k) = sum(link_mask);

    if nb_links(k) > 0
        gamma_links = gamma_pairs(link_mask);
        vrel_links = V_pairs(link_mask);

        %% Altitude absolue du milieu de chaque lien
        ri_links = positions_t(pair_i(link_mask),:);
        rj_links = positions_t(pair_j(link_mask),:);
        rmid_links = 0.5 * (ri_links + rj_links);
        rmid_norm = sqrt(sum(rmid_links.^2, 2));

        alt_links_rad = abs(asin(max(-1, min(1, rmid_links(:,3) ./ rmid_norm))));
        alt_links_deg = rad2deg(alt_links_rad);

        % Comptage EMPIRIQUE des liens par strate d'altitude.
        % Ce comptage sert uniquement au diagnostic / comparaison.
        % La theorie utilise L_alt_th_time charge depuis liens_inter.mat.
        alt_bin_idx = discretize(alt_links_deg, alt_edges_emp_deg);
        for mm = 1:n_alt_bins_emp
            link_counts_alt_time_emp(k,mm) = sum(alt_bin_idx == mm);
        end


        %% Diagnostic : moyenne brute des angles
        gamma_mean_raw(k) = mean(gamma_links, 'omitnan');
        gamma_std_raw(k)  = std(gamma_links, 'omitnan');
        gamma_med_raw(k)  = median(gamma_links, 'omitnan');

        %% Quantite pertinente pour la vitesse relative moyenne
        E_sin_half_link(k) = mean(sin(gamma_links/2), 'omitnan');
        E_sin_half_link(k) = max(0, min(1, E_sin_half_link(k)));

        gamma_eff_link(k) = 2 * asin(E_sin_half_link(k));

        %% Verification directe par les vitesses relatives
        vrel_link_direct(k) = mean(vrel_links, 'omitnan');
        vrel_link_from_gamma(k) = 2 * v_orb * E_sin_half_link(k);

        consistency_error(k) = abs(vrel_link_direct(k) - vrel_link_from_gamma(k));
    end
end

%% Conversion en degres
gamma_mean_raw_deg = rad2deg(gamma_mean_raw);
gamma_std_raw_deg  = rad2deg(gamma_std_raw);
gamma_med_raw_deg  = rad2deg(gamma_med_raw);
gamma_eff_link_deg = rad2deg(gamma_eff_link);

gamma_eff_mean_deg = mean(gamma_eff_link_deg, 'omitnan');
gamma_raw_mean_time_deg = mean(gamma_mean_raw_deg, 'omitnan');
E_sin_half_mean = mean(E_sin_half_link, 'omitnan');
vrel_direct_mean = mean(vrel_link_direct, 'omitnan');
vrel_from_gamma_mean = mean(vrel_link_from_gamma, 'omitnan');
max_consistency_error = max(consistency_error, [], 'omitnan');


%% ============================================================
%  MODELE THEORIQUE PAR STRATES D'ALTITUDE
%% ============================================================

% Reconstruction theorique de la vitesse relative a partir :
%   - du nombre de liens THEORIQUE par strate charge depuis liens_inter.mat ;
%   - du facteur angulaire analytique s_m.
[E_sin_half_th, gamma_eff_th, gamma_eff_th_deg, vrel_link_theory_strates, ...
    alt_centers_deg_th, s_alt_th, Phi_alt_th, L_alt_th_time] = ...
    gamma_lien_theorique_strates(L_alt_th_time, alt_edges_th_deg, R, dmax, v_orb);

% Reconstruction semi-analytique de diagnostic, utilisant les liens
% empiriques par strate calcules dans ce script.
[E_sin_half_emp_strates, gamma_eff_emp_strates, gamma_eff_emp_strates_deg, ...
    vrel_emp_strates, alt_centers_emp_deg_diag, s_alt_emp_diag, Phi_alt_emp_diag, ...
    link_counts_alt_time_emp] = ...
    gamma_lien_theorique_strates(link_counts_alt_time_emp, alt_edges_emp_deg, R, dmax, v_orb);

E_sin_half_th_mean = mean(E_sin_half_th, 'omitnan');
gamma_eff_th_mean_deg = mean(gamma_eff_th_deg, 'omitnan');
vrel_link_theory_mean = mean(vrel_link_theory_strates, 'omitnan');

fprintf('\n--- Resultats temporels conditionnes aux liens ---\n');
fprintf('Nombre moyen de liens                                      : %.2f\n', mean(nb_links));
fprintf('Angle effectif moyen gamma_eff                             : %.4f deg\n', gamma_eff_mean_deg);
fprintf('Angle brut moyen E[gamma | lien]                            : %.4f deg\n', gamma_raw_mean_time_deg);
fprintf('E[sin(gamma/2) | lien] moyen                                : %.6f\n', E_sin_half_mean);
fprintf('Vitesse relative directe moyenne                            : %.4f km/s\n', vrel_direct_mean);
fprintf('Vitesse deduite de E[sin(gamma/2)] moyenne                  : %.4f km/s\n', vrel_from_gamma_mean);
fprintf('Erreur max |vrel_direct - vrel_angle|                       : %.3e km/s\n', max_consistency_error);
fprintf('Vitesse isotrope 4/pi v_orb                                 : %.4f km/s\n', (4/pi)*v_orb);
fprintf('--- Modele theorique par strates ---\n');
fprintf('E[sin(gamma/2)] theorique moyen par strates                  : %.6f\n', E_sin_half_th_mean);
fprintf('Angle effectif theorique moyen par strates                   : %.4f deg\n', gamma_eff_th_mean_deg);
fprintf('Vitesse relative theorique moyenne par strates               : %.4f km/s\n', vrel_link_theory_mean);
fprintf('  -> L_m^th(t) charge depuis liens_inter.mat                  \n');
fprintf('Vitesse semi-analytique moyenne avec L_m empirique            : %.4f km/s\n', mean(vrel_emp_strates, 'omitnan'));

%% ============================================================
%  GRAPHE 1 : ANGLE EFFECTIF
%% ============================================================

figure;
plot(time_values, gamma_eff_link_deg, 'LineWidth', 1.5);
hold on;
plot(time_values, gamma_eff_th_deg, '--', 'LineWidth', 1.5);
yline(gamma_eff_mean_deg, ':', ...
    sprintf('moyenne \\gamma_{eff} = %.2f deg', gamma_eff_mean_deg), ...
    'LineWidth', 1.2);
grid on;
xlabel('Temps (s)');
ylabel('\gamma_{eff,link}(t) (degres)');
title('Angle effectif des vitesses conditionne a l''existence d''un lien');
legend('\gamma_{eff} empirique', ...
       '\gamma_{eff} theorique strates', ...
       'Moyenne temporelle empirique', ...
       'Location', 'best');

%% ============================================================
%  GRAPHE 2 : DIAGNOSTIC ANGLE BRUT
%% ============================================================

figure;
plot(time_values, gamma_mean_raw_deg, 'LineWidth', 1.2);
hold on;
plot(time_values, gamma_mean_raw_deg + gamma_std_raw_deg, ':', 'LineWidth', 1.0);
plot(time_values, gamma_mean_raw_deg - gamma_std_raw_deg, ':', 'LineWidth', 1.0);
plot(time_values, gamma_med_raw_deg, '--', 'LineWidth', 1.0);
yline(gamma_raw_mean_time_deg, '-.', ...
    sprintf('moyenne brute = %.2f deg', gamma_raw_mean_time_deg), ...
    'LineWidth', 1.0);
grid on;
xlabel('Temps (s)');
ylabel('\gamma_{link}(t) brut (degres)');
title('Diagnostic : moyenne brute de l''angle conditionne aux liens');
legend('E[\gamma | lien]', ...
       'Moyenne + ecart-type', ...
       'Moyenne - ecart-type', ...
       'Mediane', ...
       'Moyenne temporelle', ...
       'Location', 'best');

%% ============================================================
%  GRAPHE 3 : E[sin(gamma/2) | lien]
%% ============================================================

figure;
plot(time_values, E_sin_half_link, 'LineWidth', 1.5);
hold on;
plot(time_values, E_sin_half_th, '--', 'LineWidth', 1.5);
yline(E_sin_half_mean, ':', ...
    sprintf('moyenne = %.4f', E_sin_half_mean), ...
    'LineWidth', 1.2);
yline(2/pi, '--', ...
    sprintf('cas isotrope 2/\\pi = %.4f', 2/pi), ...
    'LineWidth', 1.2);
grid on;
xlabel('Temps (s)');
ylabel('E[sin(\gamma/2) | d_{ij} \leq d_{max}]');
title('Esperance conditionnelle pertinente pour la vitesse relative');
legend('Empirique conditionne aux liens', ...
       'Theorie par strates avec L_m^{th} liens\_inter', ...
       'Semi-analytique avec L_m empirique', ...
       'Moyenne temporelle empirique', ...
       'Reference isotrope', ...
       'Location', 'best');

%% ============================================================
%  GRAPHE 4 : VERIFICATION DE LA VITESSE RELATIVE
%% ============================================================

figure;
plot(time_values, vrel_link_direct, 'LineWidth', 1.5);
hold on;
plot(time_values, vrel_link_from_gamma, '--', 'LineWidth', 1.2);
plot(time_values, vrel_link_theory_strates, '-.', 'LineWidth', 1.7);
plot(time_values, vrel_emp_strates, ':', 'LineWidth', 1.2);
yline(vrel_direct_mean, ':', ...
    sprintf('moyenne = %.2f km/s', vrel_direct_mean), ...
    'LineWidth', 1.2);
yline((4/pi)*v_orb, '--', ...
    sprintf('4/\\pi v_{orb} = %.2f km/s', (4/pi)*v_orb), ...
    'LineWidth', 1.2);
grid on;
xlabel('Temps (s)');
ylabel('Vitesse relative moyenne des liens (km/s)');
title('Verification : vitesse directe et vitesse deduite des angles');
legend('Mean ||v_i - v_j|| | lien', ...
       '2 v_{orb} E[sin(\gamma/2) | lien]', ...
       'Theorie par strates avec L_m^{th} liens\_inter', ...
       'Semi-analytique avec L_m empirique', ...
       'Moyenne temporelle empirique', ...
       'Reference isotrope', ...
       'Location', 'best');


%% ============================================================
%  GRAPHE 5 : FACTEUR ANGULAIRE THEORIQUE PAR STRATE
%% ============================================================

figure;
plot(alt_centers_deg_th, s_alt_th, 'LineWidth', 1.5);
hold on;
yline(2/pi, '--', ...
    sprintf('cas isotrope 2/\\pi = %.4f', 2/pi), ...
    'LineWidth', 1.2);
grid on;
xlabel('Altitude absolue du lien sur la demi-sphere (degres)');
ylabel('s_m^{th} = E[sin(\gamma/2) | strate m]');
title('Facteur angulaire theorique par strate d''altitude');
legend('Modele geometrique par strate', ...
       'Reference isotrope', ...
       'Location', 'best');


%% ============================================================
%  GRAPHE 6 : LIENS THEORIQUES PAR STRATE CHARGES DE liens_inter.mat
%% ============================================================

% Comparaison au pic empirique de liens de la simulation courante.
[~, k_peak_emp] = max(nb_links);
[~, k_peak_th]  = max(sum(L_alt_th_time, 2));

L_emp_peak_alt = link_counts_alt_time_emp(k_peak_emp, :);
% Pour comparer aux 8 strates theoriques, on recompte empiriquement sur
% les bords de strates de liens_inter.mat.
L_emp_peak_thbins = zeros(1, M_th);

% Recalcul léger des liens empiriques au pic, sur les bords theoriques.
t = time_values(k_peak_emp);
phi_t = phi;
theta_t = theta + rotation_sign * omega * t;
theta_dot = rotation_sign * omega;

x_t = R * sin(theta_t) .* cos(phi_t);
y_t = R * sin(theta_t) .* sin(phi_t);
z_t = R * cos(theta_t);
positions_t = [x_t y_t z_t];

vx_t = R * theta_dot .* cos(theta_t) .* cos(phi_t);
vy_t = R * theta_dot .* cos(theta_t) .* sin(phi_t);
vz_t = -R * theta_dot .* sin(theta_t);
velocities_t = [vx_t vy_t vz_t]; %#ok<NASGU>

dpos = positions_t(pair_i,:) - positions_t(pair_j,:);
D_pairs = sqrt(sum(dpos.^2, 2));
link_mask = (D_pairs <= dmax);

if any(link_mask)
    ri_links = positions_t(pair_i(link_mask),:);
    rj_links = positions_t(pair_j(link_mask),:);
    rmid_links = 0.5 * (ri_links + rj_links);
    rmid_norm = sqrt(sum(rmid_links.^2, 2));
    alt_links_rad = abs(asin(max(-1, min(1, rmid_links(:,3) ./ rmid_norm))));
    alt_links_deg = rad2deg(alt_links_rad);

    alt_bin_idx_th = discretize(alt_links_deg, alt_edges_th_deg);
    for mm = 1:M_th
        L_emp_peak_thbins(mm) = sum(alt_bin_idx_th == mm);
    end
end

L_th_peak = L_alt_th_time(k_peak_th, :);

active_mask = (L_emp_peak_thbins > 0) | (L_th_peak > 0.5);
active_idx = find(active_mask);

figure;
bar(1:numel(active_idx), [L_emp_peak_thbins(active_idx).', L_th_peak(active_idx).']);
grid on;
xlabel('Indice de strate active (equateur \rightarrow pole)');
ylabel('Nombre de liens');
title(sprintf('Liens par strate : empirique au pic vs theorie liens\\_inter.mat'));
legend('Empirique simulation courante', 'Theorie liens\_inter.mat', 'Location', 'best');

%% Sauvegarde
save('angle_vitesses_lien_temp_corrige_v2_results.mat', ...
    'N', 'R', 'h', 'lambda', 'dmax', 'dt', 'Tmax', ...
    'time_values', ...
    'v_orb', ...
    'gamma_mean_raw', 'gamma_std_raw', 'gamma_med_raw', ...
    'gamma_mean_raw_deg', 'gamma_std_raw_deg', 'gamma_med_raw_deg', ...
    'gamma_eff_link', 'gamma_eff_link_deg', 'gamma_eff_mean_deg', ...
    'gamma_raw_mean_time_deg', ...
    'E_sin_half_link', 'E_sin_half_mean', ...
    'vrel_link_direct', 'vrel_link_from_gamma', ...
    'vrel_direct_mean', 'vrel_from_gamma_mean', ...
    'consistency_error', 'max_consistency_error', ...
    'nb_links', ...
    'n_alt_bins_emp', 'alt_edges_emp_deg', 'alt_centers_emp_deg', ...
    'link_counts_alt_time_emp', ...
    'alt_edges_th_deg', 'alt_centers_th_deg', ...
    'L_alt_th_time', 'L_alt_th_from_mat', 'time_values_liens', ...
    'E_sin_half_th', 'E_sin_half_th_mean', ...
    'gamma_eff_th', 'gamma_eff_th_deg', 'gamma_eff_th_mean_deg', ...
    'vrel_link_theory_strates', 'vrel_link_theory_mean', ...
    'E_sin_half_emp_strates', 'gamma_eff_emp_strates', ...
    'gamma_eff_emp_strates_deg', 'vrel_emp_strates', ...
    'alt_centers_deg_th', 's_alt_th', 'Phi_alt_th');

fprintf('\nResultats sauvegardes dans angle_vitesses_lien_temp_corrige_v2_results.mat\n');

%% ============================================================
%  FONCTION LOCALE
%% ============================================================

function y = y_from_angles(theta, phi)
    y = sin(theta) .* sin(phi);
end
