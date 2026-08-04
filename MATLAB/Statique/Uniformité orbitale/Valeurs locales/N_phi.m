%% satellites_phi.m
% Comparaison theorique / Monte-Carlo du nombre moyen de satellites
% en fonction de la latitude.
%
% Entree :
%   plink_phi_results.mat produit par plink_phi_mc_corrige.m
%
% La valeur N utilisee dans la theorie est la moyenne des N tires
% pendant les simulations Monte-Carlo.
%
% Theorie :
%   N^Delta(phi) = N f_Phi(phi)
%
% Pour une tranche b :
%   N_b^th = N P(Phi dans la tranche b)

clear; clc; close all;

%% ============================================================
%  1. Chargement
%% ============================================================

script_dir = fileparts(mfilename('fullpath'));
input_file = fullfile(script_dir,'plink_phi_results.mat');

if ~isfile(input_file)
    input_file = fullfile(script_dir,'..','plink_phi_results.mat');
end

if ~isfile(input_file)
    error('Fichier plink_phi_results.mat introuvable.');
end

S = load(input_file);

required_fields = { ...
    'N','phi_vals','f_phi','dphi_bins', ...
    'phi_bin_probability','count_sat_by_bin','n_sim'};

for k = 1:numel(required_fields)
    if ~isfield(S,required_fields{k})
        error(['Le fichier plink_phi_results.mat doit contenir %s. ', ...
               'Relancez d''abord plink_phi_mc_corrige.m.'], ...
               required_fields{k});
    end
end

N = double(S.N);
phi_vals = S.phi_vals(:).';
f_phi = S.f_phi(:).';
dphi_bins = S.dphi_bins(:).';
phi_bin_probability = S.phi_bin_probability(:).';
count_sat_by_bin = S.count_sat_by_bin(:).';
n_sim = double(S.n_sim);

%% ============================================================
%  2. Theorie
%% ============================================================

% Densite moyenne de satellites par radian de latitude
satellites_density_phi_th = N .* f_phi;

% Nombre moyen de satellites dans chaque tranche
satellites_per_bin_th = N .* phi_bin_probability;

% Verification globale
N_total_th = sum(satellites_per_bin_th);

%% ============================================================
%  3. Monte-Carlo
%% ============================================================

% count_sat_by_bin contient le nombre total d'observations cumulees
% sur les n_sim simulations.
satellites_per_bin_emp = count_sat_by_bin ./ n_sim;

% Densite empirique par radian
satellites_density_phi_emp = ...
    satellites_per_bin_emp ./ dphi_bins;

N_total_emp = sum(satellites_per_bin_emp);

%% ============================================================
%  4. Erreurs
%% ============================================================

valid = isfinite(satellites_density_phi_emp) ...
    & isfinite(satellites_density_phi_th);

rmse_satellites_density = sqrt(mean( ...
    (satellites_density_phi_emp(valid) ...
    -satellites_density_phi_th(valid)).^2));

relative_error_total = ...
    abs(N_total_emp-N_total_th) / max(abs(N_total_th),eps);

ratio_emp_th_phi = ...
    satellites_density_phi_emp ./ satellites_density_phi_th;

%% ============================================================
%  5. Figures
%% ============================================================

figure;
plot(rad2deg(phi_vals),satellites_density_phi_th,'LineWidth',2);
hold on;
plot(rad2deg(phi_vals),satellites_density_phi_emp,'--','LineWidth',1.8);
grid on;
xlabel('Latitude \phi (deg)');
ylabel('Densite moyenne de satellites par radian');
title('Nombre local de satellites : theorie et Monte-Carlo');
legend('Theorie','Monte-Carlo','Location','best');
hold off;

figure;
plot(rad2deg(phi_vals),satellites_per_bin_th,'LineWidth',2);
hold on;
plot(rad2deg(phi_vals),satellites_per_bin_emp,'--','LineWidth',1.8);
grid on;
xlabel('Latitude \phi (deg)');
ylabel('Nombre moyen de satellites dans la tranche');
title('Nombre moyen de satellites par tranche de latitude');
legend('Theorie','Monte-Carlo','Location','best');
hold off;

figure;
plot(rad2deg(phi_vals),ratio_emp_th_phi,'LineWidth',1.8);
hold on;
yline(1,'--','Accord parfait');
grid on;
xlabel('Latitude \phi (deg)');
ylabel('Empirique / theorie');
title('Qualite du modele local du nombre de satellites');
hold off;

%% ============================================================
%  6. Affichage
%% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' SATELLITES(phi) - THEORIE / MONTE-CARLO\n');
fprintf('============================================================\n');
fprintf('Nombre moyen tire N                 : %.8f\n',N);
fprintf('N total theorique                   : %.10f\n',N_total_th);
fprintf('N total empirique                   : %.10f\n',N_total_emp);
fprintf('Erreur relative totale             : %.3e\n',relative_error_total);
fprintf('RMSE densite locale                : %.3e\n', ...
    rmse_satellites_density);
fprintf('============================================================\n');

%% ============================================================
%  7. Sauvegarde
%% ============================================================

output_file = fullfile(script_dir,'N_phi_results.mat');

save(output_file, ...
    'N','phi_vals','f_phi','dphi_bins','phi_bin_probability', ...
    'satellites_density_phi_th','satellites_density_phi_emp', ...
    'satellites_per_bin_th','satellites_per_bin_emp', ...
    'N_total_th','N_total_emp', ...
    'ratio_emp_th_phi','rmse_satellites_density', ...
    'relative_error_total');

fprintf('Resultats sauvegardes dans %s\n',output_file);
