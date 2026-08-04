%% degree_phi.m
% Comparaison theorique / Monte-Carlo du degre moyen
% en fonction de la latitude.
%
% Entree :
%   plink_phi_results.mat produit par plink_phi.m
%
% La valeur N utilisee dans la theorie est la moyenne des N tires
% pendant les simulations Monte-Carlo.

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
    'N','phi_vals','f_phi','p_link_phi', ...
    'phi_bin_probability', ...
    'phi_bin_probability','degree_phi_emp','degree_global_emp'};

for k = 1:numel(required_fields)
    if ~isfield(S,required_fields{k})
        error(['Le fichier plink_phi_results.mat doit contenir %s. ', ...
               'Relancez d''abord plink_phi.m.'],required_fields{k});
    end
end

N = double(S.N);
phi_vals = S.phi_vals(:).';
f_phi = S.f_phi(:).';
p_link_phi = S.p_link_phi(:).';
phi_bin_probability = S.phi_bin_probability(:).';

degree_phi_emp = S.degree_phi_emp(:).';
degree_global_emp = double(S.degree_global_emp);

%% ============================================================
%  2. Theorie avec N = moyenne des tirages
%% ============================================================

degree_phi_th = (N-1).*p_link_phi;

degree_global_th_from_phi = ...
    sum(degree_phi_th.*phi_bin_probability);

if isfield(S,'p_link_global_from_phi')
    p_link_global_th = double(S.p_link_global_from_phi);
else
    p_link_global_th = ...
        sum(p_link_phi.*phi_bin_probability);
end

degree_global_th_direct = ...
    (N-1)*p_link_global_th;

%% ============================================================
%  3. Erreurs
%% ============================================================

valid = isfinite(degree_phi_emp) & isfinite(degree_phi_th);

rmse_degree_phi = sqrt(mean( ...
    (degree_phi_emp(valid)-degree_phi_th(valid)).^2));

relative_error_global = ...
    abs(degree_global_emp-degree_global_th_direct) ...
    / max(abs(degree_global_th_direct),eps);

ratio_emp_th_phi = degree_phi_emp./degree_phi_th;

%% ============================================================
%  4. Figures
%% ============================================================

figure;
plot(rad2deg(phi_vals),degree_phi_th,'LineWidth',2);
hold on;
plot(rad2deg(phi_vals),degree_phi_emp,'--','LineWidth',1.8);
grid on;
xlabel('Latitude \phi (deg)');
ylabel('\mu^\Delta(\phi)');
title('Degre moyen local : theorie et Monte-Carlo');
legend('Theorie','Monte-Carlo','Location','best');
hold off;

figure;
plot(rad2deg(phi_vals),ratio_emp_th_phi,'LineWidth',1.8);
hold on;
yline(1,'--','Accord parfait');
grid on;
xlabel('Latitude \phi (deg)');
ylabel('Empirique / theorie');
title('Qualite du modele local du degre moyen');
hold off;

%% ============================================================
%  5. Affichage
%% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' DEGRE(phi) - THEORIE / MONTE-CARLO\n');
fprintf('============================================================\n');
fprintf('Nombre moyen tire N                 : %.8f\n',N);
fprintf('Degre global theorique integre     : %.10f\n', ...
    degree_global_th_from_phi);
fprintf('Degre global theorique direct      : %.10f\n', ...
    degree_global_th_direct);
fprintf('Degre global empirique             : %.10f\n', ...
    degree_global_emp);
fprintf('Erreur relative globale            : %.3e\n', ...
    relative_error_global);
fprintf('RMSE locale                        : %.3e\n', ...
    rmse_degree_phi);
fprintf('============================================================\n');

%% ============================================================
%  6. Sauvegarde
%% ============================================================

output_file = fullfile(script_dir,'degree_phi_results.mat');

save(output_file, ...
    'N','phi_vals','f_phi','p_link_phi', ...
    'degree_phi_th','degree_phi_emp', ...
    'degree_global_th_from_phi','degree_global_th_direct', ...
    'degree_global_emp','ratio_emp_th_phi', ...
    'rmse_degree_phi','relative_error_global');

fprintf('Resultats sauvegardes dans %s\n',output_file);
