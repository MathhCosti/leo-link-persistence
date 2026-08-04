%% edges_phi.m
% Comparaison theorique / Monte-Carlo du nombre moyen d'aretes
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
    'dphi_bins','phi_bin_probability', ...
    'edges_per_bin_emp','edges_density_phi_emp','edges_total_emp'};

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
dphi_bins = S.dphi_bins(:).';
phi_bin_probability = S.phi_bin_probability(:).';

edges_per_bin_emp = S.edges_per_bin_emp(:).';
edges_density_phi_emp = S.edges_density_phi_emp(:).';
edges_total_emp = double(S.edges_total_emp);

%% ============================================================
%  2. Theorie avec N = moyenne des tirages
%% ============================================================

% Le nombre theorique d'aretes dans une tranche est calcule avec
% la probabilite exacte de cette tranche, et non avec
% f_phi(phi_centre)*Delta_phi. Cela traite correctement les
% singularites integrables de f_phi au voisinage de +/-inc.
edges_per_bin_th = ...
    N*(N-1)/2 .* phi_bin_probability .* p_link_phi;

edges_density_phi_th = ...
    edges_per_bin_th ./ dphi_bins;

E_total_th_from_phi = sum(edges_per_bin_th);

if isfield(S,'p_link_global_from_phi')
    p_link_global_th = double(S.p_link_global_from_phi);
else
    p_link_global_th = ...
        sum(phi_bin_probability.*p_link_phi);
end

E_total_th_direct = ...
    N*(N-1)/2*p_link_global_th;

%% ============================================================
%  3. Erreurs
%% ============================================================

valid_density = isfinite(edges_density_phi_emp) ...
    & isfinite(edges_density_phi_th);

rmse_edges_density = sqrt(mean( ...
    (edges_density_phi_emp(valid_density) ...
    -edges_density_phi_th(valid_density)).^2));

relative_error_total = ...
    abs(edges_total_emp-E_total_th_direct) ...
    / max(abs(E_total_th_direct),eps);

ratio_emp_th_phi = ...
    edges_density_phi_emp./edges_density_phi_th;

%% ============================================================
%  4. Figures
%% ============================================================

figure;
plot(rad2deg(phi_vals),edges_density_phi_th,'LineWidth',2);
hold on;
plot(rad2deg(phi_vals),edges_density_phi_emp,'--','LineWidth',1.8);
grid on;
xlabel('Latitude \phi (deg)');
ylabel('Densite moyenne d''aretes par radian');
title('Aretes locales : theorie et Monte-Carlo');
legend('Theorie','Monte-Carlo','Location','best');
hold off;

figure;
plot(rad2deg(phi_vals),edges_per_bin_th,'LineWidth',2);
hold on;
plot(rad2deg(phi_vals),edges_per_bin_emp,'--','LineWidth',1.8);
grid on;
xlabel('Latitude \phi (deg)');
ylabel('Nombre moyen d''aretes dans la tranche');
title('Nombre moyen d''aretes par tranche');
legend('Theorie','Monte-Carlo','Location','best');
hold off;

figure;
plot(rad2deg(phi_vals),ratio_emp_th_phi,'LineWidth',1.8);
hold on;
yline(1,'--','Accord parfait');
grid on;
xlabel('Latitude \phi (deg)');
ylabel('Empirique / theorie');
title('Qualite du modele local du nombre d''aretes');
hold off;

%% ============================================================
%  5. Affichage
%% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' ARETES(phi) - THEORIE / MONTE-CARLO\n');
fprintf('============================================================\n');
fprintf('Nombre moyen tire N                 : %.8f\n',N);
fprintf('E total theorique par integration  : %.10f\n', ...
    E_total_th_from_phi);
fprintf('E total theorique direct           : %.10f\n', ...
    E_total_th_direct);
fprintf('E total empirique                  : %.10f\n', ...
    edges_total_emp);
fprintf('Erreur relative totale            : %.3e\n', ...
    relative_error_total);
fprintf('RMSE densite locale               : %.3e\n', ...
    rmse_edges_density);
fprintf('============================================================\n');

%% ============================================================
%  6. Sauvegarde
%% ============================================================

output_file = fullfile(script_dir,'edges_phi_results.mat');

save(output_file, ...
    'N','phi_vals','f_phi','p_link_phi','dphi_bins', ...
    'phi_bin_probability', ...
    'edges_density_phi_th','edges_density_phi_emp', ...
    'edges_per_bin_th','edges_per_bin_emp', ...
    'E_total_th_from_phi','E_total_th_direct','edges_total_emp', ...
    'ratio_emp_th_phi','rmse_edges_density', ...
    'relative_error_total');

fprintf('Resultats sauvegardes dans %s\n',output_file);
