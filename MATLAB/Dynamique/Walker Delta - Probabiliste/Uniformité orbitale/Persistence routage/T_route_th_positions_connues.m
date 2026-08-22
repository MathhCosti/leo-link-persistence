clear; clc; close all;

%% ============================================================
% MATRICE THEORIQUE T_route(phi_A,phi_B) POUR LONGITUDES FIXEES
%
% Les latitudes phi_A et phi_B varient comme auparavant.
% Les longitudes lambda_A et lambda_B sont supposees connues.
%
% Pour chaque couple :
%   cos(gamma_AB) =
%       sin(phi_A)sin(phi_B)
%       + cos(phi_A)cos(phi_B)cos(lambda_A-lambda_B)
%
% Puis :
%   H(phi_A,phi_B) = ceil(gamma_AB/alpha_max)
%
% avec :
%   alpha_max = 2 asin(dmax/(2R)).
%
% La valeur moyenne globale de H n'est donc plus utilisee.
%% ============================================================

%% Longitudes connues
user_lon_A_deg = 0;
user_lon_B_deg = 90;

lambda_A = deg2rad(user_lon_A_deg);
lambda_B = deg2rad(user_lon_B_deg);
Delta_lambda = lambda_A-lambda_B;

%% Source des durees d'assignation
assignment_source = 'theoretical';   % 'theoretical' ou 'empirical'

%% Fichiers
script_dir = fileparts(mfilename('fullpath'));

assignment_file = fullfile(script_dir,'T_assignation_results.mat');
shortest_path_file = fullfile(script_dir,'H_jumps_results.mat');

assert(isfile(assignment_file), ...
    'Fichier d''assignation introuvable : %s',assignment_file);

assert(isfile(shortest_path_file), ...
    'Fichier H_jumps introuvable : %s',shortest_path_file);

assignment_data = load(assignment_file);
path_data = load(shortest_path_file);

%% Durees d'assignation
assert(isfield(assignment_data,'lambda_values'), ...
    'Champ lambda_values manquant.');

assert(isfield(assignment_data,'user_lat_deg'), ...
    'Champ user_lat_deg manquant.');

lambda_values = double(assignment_data.lambda_values(:));
user_lat_deg = double(assignment_data.user_lat_deg(:));

switch lower(assignment_source)
    case 'theoretical'
        field_assign = 'MeanAssign_theory';
        assignment_label = 'theorique';
    case 'empirical'
        field_assign = 'MeanAssign_emp';
        assignment_label = 'empirique';
    otherwise
        error('assignment_source doit valoir theoretical ou empirical.');
end

assert(isfield(assignment_data,field_assign), ...
    'Le fichier d''assignation ne contient pas %s.',field_assign);

Tassign_matrix = double(assignment_data.(field_assign));

%% Parametres orbitaux
Delta_t_break = 10;
mu = 398600;

assert(isfield(path_data,'R'),'Le fichier H_jumps ne contient pas R.');
assert(isfield(path_data,'dmax'),'Le fichier H_jumps ne contient pas dmax.');
assert(isfield(path_data,'inc_deg'),'Le fichier H_jumps ne contient pas inc_deg.');

R_orbit = double(path_data.R);
dmax = double(path_data.dmax);
inc_deg = double(path_data.inc_deg);
inc = deg2rad(inc_deg);

v_orb = sqrt(mu/R_orbit);

alpha_max = 2*asin(min(1,dmax/(2*R_orbit)));

%% Densite utilisee pour T_assign
if isfield(path_data,'lambda')
    lambda_target = double(path_data.lambda);
else
    lambda_target = lambda_values(1);
end

[~,il] = min(abs(lambda_values-lambda_target));
Tassign_lat = Tassign_matrix(il,:).';

%% Calcul matriciel
N_lat = numel(user_lat_deg);

Troute_matrix = NaN(N_lat,N_lat);
Gamma_matrix_deg = NaN(N_lat,N_lat);
H_matrix = NaN(N_lat,N_lat);

BetaRouteISL_matrix = NaN(N_lat,N_lat);
BetaRouteGSL_matrix = NaN(N_lat,N_lat);
BetaRouteTotal_matrix = NaN(N_lat,N_lat);
MeanPbreakISL_matrix = NaN(N_lat,N_lat);

for qA = 1:N_lat

    phi_A = deg2rad(user_lat_deg(qA));
    Tassign_A = Tassign_lat(qA);

    if ~isfinite(Tassign_A) || Tassign_A <= 0
        continue;
    end

    beta_assign_A = 1/Tassign_A;

    for qB = 1:N_lat

        phi_B = deg2rad(user_lat_deg(qB));
        Tassign_B = Tassign_lat(qB);

        if ~isfinite(Tassign_B) || Tassign_B <= 0
            continue;
        end

        beta_assign_B = 1/Tassign_B;

        %% Separation angulaire entre utilisateurs
        cos_gamma = ...
            sin(phi_A)*sin(phi_B) + ...
            cos(phi_A)*cos(phi_B)*cos(Delta_lambda);

        cos_gamma = max(-1,min(1,cos_gamma));
        gamma_AB = acos(cos_gamma);

        Gamma_matrix_deg(qA,qB) = rad2deg(gamma_AB);

        %% Nombre de liens ISL
        if gamma_AB < 1e-12
            H_links = 0;
        else
            H_links = ceil(gamma_AB/alpha_max);
        end

        H_matrix(qA,qB) = H_links;

        %% Contribution ISL
        if H_links == 0

            beta_route_ISL = 0;
            mean_p_break_ISL = 0;

        else

            phi_nodes = linspace(phi_A,phi_B,H_links+1);
            phi_links = ...
                0.5*(phi_nodes(1:end-1)+phi_nodes(2:end));

            latitude_factor = ...
                sqrt(max(sin(inc)^2-sin(phi_links).^2,0)) ...
                ./ max(cos(phi_links),eps);

            p_break_ISL_links = ...
                (4*v_orb*Delta_t_break/(pi*dmax)) ...
                .* latitude_factor;

            p_break_ISL_links = ...
                min(p_break_ISL_links,1-eps);

            beta_break_ISL_links = ...
                -log1p(-p_break_ISL_links)/Delta_t_break;

            beta_route_ISL = sum(beta_break_ISL_links);
            mean_p_break_ISL = mean(p_break_ISL_links);
        end

        %% Contribution GSL et duree totale
        beta_route_GSL = beta_assign_A+beta_assign_B;
        beta_route_total = beta_route_GSL+beta_route_ISL;
        Troute = 1/beta_route_total;

        %% Stockage
        Troute_matrix(qA,qB) = Troute;
        BetaRouteISL_matrix(qA,qB) = beta_route_ISL;
        BetaRouteGSL_matrix(qA,qB) = beta_route_GSL;
        BetaRouteTotal_matrix(qA,qB) = beta_route_total;
        MeanPbreakISL_matrix(qA,qB) = mean_p_break_ISL;
    end
end

%% Affichage
fprintf('\n============================================================\n');
fprintf('MATRICE T_route POUR LONGITUDES FIXEES\n');
fprintf('============================================================\n');
fprintf('Longitude A        : %.2f deg\n',user_lon_A_deg);
fprintf('Longitude B        : %.2f deg\n',user_lon_B_deg);
fprintf('Delta longitude    : %.2f deg\n',rad2deg(abs(Delta_lambda)));
fprintf('Source assignation : %s\n',assignment_label);
fprintf('alpha_max          : %.4f deg\n',rad2deg(alpha_max));
fprintf('============================================================\n');

disp('Matrice H :');
disp(H_matrix);

disp('Matrice T_route (s) :');
disp(Troute_matrix);

%% Tables
lat_labels = matlab.lang.makeValidName( ...
    compose('phiB_%gdeg',user_lat_deg));

Troute_table = array2table( ...
    Troute_matrix,'VariableNames',lat_labels);

Troute_table = addvars( ...
    Troute_table,user_lat_deg, ...
    'Before',1,'NewVariableNames','phiA_deg');

H_table = array2table( ...
    H_matrix,'VariableNames',lat_labels);

H_table = addvars( ...
    H_table,user_lat_deg, ...
    'Before',1,'NewVariableNames','phiA_deg');

Gamma_table = array2table( ...
    Gamma_matrix_deg,'VariableNames',lat_labels);

Gamma_table = addvars( ...
    Gamma_table,user_lat_deg, ...
    'Before',1,'NewVariableNames','phiA_deg');

%% Carte T_route
figure;
imagesc(user_lat_deg,user_lat_deg,Troute_matrix);
set(gca,'YDir','normal');
colorbar;
xlabel('\phi_B (deg)');
ylabel('\phi_A (deg)');
title(sprintf('T_{route}(\\phi_A,\\phi_B), \\lambda_A=%.0f deg, \\lambda_B=%.0f deg', ...
    user_lon_A_deg,user_lon_B_deg));

for qA = 1:N_lat
    for qB = 1:N_lat
        if isfinite(Troute_matrix(qA,qB))
            text(user_lat_deg(qB),user_lat_deg(qA), ...
                sprintf('%.1f',Troute_matrix(qA,qB)), ...
                'HorizontalAlignment','center');
        end
    end
end

%% Carte H
figure;
imagesc(user_lat_deg,user_lat_deg,H_matrix);
set(gca,'YDir','normal');
colorbar;
xlabel('\phi_B (deg)');
ylabel('\phi_A (deg)');
title(sprintf( ...
    'H(\\phi_A,\\phi_B), \\Delta\\lambda=%.0f deg', ...
    rad2deg(abs(Delta_lambda))));

for qA = 1:N_lat
    for qB = 1:N_lat
        if isfinite(H_matrix(qA,qB))
            text(user_lat_deg(qB),user_lat_deg(qA), ...
                sprintf('%d',round(H_matrix(qA,qB))), ...
                'HorizontalAlignment','center');
        end
    end
end

%% Carte gamma
figure;
imagesc(user_lat_deg,user_lat_deg,Gamma_matrix_deg);
set(gca,'YDir','normal');
colorbar;
xlabel('\phi_B (deg)');
ylabel('\phi_A (deg)');
title(sprintf( ...
    '\\gamma_{AB}(\\phi_A,\\phi_B), \\Delta\\lambda=%.0f deg', ...
    rad2deg(abs(Delta_lambda))));

%% Verification de symetrie
symmetry_error = ...
    max(abs(Troute_matrix-Troute_matrix.'),[],'all','omitnan');

fprintf('\nEcart maximal de symetrie : %.6e s\n',symmetry_error);

%% Sauvegarde
save('T_route_th_positions_connues_results.mat', ...
    'Troute_matrix','Troute_table', ...
    'Gamma_matrix_deg','Gamma_table', ...
    'H_matrix','H_table', ...
    'BetaRouteISL_matrix', ...
    'BetaRouteGSL_matrix', ...
    'BetaRouteTotal_matrix', ...
    'MeanPbreakISL_matrix', ...
    'user_lat_deg','Tassign_lat', ...
    'user_lon_A_deg','user_lon_B_deg','Delta_lambda', ...
    'lambda_values','lambda_target','il', ...
    'alpha_max','inc_deg','dmax','R_orbit','v_orb', ...
    'Delta_t_break','assignment_source', ...
    'symmetry_error');

fprintf('\nResultats sauvegardes dans :\n');
fprintf('  T_route_th_positions_connues_results.mat\n');
