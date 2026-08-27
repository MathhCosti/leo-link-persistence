clear; clc; close all;

%% ============================================================
% T_ROUTE THEORIQUE - POSITIONS UTILISATEURS CONNUES
% WALKER DELTA - UNIFORMITE SPATIALE
%
% Sortie principale :
%   Troute_matrix(t,phi_A,phi_B)
%
% Les longitudes sont fixees et connues.
% La dependance temporelle vient de T_assign(t,phi).
%% ============================================================

%% Longitudes connues
user_lon_A_deg = 0;
user_lon_B_deg = 90;

lambda_A = deg2rad(user_lon_A_deg);
lambda_B = deg2rad(user_lon_B_deg);
Delta_lambda = lambda_A-lambda_B;

%% Fichiers
script_dir = fileparts(mfilename('fullpath'));

assignment_file = fullfile(script_dir,'T_assignation_results.mat');
path_file = fullfile(script_dir,'H_jumps_results.mat');

assert(isfile(assignment_file), ...
    'Fichier introuvable : %s',assignment_file);
assert(isfile(path_file), ...
    'Fichier introuvable : %s',path_file);

Adata = load(assignment_file);
Pdata = load(path_file);

%% Parametres assignation
assert(isfield(Adata,'user_lat_deg'),'user_lat_deg absent.');
assert(isfield(Adata,'time_phase'),'time_phase absent.');
assert(isfield(Adata,'MeanAssign_theory'),'MeanAssign_theory absent.');

user_lat_deg = double(Adata.user_lat_deg(:));
time_phase = double(Adata.time_phase(:));

Tassign_raw = double(Adata.MeanAssign_theory);

% Supporte :
%   [lambda,time,lat]  ou  [time,lat].
if ndims(Tassign_raw)==3
    if isfield(Adata,'lambda_values')
        lambda_values = double(Adata.lambda_values(:));
        lambda_target = 4e-7;
        [~,il] = min(abs(lambda_values-lambda_target));
    else
        il = 1;
        lambda_target = 4e-7;
    end
    Tassign_t_lat = squeeze(Tassign_raw(il,:,:));
else
    lambda_target = 4e-7;
    Tassign_t_lat = Tassign_raw;
end

Nt = numel(time_phase);
N_lat = numel(user_lat_deg);

assert(isequal(size(Tassign_t_lat),[Nt N_lat]), ...
    'Dimensions inattendues de T_assign.');

%% Parametres ISL
R_orbit = double(Pdata.R);
dmax = double(Pdata.dmax);
inc_deg = double(Pdata.inc_deg);
inc = deg2rad(inc_deg);

mu = 398600;
v_orb = sqrt(mu/R_orbit);
Delta_t_break = 10;

alpha_max = 2*asin(min(1,dmax/(2*R_orbit)));

%% Stockage
Troute_matrix = NaN(Nt,N_lat,N_lat);
Gamma_matrix_deg = NaN(N_lat,N_lat);
H_matrix = NaN(N_lat,N_lat);

BetaRouteISL_matrix = NaN(N_lat,N_lat);
MeanPbreakISL_matrix = NaN(N_lat,N_lat);

BetaRouteGSL_matrix = NaN(Nt,N_lat,N_lat);
BetaRouteTotal_matrix = NaN(Nt,N_lat,N_lat);

%% Geometrie constante
for qA = 1:N_lat
    phi_A = deg2rad(user_lat_deg(qA));

    for qB = 1:N_lat
        phi_B = deg2rad(user_lat_deg(qB));

        cos_gamma = ...
            sin(phi_A)*sin(phi_B) + ...
            cos(phi_A)*cos(phi_B)*cos(Delta_lambda);

        cos_gamma = max(-1,min(1,cos_gamma));
        gamma_AB = acos(cos_gamma);

        Gamma_matrix_deg(qA,qB) = rad2deg(gamma_AB);

        if gamma_AB < 1e-12
            H_links = 0;
        else
            H_links = ceil(gamma_AB/alpha_max);
        end

        H_matrix(qA,qB) = H_links;

        if H_links==0
            beta_route_ISL = 0;
            mean_p_break_ISL = 0;
        else
            phi_nodes = linspace(phi_A,phi_B,H_links+1);
            phi_links = 0.5*(phi_nodes(1:end-1)+phi_nodes(2:end));

            latitude_factor = ...
                sqrt(max(sin(inc)^2-sin(phi_links).^2,0)) ...
                ./ max(cos(phi_links),eps);

            p_break_links = ...
                (4*v_orb*Delta_t_break/(pi*dmax)) ...
                .* latitude_factor;

            p_break_links = min(max(p_break_links,0),1-eps);

            beta_links = -log1p(-p_break_links)/Delta_t_break;

            beta_route_ISL = sum(beta_links);
            mean_p_break_ISL = mean(p_break_links);
        end

        BetaRouteISL_matrix(qA,qB) = beta_route_ISL;
        MeanPbreakISL_matrix(qA,qB) = mean_p_break_ISL;
    end
end

%% Dependence temporelle via T_assign(t,phi)
for it = 1:Nt
    for qA = 1:N_lat
        Tassign_A = Tassign_t_lat(it,qA);
        if ~isfinite(Tassign_A) || Tassign_A<=0, continue; end

        beta_A = 1/Tassign_A;

        for qB = 1:N_lat
            Tassign_B = Tassign_t_lat(it,qB);
            if ~isfinite(Tassign_B) || Tassign_B<=0, continue; end

            beta_B = 1/Tassign_B;

            beta_GSL = beta_A+beta_B;
            beta_total = beta_GSL+BetaRouteISL_matrix(qA,qB);

            BetaRouteGSL_matrix(it,qA,qB) = beta_GSL;
            BetaRouteTotal_matrix(it,qA,qB) = beta_total;
            Troute_matrix(it,qA,qB) = 1/beta_total;
        end
    end
end

%% Moyenne temporelle locale
Troute_time_mean_matrix = squeeze(mean(Troute_matrix,1,'omitnan'));

%% Figures
snapshot_idx = unique(round(linspace(1,ceil(Nt/2),4)));

for s = 1:numel(snapshot_idx)
    it = snapshot_idx(s);

    figure;
    imagesc(user_lat_deg,user_lat_deg,squeeze(Troute_matrix(it,:,:)));
    set(gca,'YDir','normal');
    colorbar;
    xlabel('\phi_B (deg)');
    ylabel('\phi_A (deg)');
    title(sprintf( ...
        'T_{route}^{th}(t,\\phi_A,\\phi_B), t=%.0f s, \\Delta\\lambda=%.0f deg', ...
        time_phase(it),rad2deg(abs(Delta_lambda))));
end

figure;
imagesc(user_lat_deg,user_lat_deg,Troute_time_mean_matrix);
set(gca,'YDir','normal');
colorbar;
xlabel('\phi_B (deg)');
ylabel('\phi_A (deg)');
title('Moyenne temporelle de T_{route}^{th}');

%% Sauvegarde
save('T_route_th_positions_connues_results.mat', ...
    'Troute_matrix','Troute_time_mean_matrix', ...
    'Gamma_matrix_deg','H_matrix', ...
    'BetaRouteISL_matrix','BetaRouteGSL_matrix', ...
    'BetaRouteTotal_matrix','MeanPbreakISL_matrix', ...
    'time_phase','user_lat_deg', ...
    'user_lon_A_deg','user_lon_B_deg','Delta_lambda', ...
    'lambda_target','alpha_max','inc_deg','dmax','R_orbit', ...
    'v_orb','Delta_t_break');

fprintf('\nResultats sauvegardes dans :\n');
fprintf('  T_route_th_positions_connues_spatial_results.mat\n');
