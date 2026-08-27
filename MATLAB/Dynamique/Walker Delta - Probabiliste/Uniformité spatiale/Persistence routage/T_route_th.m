clear; clc; close all;

%% ============================================================
% DUREE MOYENNE THEORIQUE DE ROUTE PAR QUADRATURE
% WALKER DELTA - UNIFORMITE SPATIALE
%
% Integration sur :
%   - phase initiale u0,A ~ f_U(u,0)=|cos(u)|/4 ;
%   - phase initiale u0,B ~ f_U(u,0)=|cos(u)|/4 ;
%   - temps t uniforme sur une demi-periode spatiale ;
%   - Delta_lambda uniforme sur [0,pi].
%
% A l'instant t :
%   u(t)=u0+omega_sat*t
%   phi(t)=asin(sin(i) sin(u(t)))
%
% T_assign = T_assign(t,phi)
%
% Conditionnement analytique spatial :
%   P_conn(t,q) ~= prod_k p_step(t,phi_k)
%
% avec
%   p_step(t,phi)=1-exp[-lambda(t,phi) A_progress].
%% ============================================================

%% Reglages
Nu = 24;
Ntq = 16;
Ndlon = 36;

Delta_t_break = 10;
mu = 398600;

%% Fichiers
script_dir = fileparts(mfilename('fullpath'));

assignment_file = fullfile(script_dir,'T_assignation_results.mat');
path_file = fullfile(script_dir,'H_jumps_results.mat');
emp_route_file = fullfile(script_dir,'T_route_emp_results.mat');

assert(isfile(assignment_file), ...
    'Fichier introuvable : %s',assignment_file);
assert(isfile(path_file), ...
    'Fichier introuvable : %s',path_file);

Adata = load(assignment_file);
Pdata = load(path_file);

has_emp_route = isfile(emp_route_file);
if has_emp_route
    Edata = load(emp_route_file);
else
    Edata = struct();
end

%% Parametres orbitaux
R_orbit = double(Pdata.R);
dmax = double(Pdata.dmax);
inc_deg = double(Pdata.inc_deg);
inc = deg2rad(inc_deg);

omega_sat = sqrt(mu/R_orbit^3);
T_orbit = 2*pi/omega_sat;
T_spatial = T_orbit/2;

v_orb = sqrt(mu/R_orbit);

alpha_max = 2*asin(min(1,dmax/(2*R_orbit)));
Hmax = ceil(pi/alpha_max);

%% Densite
lambda = 4e-7;
if isfield(Pdata,'lambda')
    lambda = double(Pdata.lambda);
end

N_mean = lambda*4*pi*R_orbit^2;

A_cap_ISL = 2*pi*R_orbit^2*(1-cos(alpha_max));
A_progress = 0.5*A_cap_ISL;

%% T_assign(t,phi)
assert(isfield(Adata,'user_lat_deg'),'user_lat_deg absent.');
assert(isfield(Adata,'time_phase'),'time_phase absent.');
assert(isfield(Adata,'MeanAssign_theory'),'MeanAssign_theory absent.');

assign_lat_deg = double(Adata.user_lat_deg(:));
assign_time = double(Adata.time_phase(:));
Tassign_raw = double(Adata.MeanAssign_theory);

if ndims(Tassign_raw)==3
    if isfield(Adata,'lambda_values')
        lambda_values = double(Adata.lambda_values(:));
        [~,il] = min(abs(lambda_values-lambda));
    else
        il = 1;
    end
    Tassign_grid = squeeze(Tassign_raw(il,:,:));
else
    Tassign_grid = Tassign_raw;
end

assert(isequal(size(Tassign_grid), ...
    [numel(assign_time),numel(assign_lat_deg)]), ...
    'Dimensions inattendues de MeanAssign_theory.');

%% Repli de T_assign sur la demi-periode fondamentale
%
% T_assignation_spatial peut avoir ete sauvegarde sur une periode orbitale
% complete, alors que la densite spatiale verifie une periodicite T_orbit/2.
% On replie donc tous les instants modulo T_spatial, puis on moyenne les
% valeurs correspondant a une meme phase temporelle.

assign_time_fold = mod(assign_time,T_spatial);

% Tri croissant obligatoire pour interp2 / griddedInterpolant.
[assign_time_fold_sorted,idx_time_sort] = sort(assign_time_fold);
Tassign_fold_sorted = Tassign_grid(idx_time_sort,:);

% Regroupement des phases identiques ou numeriquement quasi identiques.
tol_time = 1e-8*T_spatial;
time_unique = [];
Tassign_unique = [];

k0 = 1;
while k0 <= numel(assign_time_fold_sorted)

    k1 = k0;

    while k1 < numel(assign_time_fold_sorted) && ...
            abs(assign_time_fold_sorted(k1+1) - ...
                assign_time_fold_sorted(k0)) <= tol_time
        k1 = k1 + 1;
    end

    time_unique(end+1,1) = ...
        mean(assign_time_fold_sorted(k0:k1)); %#ok<SAGROW>

    Tassign_unique(end+1,:) = ...
        mean(Tassign_fold_sorted(k0:k1,:),1,'omitnan'); %#ok<SAGROW>

    k0 = k1 + 1;
end

% Securite : ordre strictement croissant.
[time_unique,idx_unique] = sort(time_unique);
Tassign_unique = Tassign_unique(idx_unique,:);

% Extension periodique autour de [0,T_spatial).
time_ext = [ ...
    time_unique(end)-T_spatial; ...
    time_unique; ...
    time_unique(1)+T_spatial];

Tassign_ext = [ ...
    Tassign_unique(end,:); ...
    Tassign_unique; ...
    Tassign_unique(1,:)];

assert(all(diff(time_ext) > 0), ...
    'La grille temporelle de T_assign n''est pas strictement croissante.');

fprintf('\nInterpolation T_assign spatiale :\n');
fprintf('  Points initiaux       : %d\n',numel(assign_time));
fprintf('  Points apres repli    : %d\n',numel(time_unique));
fprintf('  Intervalle fondamental: [0, %.3f] s\n',T_spatial);

%% Donnees empiriques facultatives pour diagnostics
has_emp_conn_tg = isfield(Pdata,'GammaStatisticsTime');
has_emp_H_tg = has_emp_route && isfield(Edata,'HRouteStatisticsTimeGamma');

beta_link_emp = NaN;
Troute_emp_reference = NaN;

if has_emp_route
    if isfield(Edata,'beta_link_emp_direct')
        beta_link_emp = double(Edata.beta_link_emp_direct);
    end
    if isfield(Edata,'Troute_emp_mean')
        Troute_emp_reference = double(Edata.Troute_emp_mean);
    end
end

%% Quadratures
[xu,wu0] = gauss_legendre(Nu,0,2*pi);

% Mesure f_U(u,0)du.
wu = wu0 .* abs(cos(xu))/4;
wu = wu/sum(wu);

[xt,wt] = gauss_legendre(Ntq,0,T_spatial);
wt = wt/T_spatial;

[xd,wd] = gauss_legendre(Ndlon,0,pi);
wd = wd/pi;

%% Accumulateurs
P_H = zeros(Hmax,1);

weight_sum = 0;
Troute_sum = 0;

weight_cond_ana = 0;
Troute_cond_ana_sum = 0;

H_sum = 0;
Gamma_sum = 0;
H_cond_ana_sum = 0;
Gamma_cond_ana_sum = 0;

BetaGSL_sum = 0;
BetaISL_sum = 0;
BetaGSL_cond_ana_sum = 0;
BetaISL_cond_ana_sum = 0;

% Conditionnement empirique facultatif
weight_cond_emp = 0;
Troute_cond_emp_sum = 0;
Troute_cond_emp_Hemp_sum = 0;
Troute_cond_emp_Hemp_PbreakEmp_sum = 0;
Hemp_cond_emp_sum = 0;

%% Integration
for it = 1:Ntq

    t = xt(it);

    for ia = 1:Nu

        uA_t = mod(xu(ia)+omega_sat*t,2*pi);
        phi_A = asin(sin(inc)*sin(uA_t));
        phi_A_deg = rad2deg(phi_A);

        Tassign_A = interp_Tassign( ...
            t,phi_A_deg,time_ext,assign_lat_deg,Tassign_ext,T_spatial);

        if ~isfinite(Tassign_A) || Tassign_A<=0
            continue;
        end

        beta_A = 1/Tassign_A;

        for ib = 1:Nu

            uB_t = mod(xu(ib)+omega_sat*t,2*pi);
            phi_B = asin(sin(inc)*sin(uB_t));
            phi_B_deg = rad2deg(phi_B);

            Tassign_B = interp_Tassign( ...
                t,phi_B_deg,time_ext,assign_lat_deg,Tassign_ext,T_spatial);

            if ~isfinite(Tassign_B) || Tassign_B<=0
                continue;
            end

            beta_B = 1/Tassign_B;
            beta_GSL = beta_A+beta_B;

            w_base = wt(it)*wu(ia)*wu(ib);

            for id = 1:Ndlon

                dl = xd(id);

                %% Separation
                cg = ...
                    sin(phi_A)*sin(phi_B) + ...
                    cos(phi_A)*cos(phi_B)*cos(dl);

                cg = max(-1,min(1,cg));
                gamma = acos(cg);
                gamma_deg = rad2deg(gamma);

                H = max(1,ceil(gamma/alpha_max));
                H = min(H,Hmax);

                %% Latitudes des liens
                phi_nodes = linspace(phi_A,phi_B,H+1);
                phi_links = ...
                    0.5*(phi_nodes(1:end-1)+phi_nodes(2:end));

                %% Rupture ISL
                latitude_factor = ...
                    sqrt(max(sin(inc)^2-sin(phi_links).^2,0)) ...
                    ./ max(cos(phi_links),eps);

                p_break_links = ...
                    (4*v_orb*Delta_t_break/(pi*dmax)) ...
                    .* latitude_factor;

                p_break_links = ...
                    min(max(p_break_links,0),1-eps);

                beta_links = ...
                    -log1p(-p_break_links)/Delta_t_break;

                beta_ISL = sum(beta_links);
                beta_link_th_local = beta_ISL/max(H,1);

                Troute = 1/(beta_GSL+beta_ISL);

                %% Conditionnement analytique spatial
                %
                % H sauts -> H-1 relais intermediaires.
                if H<=1
                    Pconn_ana = 1;
                else
                    phi_relays = phi_nodes(2:end-1);

                    p_step = zeros(size(phi_relays));

                    for jr = 1:numel(phi_relays)
                        lambda_local = spatial_surface_density( ...
                            t,phi_relays(jr), ...
                            N_mean,R_orbit,inc,omega_sat);

                        p_step(jr) = ...
                            1-exp(-lambda_local*A_progress);
                    end

                    Pconn_ana = prod(p_step);
                end

                Pconn_ana = min(max(Pconn_ana,0),1);

                %% Poids
                w = w_base*wd(id);
                wc_ana = w*Pconn_ana;

                weight_sum = weight_sum+w;
                Troute_sum = Troute_sum+w*Troute;

                P_H(H) = P_H(H)+w;
                H_sum = H_sum+w*H;
                Gamma_sum = Gamma_sum+w*gamma;

                BetaGSL_sum = BetaGSL_sum+w*beta_GSL;
                BetaISL_sum = BetaISL_sum+w*beta_ISL;

                weight_cond_ana = weight_cond_ana+wc_ana;
                Troute_cond_ana_sum = ...
                    Troute_cond_ana_sum+wc_ana*Troute;

                H_cond_ana_sum = ...
                    H_cond_ana_sum+wc_ana*H;

                Gamma_cond_ana_sum = ...
                    Gamma_cond_ana_sum+wc_ana*gamma;

                BetaGSL_cond_ana_sum = ...
                    BetaGSL_cond_ana_sum+wc_ana*beta_GSL;

                BetaISL_cond_ana_sum = ...
                    BetaISL_cond_ana_sum+wc_ana*beta_ISL;

                %% Diagnostics avec conditionnement empirique
                if has_emp_conn_tg
                    Pconn_emp = lookup_time_gamma( ...
                        Pdata.GammaStatisticsTime, ...
                        t,gamma_deg,'ConnectionProbability',T_spatial);

                    if isfinite(Pconn_emp)
                        Pconn_emp = min(max(Pconn_emp,0),1);
                        wc_emp = w*Pconn_emp;

                        weight_cond_emp = weight_cond_emp+wc_emp;
                        Troute_cond_emp_sum = ...
                            Troute_cond_emp_sum+wc_emp*Troute;

                        if has_emp_H_tg
                            H_emp = lookup_time_gamma( ...
                                Edata.HRouteStatisticsTimeGamma, ...
                                t,gamma_deg,'MeanHopCount',T_spatial);

                            if isfinite(H_emp)
                                H_emp = max(H_emp,0);

                                T_Hemp = ...
                                    1/(beta_GSL+H_emp*beta_link_th_local);

                                Troute_cond_emp_Hemp_sum = ...
                                    Troute_cond_emp_Hemp_sum+wc_emp*T_Hemp;

                                Hemp_cond_emp_sum = ...
                                    Hemp_cond_emp_sum+wc_emp*H_emp;

                                if isfinite(beta_link_emp) && beta_link_emp>=0
                                    T_Hemp_Pemp = ...
                                        1/(beta_GSL+H_emp*beta_link_emp);

                                    Troute_cond_emp_Hemp_PbreakEmp_sum = ...
                                        Troute_cond_emp_Hemp_PbreakEmp_sum + ...
                                        wc_emp*T_Hemp_Pemp;
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    fprintf('Quadrature temps %d/%d terminee\n',it,Ntq);
end

%% Resultats
P_H = P_H/weight_sum;

Troute_mean = Troute_sum/weight_sum;
H_mean_from_distribution = H_sum/weight_sum;
Gamma_mean_deg = rad2deg(Gamma_sum/weight_sum);

MeanBetaRouteGSL = BetaGSL_sum/weight_sum;
MeanBetaRouteISL = BetaISL_sum/weight_sum;

assert(weight_cond_ana>0,'Poids analytique conditionne nul.');

P_connected_analytic = weight_cond_ana/weight_sum;
Troute_mean_cond_ana = Troute_cond_ana_sum/weight_cond_ana;

H_mean_cond_ana = H_cond_ana_sum/weight_cond_ana;
Gamma_mean_cond_ana_deg = ...
    rad2deg(Gamma_cond_ana_sum/weight_cond_ana);

MeanBetaRouteGSL_cond_ana = ...
    BetaGSL_cond_ana_sum/weight_cond_ana;

MeanBetaRouteISL_cond_ana = ...
    BetaISL_cond_ana_sum/weight_cond_ana;

%% Empirique facultatif
P_connected_empirical = NaN;
Troute_mean_cond_emp = NaN;
Troute_mean_cond_emp_HgammaEmp = NaN;
Troute_mean_cond_emp_HgammaEmp_PbreakEmp = NaN;
Hgamma_emp_mean_under_emp_conditioning = NaN;

if weight_cond_emp>0
    P_connected_empirical = weight_cond_emp/weight_sum;
    Troute_mean_cond_emp = Troute_cond_emp_sum/weight_cond_emp;

    if has_emp_H_tg
        Troute_mean_cond_emp_HgammaEmp = ...
            Troute_cond_emp_Hemp_sum/weight_cond_emp;

        Hgamma_emp_mean_under_emp_conditioning = ...
            Hemp_cond_emp_sum/weight_cond_emp;

        if isfinite(beta_link_emp)
            Troute_mean_cond_emp_HgammaEmp_PbreakEmp = ...
                Troute_cond_emp_Hemp_PbreakEmp_sum/weight_cond_emp;
        end
    end
end

%% Loi de H
h_values = (1:Hmax).';
cdf_H = cumsum(P_H);

HopDistribution = table(h_values,P_H,cdf_H, ...
    'VariableNames',{'HopCount','Probability','CDF'});

%% Affichage
fprintf('\n============================================================\n');
fprintf('T_ROUTE - QUADRATURE DELTA UNIFORMITE SPATIALE\n');
fprintf('============================================================\n');
fprintf('Nu / Nt / Ndlon                    : %d / %d / %d\n', ...
    Nu,Ntq,Ndlon);
fprintf('lambda                             : %.4e km^-2\n',lambda);
fprintf('T_spatial                          : %.3f s\n',T_spatial);
fprintf('E[H] brute                         : %.6f\n',H_mean_from_distribution);
fprintf('E[gamma] brute                     : %.6f deg\n',Gamma_mean_deg);
fprintf('E[T_route] brute                   : %.6f s\n',Troute_mean);

fprintf('\n--- CONDITIONNEMENT ANALYTIQUE SPATIAL ---\n');
fprintf('P(route existe)                    : %.6f\n', ...
    P_connected_analytic);
fprintf('E[T_route | route]                 : %.6f s\n', ...
    Troute_mean_cond_ana);
fprintf('E[H | route]                       : %.6f\n', ...
    H_mean_cond_ana);
fprintf('E[gamma | route]                   : %.6f deg\n', ...
    Gamma_mean_cond_ana_deg);
fprintf('E[beta_GSL | route]                : %.10e s^-1\n', ...
    MeanBetaRouteGSL_cond_ana);
fprintf('E[beta_ISL | route]                : %.10e s^-1\n', ...
    MeanBetaRouteISL_cond_ana);

if weight_cond_emp>0
    fprintf('\n--- COMPARAISON DES APPROXIMATIONS ---\n');
    fprintf('1) Cond.emp + H_th + pbreak_th      : %.6f s\n', ...
        Troute_mean_cond_emp);

    if isfinite(Troute_mean_cond_emp_HgammaEmp)
        fprintf('2) Cond.emp + H_emp(t,gamma)+pbreak_th : %.6f s\n', ...
            Troute_mean_cond_emp_HgammaEmp);
    end

    if isfinite(Troute_mean_cond_emp_HgammaEmp_PbreakEmp)
        fprintf('3) Cond.emp + H_emp(t,gamma)+pbreak_emp: %.6f s\n', ...
            Troute_mean_cond_emp_HgammaEmp_PbreakEmp);
    end

    if isfinite(Troute_emp_reference)
        fprintf('Reference simulation complete         : %.6f s\n', ...
            Troute_emp_reference);
    end
end

fprintf('============================================================\n');

%% Figures
figure;
bar(h_values,P_H);
grid on;
xlabel('Nombre de sauts H');
ylabel('P(H=h)');
title('Distribution theorique de H - Delta spatial');

if weight_cond_emp>0
    values = Troute_mean_cond_emp;
    labels = {'Cond.emp + H_{th} + pbreak_{th}'};

    if isfinite(Troute_mean_cond_emp_HgammaEmp)
        values(end+1)=Troute_mean_cond_emp_HgammaEmp;
        labels{end+1}='Cond.emp + H_{emp}(t,\gamma) + pbreak_{th}';
    end

    if isfinite(Troute_mean_cond_emp_HgammaEmp_PbreakEmp)
        values(end+1)=Troute_mean_cond_emp_HgammaEmp_PbreakEmp;
        labels{end+1}='Cond.emp + H_{emp}(t,\gamma) + pbreak_{emp}';
    end

    if isfinite(Troute_emp_reference)
        values(end+1)=Troute_emp_reference;
        labels{end+1}='Simulation empirique';
    end

    figure;
    bar(values);
    set(gca,'XTick',1:numel(labels),'XTickLabel',labels);
    xtickangle(18);
    ylabel('T_{route} moyen (s)');
    title('Impact successif des approximations - Delta spatial');
    grid on;
end

%% Table globale
Results = table( ...
    lambda,inc_deg,dmax,R_orbit,Nu,Ntq,Ndlon, ...
    H_mean_from_distribution,Gamma_mean_deg,Troute_mean, ...
    P_connected_analytic,Troute_mean_cond_ana, ...
    H_mean_cond_ana,Gamma_mean_cond_ana_deg, ...
    MeanBetaRouteGSL_cond_ana,MeanBetaRouteISL_cond_ana, ...
    P_connected_empirical,Troute_mean_cond_emp, ...
    Troute_mean_cond_emp_HgammaEmp, ...
    Troute_mean_cond_emp_HgammaEmp_PbreakEmp, ...
    Hgamma_emp_mean_under_emp_conditioning,beta_link_emp, ...
    'VariableNames',{ ...
    'Lambda','Inclination_deg','Dmax_km','OrbitalRadius_km', ...
    'QuadratureOrderU','QuadratureOrderTime','QuadratureOrderDeltaLongitude', ...
    'MeanHopCountRaw','MeanAngularSeparationRaw_deg','MeanRouteLifetimeRaw_s', ...
    'AnalyticConnectionProbability','MeanRouteLifetimeAnalyticCond_s', ...
    'MeanHopCountAnalyticCond','MeanAngularSeparationAnalyticCond_deg', ...
    'MeanBetaGSLAnalyticCond_per_s','MeanBetaISLAnalyticCond_per_s', ...
    'EmpiricalConnectionProbability','MeanRouteLifetimeEmpiricalCond_Hth_s', ...
    'MeanRouteLifetimeEmpiricalCond_Hemp_s', ...
    'MeanRouteLifetimeEmpiricalCond_Hemp_PbreakEmp_s', ...
    'MeanEmpiricalHopCountUnderEmpiricalConditioning', ...
    'EmpiricalBetaLink_per_s'});

disp(Results);

%% Sauvegarde
save('T_route_th_results.mat', ...
    'Results','HopDistribution','h_values','P_H','cdf_H', ...
    'Troute_mean','Troute_mean_cond_ana', ...
    'P_connected_analytic','H_mean_cond_ana', ...
    'Gamma_mean_cond_ana_deg', ...
    'MeanBetaRouteGSL','MeanBetaRouteISL', ...
    'MeanBetaRouteGSL_cond_ana','MeanBetaRouteISL_cond_ana', ...
    'P_connected_empirical','Troute_mean_cond_emp', ...
    'Troute_mean_cond_emp_HgammaEmp', ...
    'Troute_mean_cond_emp_HgammaEmp_PbreakEmp', ...
    'Hgamma_emp_mean_under_emp_conditioning', ...
    'beta_link_emp','Troute_emp_reference', ...
    'Nu','Ntq','Ndlon','lambda','inc_deg','dmax','R_orbit', ...
    'alpha_max','A_cap_ISL','A_progress', ...
    'T_orbit','T_spatial','Delta_t_break');

%% ============================================================
% Fonctions locales
%% ============================================================

function T = interp_Tassign( ...
    t,phi_deg,time_ext,lat_deg,T_ext,Tsp)

    tq = mod(t,Tsp);

    phi_q = min(max(phi_deg,min(lat_deg)),max(lat_deg));

    T = interp2( ...
        lat_deg,time_ext,T_ext, ...
        phi_q,tq,'linear');

    if ~isfinite(T)
        T = interp2( ...
            lat_deg,time_ext,T_ext, ...
            phi_q,tq,'nearest');
    end
end

function rho = spatial_surface_density( ...
    t,phi,Nmean,R,inc,omega)

    % Densite surfacique sur la sphere orbitale :
    %
    % lambda(t,phi)
    %   = N f_phi(t,phi)/(2*pi*R^2 cos(phi)).

    if abs(phi) >= inc-1e-12
        % Limite tres dense au bord ; p_step saturera vers 1.
        rho = realmax('double')^0.25;
        return;
    end

    x = max(-1,min(1,sin(phi)/sin(inc)));

    up = asin(x);
    um = pi-up;

    fUp = abs(cos(up-omega*t))/4;
    fUm = abs(cos(um-omega*t))/4;

    denom = sqrt(max(sin(inc)^2-sin(phi)^2,eps));

    fphi = cos(phi)/denom*(fUp+fUm);

    rho = ...
        Nmean*fphi/(2*pi*R^2*max(cos(phi),eps));

    rho = max(rho,0);
end

function value = lookup_time_gamma( ...
    T,t,gamma_deg,field_name,Tsp)

    times = unique(double(T.Time_s(:)));
    tq = mod(t,Tsp);

    % Distance periodique au temps tabule.
    d = abs(times-tq);
    d = min(d,Tsp-d);

    [~,it] = min(d);
    tsel = times(it);

    mask = abs(double(T.Time_s)-tsel)<1e-9;

    g = double(T.GammaCenter_deg(mask));
    y = double(T.(field_name)(mask));

    valid = isfinite(g)&isfinite(y);

    g = g(valid);
    y = y(valid);

    if numel(g)<2
        value = NaN;
        return;
    end

    [g,idx] = sort(g);
    y = y(idx);

    gq = min(max(gamma_deg,min(g)),max(g));

    value = interp1(g,y,gq,'linear');
end

function [x,w] = gauss_legendre(n,a,b)

    k = (1:n-1).';
    beta = k./sqrt(4*k.^2-1);

    J = diag(beta,1)+diag(beta,-1);

    [V,D] = eig(J);

    x0 = diag(D);
    [x0,idx] = sort(x0);
    V = V(:,idx);

    w0 = 2*(V(1,:).^2).';

    x = (a+b)/2+(b-a)/2*x0;
    w = (b-a)/2*w0;
end
