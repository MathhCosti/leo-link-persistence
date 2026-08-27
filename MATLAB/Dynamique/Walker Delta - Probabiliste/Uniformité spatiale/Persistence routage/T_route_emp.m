clear; clc; close all;

%% ============================================================
% DUREE DE VIE EMPIRIQUE D'UNE ROUTE
% WALKER DELTA - UNIFORMITE SPATIALE
%
% Les phases initiales suivent :
%   f_U(u,0)=|cos(u)|/4.
%
% Les episodes sont en plus etiquetes par leur phase temporelle
% modulo T_spatial=T_orbit/2.
%% ============================================================

rng(8);

%% Parametres physiques
R_earth = 6371;
h = 550;
R = R_earth+h;
mu = 398600;

omega_sat = sqrt(mu/R^3);
omega_earth = 2*pi/86164;
T_orbit = 2*pi/omega_sat;
T_spatial = T_orbit/2;

%% Constellation
lambda = 4e-7;
inc_deg = 58;
inc = deg2rad(inc_deg);

%% Liens
dmax = 1500;
elevation_min_deg = 20;
elevation_min = deg2rad(elevation_min_deg);

%% Simulation
N_realizations = 50;
dt = 2;
Tmax = 15000;
time_values = (0:dt:Tmax-dt).';
Nt = numel(time_values);

restrict_users_to_orbital_band = true;
min_episode_samples = 1;

%% Stockage
route_durations = [];
route_hops_initial = [];
route_gamma_ground_deg = [];
route_end_cause = strings(0,1);
route_realization = [];
route_user_lat_A_deg = [];
route_user_lat_B_deg = [];
route_start_phase_s = [];

total_samples = N_realizations*Nt;
samples_both_assigned = 0;
samples_route_exists = 0;

direct_link_exposure = 0;
direct_ISL_break_count = 0;
direct_route_steps_with_ISL_break = 0;

%% Monte-Carlo
for r = 1:N_realizations

    N_mean = lambda*4*pi*R^2;
    N_sat = poissrnd(N_mean);

    if N_sat < 2, continue; end

    Omega = 2*pi*rand(N_sat,1);
    u0 = sample_u_spatial(N_sat);

    %% Utilisateurs fixes sur la Terre
    if restrict_users_to_orbital_band
        sin_lat = (2*rand(2,1)-1)*sin(inc);
        user_lat = asin(sin_lat);
    else
        user_lat = asin(2*rand(2,1)-1);
    end

    user_lon0 = 2*pi*rand(2,1)-pi;

    cos_gamma_ground = ...
        sin(user_lat(1))*sin(user_lat(2)) + ...
        cos(user_lat(1))*cos(user_lat(2))* ...
        cos(user_lon0(1)-user_lon0(2));

    cos_gamma_ground = max(-1,min(1,cos_gamma_ground));
    gamma_ground_deg = rad2deg(acos(cos_gamma_ground));

    %% Route active
    route_active = false;
    route_start_time = NaN;
    route_start_phase = NaN;
    route_path = [];
    access_A_initial = 0;
    access_B_initial = 0;
    route_hops = NaN;

    for k = 1:Nt
        t = time_values(k);

        %% Satellites
        u_t = mod(u0+omega_sat*t,2*pi);
        sat_pos = walker_delta_positions(R,inc,Omega,u_t);

        %% Graphe ISL
        D = squareform(pdist(sat_pos));
        A = (D<=dmax)&(D>0);
        G = graph(A);

        %% Utilisateurs
        lon_t = mod(user_lon0+omega_earth*t+pi,2*pi)-pi;
        user_pos = ground_user_positions(R_earth,user_lat,lon_t);

        %% Assignations
        current_assignment = zeros(2,1);

        for q=1:2
            current_assignment(q) = ...
                closest_visible_satellite( ...
                sat_pos,user_pos(q,:),R_earth,elevation_min);
        end

        both_assigned = all(current_assignment>0);
        samples_both_assigned = samples_both_assigned+both_assigned;

        %% Persistance route active
        if route_active

            H_current = max(0,numel(route_path)-1);

            if H_current>0
                direct_link_exposure = ...
                    direct_link_exposure+H_current*dt;

                idx_links = sub2ind(size(A), ...
                    route_path(1:end-1),route_path(2:end));

                n_broken = nnz(~A(idx_links));

                if n_broken>0
                    direct_ISL_break_count = ...
                        direct_ISL_break_count+n_broken;
                    direct_route_steps_with_ISL_break = ...
                        direct_route_steps_with_ISL_break+1;
                end
            end

            cause = "";

            if current_assignment(1) ~= access_A_initial
                cause = "assignment_A";
            elseif current_assignment(2) ~= access_B_initial
                cause = "assignment_B";
            elseif ~fixed_path_is_alive(A,route_path)
                cause = "ISL_break";
            end

            if strlength(cause)>0
                duration = t-route_start_time;

                if duration>=min_episode_samples*dt
                    route_durations(end+1,1)=duration; %#ok<SAGROW>
                    route_hops_initial(end+1,1)=route_hops; %#ok<SAGROW>
                    route_gamma_ground_deg(end+1,1)=gamma_ground_deg; %#ok<SAGROW>
                    route_end_cause(end+1,1)=cause; %#ok<SAGROW>
                    route_realization(end+1,1)=r; %#ok<SAGROW>
                    route_user_lat_A_deg(end+1,1)=rad2deg(user_lat(1)); %#ok<SAGROW>
                    route_user_lat_B_deg(end+1,1)=rad2deg(user_lat(2)); %#ok<SAGROW>
                    route_start_phase_s(end+1,1)=route_start_phase; %#ok<SAGROW>
                end

                route_active=false;
                route_path=[];
            end
        end

        %% Nouvelle route
        if ~route_active && both_assigned
            a=current_assignment(1);
            b=current_assignment(2);

            if a==b
                path_now=a;
            else
                path_now=shortestpath(G,a,b,'Method','unweighted');
            end

            if ~isempty(path_now)
                route_active=true;
                route_start_time=t;
                route_start_phase=mod(t,T_spatial);
                route_path=path_now;
                access_A_initial=a;
                access_B_initial=b;
                route_hops=max(0,numel(path_now)-1);
            end
        end

        if route_active
            samples_route_exists=samples_route_exists+1;
        end
    end

    fprintf('Realisation %d/%d terminee, episodes=%d\n', ...
        r,N_realizations,numel(route_durations));
end

%% Statistiques globales
if isempty(route_durations)
    error('Aucun episode complet observe.');
end

Troute_emp_mean = mean(route_durations);
Troute_emp_median = median(route_durations);
Troute_emp_std = std(route_durations);
Troute_emp_sem = Troute_emp_std/sqrt(numel(route_durations));

ci95_low = Troute_emp_mean-1.96*Troute_emp_sem;
ci95_high = Troute_emp_mean+1.96*Troute_emp_sem;

H_route_emp_mean = mean(route_hops_initial);

P_both_assigned = samples_both_assigned/total_samples;
P_route_available = samples_route_exists/total_samples;

%% Taux empiriques
total_route_exposure = sum(route_durations);

N_end_A = nnz(route_end_cause=="assignment_A");
N_end_B = nnz(route_end_cause=="assignment_B");
N_end_ISL = nnz(route_end_cause=="ISL_break");

beta_A_emp = N_end_A/total_route_exposure;
beta_B_emp = N_end_B/total_route_exposure;
beta_GSL_emp = beta_A_emp+beta_B_emp;
beta_ISL_emp = N_end_ISL/total_route_exposure;
beta_total_emp_from_causes = beta_GSL_emp+beta_ISL_emp;

if direct_link_exposure>0
    beta_link_emp_direct = direct_ISL_break_count/direct_link_exposure;
else
    beta_link_emp_direct = NaN;
end

beta_link_emp = beta_link_emp_direct;

%% Statistiques selon la phase de debut
N_time_bins = 24;
time_edges = linspace(0,T_spatial,N_time_bins+1).';
time_phase = 0.5*(time_edges(1:end-1)+time_edges(2:end));

Troute_emp_by_time = NaN(N_time_bins,1);
H_route_emp_by_time = NaN(N_time_bins,1);
EpisodeCount_by_time = zeros(N_time_bins,1);

for b=1:N_time_bins
    mask = route_start_phase_s>=time_edges(b) & ...
           route_start_phase_s<time_edges(b+1);

    EpisodeCount_by_time(b)=nnz(mask);

    if any(mask)
        Troute_emp_by_time(b)=mean(route_durations(mask));
        H_route_emp_by_time(b)=mean(route_hops_initial(mask));
    end
end

%% Statistiques H(t,gamma)
gamma_edges = (0:10:180).';
gamma_centers = 0.5*(gamma_edges(1:end-1)+gamma_edges(2:end));

rows = N_time_bins*numel(gamma_centers);
Time_s = NaN(rows,1);
GammaCenter_deg = NaN(rows,1);
EpisodeCount_tg = zeros(rows,1);
MeanHopCount_tg = NaN(rows,1);

ii=0;
for bt=1:N_time_bins
    for bg=1:numel(gamma_centers)
        ii=ii+1;
        Time_s(ii)=time_phase(bt);
        GammaCenter_deg(ii)=gamma_centers(bg);

        mask = route_start_phase_s>=time_edges(bt) & ...
               route_start_phase_s<time_edges(bt+1) & ...
               route_gamma_ground_deg>=gamma_edges(bg) & ...
               route_gamma_ground_deg<gamma_edges(bg+1);

        EpisodeCount_tg(ii)=nnz(mask);

        if any(mask)
            MeanHopCount_tg(ii)=mean(route_hops_initial(mask));
        end
    end
end

HRouteStatisticsTimeGamma = table( ...
    Time_s,GammaCenter_deg,EpisodeCount_tg,MeanHopCount_tg, ...
    'VariableNames',{'Time_s','GammaCenter_deg','EpisodeCount','MeanHopCount'});

%% Tables
EpisodeResults = table( ...
    route_realization,route_durations,route_hops_initial, ...
    route_gamma_ground_deg,route_user_lat_A_deg, ...
    route_user_lat_B_deg,route_start_phase_s,route_end_cause, ...
    'VariableNames',{'Realization','RouteLifetime_s','InitialHopCount', ...
    'GroundAngularSeparation_deg','UserLatitudeA_deg', ...
    'UserLatitudeB_deg','StartPhase_s','EndCause'});

Summary = table( ...
    numel(route_durations),Troute_emp_mean,Troute_emp_median, ...
    Troute_emp_std,ci95_low,ci95_high,H_route_emp_mean, ...
    P_both_assigned,P_route_available, ...
    beta_GSL_emp,beta_ISL_emp,beta_link_emp_direct, ...
    'VariableNames',{'NumberOfCompleteEpisodes', ...
    'MeanEmpiricalRouteLifetime_s','MedianEmpiricalRouteLifetime_s', ...
    'StdEmpiricalRouteLifetime_s','CI95Low_s','CI95High_s', ...
    'MeanInitialHopCount','ProbabilityBothUsersAssigned', ...
    'ProbabilityRouteAvailable','BetaGSL_Emp_per_s', ...
    'BetaISL_Emp_per_s','BetaLinkDirect_Emp_per_s'});

disp(Summary);

%% Figures
figure;
histogram(route_durations,'Normalization','probability');
grid on;
xlabel('Duree de vie de la route (s)');
ylabel('Probabilite');
title('Distribution empirique de T_{route} - Delta spatial');

figure;
plot(time_phase,Troute_emp_by_time,'o-','LineWidth',1.5);
grid on;
xlabel('Phase temporelle (s)');
ylabel('T_{route}^{emp} moyen (s)');
title('T_{route}^{emp} selon la phase spatiale');

%% Sauvegarde
save('T_route_emp_results.mat', ...
    'EpisodeResults','Summary','HRouteStatisticsTimeGamma', ...
    'route_durations','route_hops_initial', ...
    'route_gamma_ground_deg','route_start_phase_s','route_end_cause', ...
    'Troute_emp_mean','Troute_emp_median','Troute_emp_std', ...
    'ci95_low','ci95_high','H_route_emp_mean', ...
    'P_both_assigned','P_route_available', ...
    'beta_A_emp','beta_B_emp','beta_GSL_emp','beta_ISL_emp', ...
    'beta_link_emp','beta_link_emp_direct', ...
    'direct_link_exposure','direct_ISL_break_count', ...
    'Troute_emp_by_time','H_route_emp_by_time', ...
    'EpisodeCount_by_time','time_phase','time_edges', ...
    'lambda','inc_deg','dmax','elevation_min_deg', ...
    'N_realizations','dt','Tmax','T_orbit','T_spatial');

function u=sample_u_spatial(N)
    s=2*rand(N,1)-1;
    a=asin(s);
    branch=rand(N,1)<0.5;
    u=a;
    u(~branch)=pi-a(~branch);
    u=mod(u,2*pi);
end

function positions=walker_delta_positions(R,inc,Omega,u)
    x=R*(cos(Omega).*cos(u)-sin(Omega).*sin(u).*cos(inc));
    y=R*(sin(Omega).*cos(u)+cos(Omega).*sin(u).*cos(inc));
    z=R*(sin(u).*sin(inc));
    positions=[x y z];
end

function positions=ground_user_positions(R,lat,lon)
    positions=[R*cos(lat).*cos(lon), ...
               R*cos(lat).*sin(lon), ...
               R*sin(lat)];
end

function sat_id=closest_visible_satellite( ...
    sat_pos,user_pos,R_earth,elevation_min)

    rho=sat_pos-user_pos;
    dist=sqrt(sum(rho.^2,2));
    zenith=user_pos/R_earth;
    sin_el=(rho*zenith.')./dist;
    el=asin(max(-1,min(1,sin_el)));
    visible=el>=elevation_min;

    if any(visible)
        ids=find(visible);
        [~,j]=min(dist(ids));
        sat_id=ids(j);
    else
        sat_id=0;
    end
end

function alive=fixed_path_is_alive(A,path)
    if numel(path)<=1
        alive=true;
        return;
    end

    idx=sub2ind(size(A),path(1:end-1),path(2:end));
    alive=all(A(idx));
end
