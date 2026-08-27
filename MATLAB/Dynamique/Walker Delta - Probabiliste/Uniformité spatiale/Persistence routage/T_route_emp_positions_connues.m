clear; clc; close all;

%% ============================================================
% T_ROUTE EMPIRIQUE - POSITIONS UTILISATEURS CONNUES
% WALKER DELTA - UNIFORMITE SPATIALE
%
% Sortie principale :
%   Troute_emp_matrix(t,phi_A,phi_B)
%
% Chaque episode est classe selon sa phase de debut modulo T_orbit/2.
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

%% Latitudes
user_lat_deg = (-55:5:55).';
user_lat = deg2rad(user_lat_deg);
N_lat = numel(user_lat_deg);

%% Longitudes imposees
user_lon_A_deg = 0;
user_lon_B_deg = 90;

user_lon0_A_fixed = deg2rad(user_lon_A_deg);
user_lon0_B_fixed = deg2rad(user_lon_B_deg);

Delta_lon_deg = abs(user_lon_A_deg-user_lon_B_deg);

%% Simulation
N_realizations = 40;
dt = 5;
Tmax = 15000;
time_values = (0:dt:Tmax-dt).';
Nt = numel(time_values);

N_time_bins = 24;
time_edges = linspace(0,T_spatial,N_time_bins+1).';
time_phase = 0.5*(time_edges(1:end-1)+time_edges(2:end));

min_episode_samples = 1;

%% Stockage
route_durations = cell(N_time_bins,N_lat,N_lat);
route_hops = cell(N_time_bins,N_lat,N_lat);
route_end_cause = cell(N_time_bins,N_lat,N_lat);

samples_both_assigned = zeros(N_time_bins,N_lat,N_lat);
samples_route_exists = zeros(N_time_bins,N_lat,N_lat);
sample_count = zeros(N_time_bins,N_lat,N_lat);

%% Monte-Carlo
for r = 1:N_realizations

    N_mean = lambda*4*pi*R^2;
    N_sat = poissrnd(N_mean);

    if N_sat < 2, continue; end

    Omega = 2*pi*rand(N_sat,1);
    u0 = sample_u_spatial(N_sat);

    user_lon0_A = user_lon0_A_fixed*ones(N_lat,1);
    user_lon0_B = user_lon0_B_fixed*ones(N_lat,1);

    route_active = false(N_lat,N_lat);
    route_start_time = NaN(N_lat,N_lat);
    route_start_bin = NaN(N_lat,N_lat);

    access_A_initial = zeros(N_lat,N_lat);
    access_B_initial = zeros(N_lat,N_lat);

    route_hops_current = NaN(N_lat,N_lat);
    route_path = cell(N_lat,N_lat);

    for k = 1:Nt
        t = time_values(k);

        %% Satellites
        u_t = mod(u0+omega_sat*t,2*pi);

        x_sat = R*(cos(Omega).*cos(u_t)-sin(Omega).*sin(u_t).*cos(inc));
        y_sat = R*(sin(Omega).*cos(u_t)+cos(Omega).*sin(u_t).*cos(inc));
        z_sat = R*(sin(u_t).*sin(inc));
        sat_pos = [x_sat y_sat z_sat];

        %% Graphe ISL
        D = squareform(pdist(sat_pos));
        Adj = (D<=dmax) & (D>0);
        G = graph(Adj);

        %% Utilisateurs
        lon_A = mod(user_lon0_A+omega_earth*t+pi,2*pi)-pi;
        lon_B = mod(user_lon0_B+omega_earth*t+pi,2*pi)-pi;

        user_pos_A = [ ...
            R_earth*cos(user_lat).*cos(lon_A), ...
            R_earth*cos(user_lat).*sin(lon_A), ...
            R_earth*sin(user_lat)];

        user_pos_B = [ ...
            R_earth*cos(user_lat).*cos(lon_B), ...
            R_earth*cos(user_lat).*sin(lon_B), ...
            R_earth*sin(user_lat)];

        assignment_A = zeros(N_lat,1);
        assignment_B = zeros(N_lat,1);

        for q = 1:N_lat
            assignment_A(q) = closest_visible_satellite( ...
                sat_pos,user_pos_A(q,:),R_earth,elevation_min);

            assignment_B(q) = closest_visible_satellite( ...
                sat_pos,user_pos_B(q,:),R_earth,elevation_min);
        end

        phase_t = mod(t,T_spatial);
        b_now = discretize(phase_t,time_edges);
        if isnan(b_now), b_now=N_time_bins; end

        for qA = 1:N_lat
            for qB = 1:N_lat

                a_now = assignment_A(qA);
                b_now_sat = assignment_B(qB);

                both_assigned = (a_now>0)&&(b_now_sat>0);

                sample_count(b_now,qA,qB) = ...
                    sample_count(b_now,qA,qB)+1;

                samples_both_assigned(b_now,qA,qB) = ...
                    samples_both_assigned(b_now,qA,qB)+both_assigned;

                %% Persistance
                if route_active(qA,qB)
                    cause = "";

                    if a_now ~= access_A_initial(qA,qB)
                        cause = "assignment_A";

                    elseif b_now_sat ~= access_B_initial(qA,qB)
                        cause = "assignment_B";

                    else
                        path_current = route_path{qA,qB};

                        if numel(path_current)<=1
                            alive = true;
                        else
                            idx = sub2ind(size(Adj), ...
                                path_current(1:end-1), ...
                                path_current(2:end));
                            alive = all(Adj(idx));
                        end

                        if ~alive
                            cause = "ISL_break";
                        end
                    end

                    if strlength(cause)>0
                        duration = t-route_start_time(qA,qB);

                        if duration >= min_episode_samples*dt
                            b0 = route_start_bin(qA,qB);

                            route_durations{b0,qA,qB}(end+1,1) = duration;
                            route_hops{b0,qA,qB}(end+1,1) = ...
                                route_hops_current(qA,qB);
                            route_end_cause{b0,qA,qB}(end+1,1) = cause;
                        end

                        route_active(qA,qB)=false;
                        route_path{qA,qB}=[];
                    end
                end

                %% Nouvelle route
                if ~route_active(qA,qB) && both_assigned

                    if a_now==b_now_sat
                        path_now = a_now;
                    else
                        path_now = shortestpath(G,a_now,b_now_sat, ...
                            'Method','unweighted');
                    end

                    if ~isempty(path_now)
                        route_active(qA,qB)=true;
                        route_start_time(qA,qB)=t;
                        route_start_bin(qA,qB)=b_now;
                        route_path{qA,qB}=path_now;
                        access_A_initial(qA,qB)=a_now;
                        access_B_initial(qA,qB)=b_now_sat;
                        route_hops_current(qA,qB)=max(0,numel(path_now)-1);
                    end
                end

                if route_active(qA,qB)
                    samples_route_exists(b_now,qA,qB) = ...
                        samples_route_exists(b_now,qA,qB)+1;
                end
            end
        end
    end

    fprintf('Realisation %d/%d terminee, N=%d\n',r,N_realizations,N_sat);
end

%% Matrices
Troute_emp_matrix = NaN(N_time_bins,N_lat,N_lat);
Troute_median_matrix = NaN(N_time_bins,N_lat,N_lat);
Troute_std_matrix = NaN(N_time_bins,N_lat,N_lat);
EpisodeCount_matrix = zeros(N_time_bins,N_lat,N_lat);
MeanHopCount_matrix = NaN(N_time_bins,N_lat,N_lat);

for bt = 1:N_time_bins
    for qA = 1:N_lat
        for qB = 1:N_lat
            d = route_durations{bt,qA,qB};
            EpisodeCount_matrix(bt,qA,qB)=numel(d);

            if ~isempty(d)
                Troute_emp_matrix(bt,qA,qB)=mean(d);
                Troute_median_matrix(bt,qA,qB)=median(d);
                Troute_std_matrix(bt,qA,qB)=std(d);
            end

            hh = route_hops{bt,qA,qB};
            if ~isempty(hh)
                MeanHopCount_matrix(bt,qA,qB)=mean(hh);
            end
        end
    end
end

P_both_assigned_matrix = ...
    samples_both_assigned./max(sample_count,1);

P_route_available_matrix = ...
    samples_route_exists./max(sample_count,1);

Troute_time_mean_matrix = ...
    squeeze(mean(Troute_emp_matrix,1,'omitnan'));

%% Figures
snapshot_idx = unique(round(linspace(1,ceil(N_time_bins/2),4)));

for s=1:numel(snapshot_idx)
    bt=snapshot_idx(s);

    figure;
    imagesc(user_lat_deg,user_lat_deg,squeeze(Troute_emp_matrix(bt,:,:)));
    set(gca,'YDir','normal');
    colorbar;
    xlabel('\phi_B (deg)');
    ylabel('\phi_A (deg)');
    title(sprintf( ...
        'T_{route}^{emp}(t,\\phi_A,\\phi_B), t=%.0f s', ...
        time_phase(bt)));
end

figure;
imagesc(user_lat_deg,user_lat_deg,Troute_time_mean_matrix);
set(gca,'YDir','normal');
colorbar;
xlabel('\phi_B (deg)');
ylabel('\phi_A (deg)');
title('Moyenne temporelle empirique de T_{route}');

%% Sauvegarde
save('T_route_emp_positions_connues_results.mat', ...
    'Troute_emp_matrix','Troute_time_mean_matrix', ...
    'Troute_median_matrix','Troute_std_matrix', ...
    'EpisodeCount_matrix','MeanHopCount_matrix', ...
    'P_both_assigned_matrix','P_route_available_matrix', ...
    'route_durations','route_hops','route_end_cause', ...
    'time_phase','time_edges','user_lat_deg', ...
    'user_lon_A_deg','user_lon_B_deg','Delta_lon_deg', ...
    'lambda','inc_deg','dmax','elevation_min_deg', ...
    'N_realizations','dt','Tmax','T_orbit','T_spatial');

function u = sample_u_spatial(N)
    s=2*rand(N,1)-1;
    a=asin(s);
    branch=rand(N,1)<0.5;
    u=a;
    u(~branch)=pi-a(~branch);
    u=mod(u,2*pi);
end

function sat_id = closest_visible_satellite( ...
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
