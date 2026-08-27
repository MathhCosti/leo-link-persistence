clear; clc; close all;

%% ============================================================
% DISTRIBUTION DU NOMBRE DE SAUTS - DELTA UNIFORMITE SPATIALE
%
% La constellation est initialement uniforme en surface :
%   f_U(u,0)=|cos(u)|/4.
%
% La topologie du graphe varie ensuite periodiquement. On echantillonne
% plusieurs instants dans la demi-periode fondamentale T_orbit/2.
%
% Sorties :
%   - statistiques globales agregees ;
%   - statistiques conditionnelles selon gamma ;
%   - statistiques selon (t,gamma), utiles pour T_route spatial.
%% ============================================================

rng(4);

%% Parametres
R_earth = 6371;
h = 550;
R = R_earth+h;
mu = 398600;
omega_sat = sqrt(mu/R^3);
T_orbit = 2*pi/omega_sat;
T_spatial = T_orbit/2;

lambda = 4e-7;
inc_deg = 58;
inc = deg2rad(inc_deg);

dmax = 1500;
alpha_max = 2*asin(min(1,dmax/(2*R)));

%% Instants testes
N_time_samples = 8;
time_samples = linspace(0,T_spatial,N_time_samples+1).';
time_samples(end) = []; % evite de dupliquer 0 et T_spatial

%% Monte-Carlo
N_realizations = 20;
N_pairs_per_realization = 3000;
H_display_max = 30;

%% Stockage
all_H = [];
all_gamma_deg = [];
all_Hmin = [];
all_detour = [];
all_connected = false(0,1);
all_time_s = [];
realization_id = [];

N_sat_realizations = [];
N_edges_realizations = [];
largest_component_fraction = [];

%% Boucles temps / realisations
for it = 1:N_time_samples
    t = time_samples(it);

    for r = 1:N_realizations
        N_mean = lambda*4*pi*R^2;
        N_sat = poissrnd(N_mean);

        if N_sat < 2
            continue;
        end

        Omega = 2*pi*rand(N_sat,1);
        u0 = sample_u_spatial(N_sat);
        u_t = mod(u0+omega_sat*t,2*pi);

        positions = walker_delta_positions(R,inc,Omega,u_t);

        D = squareform(pdist(positions));
        A = (D <= dmax) & (D > 0);
        G = graph(A);

        bins = conncomp(G);
        component_sizes = accumarray(bins.',1);

        N_sat_realizations(end+1,1) = N_sat; %#ok<SAGROW>
        N_edges_realizations(end+1,1) = numedges(G); %#ok<SAGROW>
        largest_component_fraction(end+1,1) = ...
            max(component_sizes)/N_sat; %#ok<SAGROW>

        s = randi(N_sat,N_pairs_per_realization,1);
        d = randi(N_sat,N_pairs_per_realization,1);

        same = (s==d);
        while any(same)
            d(same) = randi(N_sat,nnz(same),1);
            same = (s==d);
        end

        H = inf(N_pairs_per_realization,1);
        gamma = zeros(N_pairs_per_realization,1);
        Hmin = zeros(N_pairs_per_realization,1);

        for p = 1:N_pairs_per_realization
            a = s(p); b = d(p);

            cg = dot(positions(a,:),positions(b,:))/R^2;
            cg = max(-1,min(1,cg));
            gamma(p) = acos(cg);

            Hmin(p) = max(1,ceil(gamma(p)/alpha_max));

            path = shortestpath(G,a,b,'Method','unweighted');
            if ~isempty(path)
                H(p) = numel(path)-1;
            end
        end

        connected = isfinite(H);

        all_H = [all_H;H]; %#ok<AGROW>
        all_gamma_deg = [all_gamma_deg;rad2deg(gamma)]; %#ok<AGROW>
        all_Hmin = [all_Hmin;Hmin]; %#ok<AGROW>
        all_detour = [all_detour;H-Hmin]; %#ok<AGROW>
        all_connected = [all_connected;connected]; %#ok<AGROW>
        all_time_s = [all_time_s;t*ones(N_pairs_per_realization,1)]; %#ok<AGROW>
        realization_id = [realization_id;r*ones(N_pairs_per_realization,1)]; %#ok<AGROW>
    end

    fprintf('Instant %d/%d termine : t=%.1f s\n',it,N_time_samples,t);
end

%% Statistiques globales
connected_mask = logical(all_connected);

H_connected = all_H(connected_mask);
Hmin_connected = all_Hmin(connected_mask);
detour_connected = all_detour(connected_mask);

P_connected = mean(connected_mask);
P_no_path = 1-P_connected;

mean_H = mean(H_connected);
median_H = median(H_connected);
mean_Hmin = mean(Hmin_connected);
mean_detour = mean(detour_connected);

%% Distribution empirique globale de H
max_H = max(H_connected);
h_values = (1:max_H).';
h_edges = [h_values.'-0.5,h_values(end)+0.5];

counts_H = histcounts(H_connected,h_edges).';
prob_H = counts_H/sum(counts_H);
cdf_H = cumsum(prob_H);

HopDistribution = table(h_values,counts_H,prob_H,cdf_H, ...
    'VariableNames',{'HopCount','EmpiricalCount','EmpiricalProbability','EmpiricalCDF'});

%% Statistiques selon gamma
gamma_edges = (0:10:180).';
gamma_centers = 0.5*(gamma_edges(1:end-1)+gamma_edges(2:end));
N_gamma_bins = numel(gamma_centers);

PairCount = zeros(N_gamma_bins,1);
ConnectionProbability = NaN(N_gamma_bins,1);
MeanHopCount = NaN(N_gamma_bins,1);
MedianHopCount = NaN(N_gamma_bins,1);
MeanTheoreticalHopCount = NaN(N_gamma_bins,1);

for b = 1:N_gamma_bins
    mask = all_gamma_deg >= gamma_edges(b) & ...
           all_gamma_deg < gamma_edges(b+1);

    PairCount(b) = nnz(mask);
    if PairCount(b)==0, continue; end

    ConnectionProbability(b) = mean(all_connected(mask));

    cmask = mask & all_connected;
    if any(cmask)
        MeanHopCount(b) = mean(all_H(cmask));
        MedianHopCount(b) = median(all_H(cmask));
        MeanTheoreticalHopCount(b) = mean(all_Hmin(cmask));
    end
end

GammaStatistics = table( ...
    gamma_edges(1:end-1),gamma_edges(2:end),gamma_centers, ...
    PairCount,ConnectionProbability,MeanHopCount,MedianHopCount, ...
    MeanTheoreticalHopCount, ...
    'VariableNames',{'GammaMin_deg','GammaMax_deg','GammaCenter_deg', ...
    'PairCount','ConnectionProbability','MeanHopCount', ...
    'MedianHopCount','MeanTheoreticalHopCount'});

%% Statistiques selon (t,gamma)
rows = N_time_samples*N_gamma_bins;
Time_s = NaN(rows,1);
GammaCenter_deg = NaN(rows,1);
PairCount_tg = zeros(rows,1);
ConnectionProbability_tg = NaN(rows,1);
MeanHopCount_tg = NaN(rows,1);
MeanTheoreticalHopCount_tg = NaN(rows,1);

ii = 0;
for it = 1:N_time_samples
    for b = 1:N_gamma_bins
        ii = ii+1;
        Time_s(ii) = time_samples(it);
        GammaCenter_deg(ii) = gamma_centers(b);

        mask = abs(all_time_s-time_samples(it))<1e-9 & ...
               all_gamma_deg >= gamma_edges(b) & ...
               all_gamma_deg < gamma_edges(b+1);

        PairCount_tg(ii) = nnz(mask);
        if PairCount_tg(ii)==0, continue; end

        ConnectionProbability_tg(ii) = mean(all_connected(mask));

        cmask = mask & all_connected;
        if any(cmask)
            MeanHopCount_tg(ii) = mean(all_H(cmask));
            MeanTheoreticalHopCount_tg(ii) = mean(all_Hmin(cmask));
        end
    end
end

GammaStatisticsTime = table( ...
    Time_s,GammaCenter_deg,PairCount_tg, ...
    ConnectionProbability_tg,MeanHopCount_tg,MeanTheoreticalHopCount_tg, ...
    'VariableNames',{'Time_s','GammaCenter_deg','PairCount', ...
    'ConnectionProbability','MeanHopCount','MeanTheoreticalHopCount'});

%% Statistiques globales
GlobalStatistics = table( ...
    P_connected,P_no_path,mean_H,median_H,mean_Hmin,mean_detour, ...
    'VariableNames',{'ConnectionProbability','NoPathProbability', ...
    'MeanShortestPathHopCount','MedianShortestPathHopCount', ...
    'MeanTheoreticalHopCount','MeanExtraHops'});

disp(GlobalStatistics);

%% Figures
figure;
bar(h_values,prob_H);
grid on;
xlabel('Nombre de sauts H');
ylabel('Probabilite');
title('Distribution empirique de H - Delta spatial');
xlim([0.5,min(H_display_max,h_values(end))+0.5]);

figure;
plot(gamma_centers,ConnectionProbability,'o-','LineWidth',1.5);
grid on; ylim([0 1.05]);
xlabel('\gamma (deg)');
ylabel('P_{conn}');
title('Probabilite de connexion moyenne selon \gamma');

figure;
hold on;
for it = 1:N_time_samples
    mask = abs(GammaStatisticsTime.Time_s-time_samples(it))<1e-9;
    plot(GammaStatisticsTime.GammaCenter_deg(mask), ...
         GammaStatisticsTime.ConnectionProbability(mask), ...
         'LineWidth',1.2, ...
         'DisplayName',sprintf('t=%.0f s',time_samples(it)));
end
grid on; ylim([0 1.05]);
xlabel('\gamma (deg)');
ylabel('P_{conn}(t,\gamma)');
title('Connectivite selon le temps - Delta spatial');
legend('Location','best');

%% Sauvegarde
save('H_jumps_results.mat', ...
    'GlobalStatistics','HopDistribution','GammaStatistics', ...
    'GammaStatisticsTime', ...
    'all_H','all_gamma_deg','all_Hmin','all_detour', ...
    'all_connected','all_time_s', ...
    'lambda','inc_deg','dmax','R','alpha_max', ...
    'P_connected','P_no_path','mean_H','median_H', ...
    'mean_Hmin','mean_detour','time_samples','T_orbit','T_spatial');

function u = sample_u_spatial(N)
    s = 2*rand(N,1)-1;
    a = asin(s);
    branch = rand(N,1)<0.5;
    u = a;
    u(~branch) = pi-a(~branch);
    u = mod(u,2*pi);
end

function positions = walker_delta_positions(R,inc,Omega,u)
    x = R*(cos(Omega).*cos(u)-sin(Omega).*sin(u).*cos(inc));
    y = R*(sin(Omega).*cos(u)+cos(Omega).*sin(u).*cos(inc));
    z = R*(sin(u).*sin(inc));
    positions = [x y z];
end
