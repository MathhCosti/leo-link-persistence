clear; clc; close all;

%% ============================================================
% DISTRIBUTION EMPIRIQUE DU NOMBRE DE SAUTS DU PLUS COURT CHEMIN
% MODELE DELTA STOCHASTIQUE EN UNIFORMITE ORBITALE
%
% Le script :
%   1) genere une constellation Delta stochastique ;
%   2) construit le graphe ISL avec une distance maximale dmax ;
%   3) tire N_pairs paires de satellites distincts ;
%   4) calcule le plus court chemin en nombre de sauts ;
%   5) estime la distribution P(H=h | H<Inf) ;
%   6) mesure aussi la probabilite d'absence de chemin ;
%   7) compare H a la borne geometrique
%          Hmin = ceil(gamma/alpha_max).
%% ============================================================

rng(4);

%% Parametres physiques et constellation
R_earth = 6371;              % km
h = 550;                     % km
R = R_earth + h;             % rayon orbital, km
mu = 398600;                 % km^3/s^2
omega_sat = sqrt(mu/R^3);    % rad/s

lambda = 4e-7;               % satellites/km^2 sur la sphere orbitale
inc_deg = 58;                % deg
inc = deg2rad(inc_deg);

dmax = 1500;                 % km
t = 0;                       % instant observe, s

%% Parametres Monte-Carlo
N_realizations = 30;         % nombre de graphes independants
N_pairs_per_realization = 5000;

% Nombre maximal de sauts affiche dans l'histogramme.
% Les valeurs plus grandes restent conservees dans les donnees.
H_display_max = 30;

N_pairs_theory = 1e6;

%% Stockage global
all_H = [];
all_gamma_deg = [];
all_Hmin = [];
all_detour = [];
all_connected = false(0,1);

realization_id = [];
source_id = [];
target_id = [];

N_sat_realizations = zeros(N_realizations,1);
N_edges_realizations = zeros(N_realizations,1);
largest_component_fraction = zeros(N_realizations,1);

%% Angle maximal couvert par un lien ISL
alpha_max = 2*asin(min(1,dmax/(2*R)));

%% Loi theorique Delta du nombre de sauts
Omega1_th=2*pi*rand(N_pairs_theory,1); Omega2_th=2*pi*rand(N_pairs_theory,1);
u1_th=2*pi*rand(N_pairs_theory,1); u2_th=2*pi*rand(N_pairs_theory,1);
pos1_th=walker_delta_positions(R,inc,Omega1_th,u1_th);
pos2_th=walker_delta_positions(R,inc,Omega2_th,u2_th);
cos_gamma_th=sum(pos1_th.*pos2_th,2)/R^2;
cos_gamma_th=max(-1,min(1,cos_gamma_th));
gamma_th=acos(cos_gamma_th);
H_theory=max(1,ceil(gamma_th/alpha_max));
mean_H_theory=mean(H_theory); median_H_theory=median(H_theory);
max_H_theory=max(H_theory);

%% Boucle sur les realisations
for r = 1:N_realizations

    %% Generation Delta stochastique
    N_mean = lambda*4*pi*R^2;
    N_sat = poissrnd(N_mean);

    if N_sat < 2
        warning('Realisation %d ignoree : moins de deux satellites.',r);
        continue;
    end

    Omega = 2*pi*rand(N_sat,1);
    u0 = 2*pi*rand(N_sat,1);
    u_t = mod(u0 + omega_sat*t,2*pi);

    positions = walker_delta_positions(R,inc,Omega,u_t);

    %% Graphe ISL
    D = squareform(pdist(positions));
    A = (D <= dmax) & (D > 0);
    G = graph(A);

    N_sat_realizations(r) = N_sat;
    N_edges_realizations(r) = numedges(G);

    bins = conncomp(G);
    component_sizes = accumarray(bins.',1);
    largest_component_fraction(r) = max(component_sizes)/N_sat;

    %% Tirage aleatoire de paires distinctes
    s = randi(N_sat,N_pairs_per_realization,1);
    d = randi(N_sat,N_pairs_per_realization,1);

    same = (s == d);
    while any(same)
        d(same) = randi(N_sat,nnz(same),1);
        same = (s == d);
    end

    %% Calcul des plus courts chemins
    H = inf(N_pairs_per_realization,1);
    gamma = zeros(N_pairs_per_realization,1);
    Hmin = zeros(N_pairs_per_realization,1);

    for p = 1:N_pairs_per_realization
        a = s(p);
        b = d(p);

        cos_gamma = dot(positions(a,:),positions(b,:))/R^2;
        cos_gamma = max(-1,min(1,cos_gamma));
        gamma(p) = acos(cos_gamma);

        Hmin(p) = ceil(gamma(p)/alpha_max);

        path = shortestpath(G,a,b,'Method','unweighted');

        if ~isempty(path)
            H(p) = numel(path)-1;

        end
    end

    connected = isfinite(H);

    all_H = [all_H; H]; %#ok<AGROW>
    all_gamma_deg = [all_gamma_deg; rad2deg(gamma)]; %#ok<AGROW>
    all_Hmin = [all_Hmin; Hmin]; %#ok<AGROW>
    all_detour = [all_detour; H-Hmin]; %#ok<AGROW>
    all_connected = [all_connected; connected]; %#ok<AGROW>

    realization_id = [realization_id; r*ones(N_pairs_per_realization,1)]; %#ok<AGROW>
    source_id = [source_id; s]; %#ok<AGROW>
    target_id = [target_id; d]; %#ok<AGROW>

    fprintf(['Realisation %2d/%2d : N=%d, E=%d, ' ...
             'paires connectees=%.3f\n'], ...
             r,N_realizations,N_sat,numedges(G),mean(connected));
end

%% Nettoyage si certaines realisations ont ete ignorees
valid_rows = ~isnan(realization_id);
realization_id = realization_id(valid_rows);
source_id = source_id(valid_rows);
target_id = target_id(valid_rows);

%% Statistiques globales
connected_mask = logical(all_connected);
H_connected = all_H(connected_mask);
gamma_connected = all_gamma_deg(connected_mask);
Hmin_connected = all_Hmin(connected_mask);
detour_connected = all_detour(connected_mask);

P_no_path = mean(~connected_mask);
P_connected = mean(connected_mask);

mean_H = mean(H_connected);
median_H = median(H_connected);
min_H = min(H_connected);
max_H = max(H_connected);

mean_Hmin = mean(Hmin_connected);
mean_detour = mean(detour_connected);

fprintf('\n============================================================\n');
fprintf('Nombre total de paires testees : %d\n',numel(all_H));
fprintf('Probabilite empirique de connexion : %.4f\n',P_connected);
fprintf('Probabilite empirique d''absence de chemin : %.4f\n',P_no_path);
fprintf('Nombre moyen de sauts conditionnel : %.3f\n',mean_H);
fprintf('Mediane conditionnelle : %.3f\n',median_H);
fprintf('Minimum / maximum observes : %d / %d\n',min_H,max_H);
fprintf('Borne geometrique moyenne Hmin : %.3f\n',mean_Hmin);
fprintf('Surplus moyen H-Hmin : %.3f\n',mean_detour);
fprintf('============================================================\n');

% Tableau recapitulatif des statistiques du plus court chemin.
GlobalStatistics = table( ...
    P_connected,P_no_path,mean_H,median_H,min_H,max_H, ...
    mean_Hmin,mean_detour, ...
    'VariableNames',{ ...
    'ConnectionProbability','NoPathProbability', ...
    'MeanShortestPathHopCount','MedianShortestPathHopCount', ...
    'MinimumHopCount','MaximumHopCount', ...
    'MeanTheoreticalHopCount','MeanExtraHops'});

disp(GlobalStatistics);

%% Distributions empirique et theorique de H
h_values=(1:max(max_H,max_H_theory)).';
h_edges=[h_values.'-0.5,h_values(end)+0.5];
counts_H=histcounts(H_connected,h_edges).'; prob_H=counts_H/sum(counts_H); cdf_H=cumsum(prob_H);
counts_H_theory=histcounts(H_theory,h_edges).'; prob_H_theory=counts_H_theory/sum(counts_H_theory); cdf_H_theory=cumsum(prob_H_theory);
mean_error_theory=mean_H_theory-mean_H;
L1_probability_error=sum(abs(prob_H_theory-prob_H));
KS_CDF_error=max(abs(cdf_H_theory-cdf_H));
HopDistribution=table(h_values,counts_H,prob_H,cdf_H,counts_H_theory,prob_H_theory,cdf_H_theory, ...
 'VariableNames',{'HopCount','EmpiricalCount','EmpiricalProbability','EmpiricalCDF','TheoreticalCount','TheoreticalProbability','TheoreticalCDF'});
disp(HopDistribution);
fprintf('Moyenne theorique Delta : %.3f\n',mean_H_theory);
fprintf('Ecart moyenne theorie - empirique : %.3f\n',mean_error_theory);
fprintf('Erreur L1 sur les probabilites : %.3f\n',L1_probability_error);
fprintf('Ecart maximal entre les CDF : %.3f\n',KS_CDF_error);

%% Table complete des paires
PairResults = table( ...
    realization_id,source_id,target_id, ...
    all_gamma_deg,all_Hmin,all_H,all_detour,all_connected, ...
    'VariableNames',{ ...
    'Realization','SourceID','TargetID', ...
    'AngularSeparation_deg','TheoreticalHopCount', ...
    'ShortestPathHopCount','ExtraHops','Connected'});

%% Distribution conditionnelle selon la distance angulaire
gamma_edges = 0:15:180;
gamma_centers = (gamma_edges(1:end-1)+gamma_edges(2:end))/2;
N_gamma_bins = numel(gamma_centers);

MeanH_by_gamma = NaN(N_gamma_bins,1);
MedianH_by_gamma = NaN(N_gamma_bins,1);
Pconnected_by_gamma = NaN(N_gamma_bins,1);
MeanHmin_by_gamma = NaN(N_gamma_bins,1);
Npairs_by_gamma = zeros(N_gamma_bins,1);

for b = 1:N_gamma_bins
    in_bin = all_gamma_deg >= gamma_edges(b) & ...
             all_gamma_deg < gamma_edges(b+1);

    Npairs_by_gamma(b) = nnz(in_bin);

    if Npairs_by_gamma(b) == 0
        continue;
    end

    Pconnected_by_gamma(b) = mean(all_connected(in_bin));

    in_bin_connected = in_bin & all_connected;
    if any(in_bin_connected)
        MeanH_by_gamma(b) = mean(all_H(in_bin_connected));
        MedianH_by_gamma(b) = median(all_H(in_bin_connected));
        MeanHmin_by_gamma(b) = mean(all_Hmin(in_bin_connected));
    end
end

GammaStatistics = table( ...
    gamma_edges(1:end-1).',gamma_edges(2:end).',gamma_centers.', ...
    Npairs_by_gamma,Pconnected_by_gamma,MeanH_by_gamma, ...
    MedianH_by_gamma,MeanHmin_by_gamma, ...
    'VariableNames',{ ...
    'GammaMin_deg','GammaMax_deg','GammaCenter_deg','PairCount', ...
    'ConnectionProbability','MeanHopCount','MedianHopCount', ...
    'MeanTheoreticalHopCount'});

disp(GammaStatistics);

%% Graphique 1 : probabilites empirique et theorique Delta
figure;
bar(h_values,prob_H,0.85,'DisplayName','Empirique'); hold on;
plot(h_values,prob_H_theory,'o-','LineWidth',1.8,'DisplayName','Theorie : H=\lceil\gamma/\alpha_{max}\rceil');
xline(mean_H,'--','LineWidth',1.5,'DisplayName',sprintf('Moyenne empirique = %.2f',mean_H));
xline(mean_H_theory,'-.','LineWidth',1.5,'DisplayName',sprintf('Moyenne theorique = %.2f',mean_H_theory));
grid on; xlabel('Nombre de sauts H'); ylabel('Probabilite');
title(sprintf('Distribution de probabilités | i=%.0f deg, d_{max}=%.0f km',inc_deg,dmax));
xlim([0.5,min(h_values(end)+0.5,H_display_max)]); legend('Location','best'); hold off;

%% Graphique 2 : CDF empirique et theorique Delta
figure;
stairs(h_values,cdf_H,'LineWidth',1.8,'DisplayName','Empirique'); hold on;
stairs(h_values,cdf_H_theory,'--','LineWidth',1.8,'DisplayName','Theorie : H=\lceil\gamma/\alpha_{max}\rceil');
grid on; xlabel('Nombre de sauts h'); ylabel('P(H\leq h)');
title(sprintf('CDF | ecart max = %.3f',KS_CDF_error));
legend('Location','best'); hold off;

%% Graphique 3 : nombre de sauts selon gamma
figure;
plot(gamma_centers,MeanH_by_gamma,'o-','LineWidth',1.8, ...
    'DisplayName','Empirique');
hold on;
plot(gamma_centers,MeanHmin_by_gamma,'--','LineWidth',1.5, ...
    'DisplayName','Theorie : H=\lceil\gamma/\alpha_{max}\rceil');

grid on;
xlabel('Séparation angulaire \gamma (deg)');
ylabel('Nombre moyen de sauts');
title('Nombre de sauts du plus court chemin');
legend('Location','best');

%% Graphique 4 : probabilite de connexion selon gamma
figure;
plot(gamma_centers,Pconnected_by_gamma,'o-','LineWidth',1.5);
grid on;
ylim([0 1.05]);
xlabel('Separation angulaire \gamma (deg)');
ylabel('P(H<\infty | \gamma)');
title('Probabilite empirique d''existence d''un chemin');

%% Graphique 5 : surplus de sauts par rapport a la borne
detour_values = (min(detour_connected):max(detour_connected)).';
detour_counts = histcounts(detour_connected, ...
    [detour_values.'-0.5,detour_values(end)+0.5]).';
detour_prob = detour_counts/sum(detour_counts);

figure;
bar(detour_values,detour_prob);
grid on;
xlabel('Surplus K = H-H_{min}');
ylabel('P(K=k | H<\infty)');
title('Distribution des detours par rapport a la borne geometrique');

%% Sauvegarde
RealizationStatistics = table( ...
    (1:N_realizations).',N_sat_realizations,N_edges_realizations, ...
    largest_component_fraction, ...
    'VariableNames',{ ...
    'Realization','NumberOfSatellites','NumberOfEdges', ...
    'LargestComponentFraction'});

% writetable(GlobalStatistics,'shortest_path_global_statistics.csv');
% writetable(HopDistribution,'distribution_shortest_path_hops.csv');
% writetable(PairResults,'shortest_path_pairs.csv');
% writetable(GammaStatistics,'shortest_path_gamma_statistics.csv');
% writetable(RealizationStatistics,'shortest_path_realizations.csv');

save('H_jumps_results.mat', ...
    'GlobalStatistics','HopDistribution','PairResults', ...
    'GammaStatistics','RealizationStatistics', ...
    'lambda','inc_deg','dmax','R','alpha_max', ...
    'P_connected','P_no_path','mean_H','mean_H_theory','median_H', ...
    'mean_Hmin','mean_detour');

%% Fonctions
function positions = walker_delta_positions(R,inc,Omega,u)
    x = R*(cos(Omega).*cos(u) ...
        - sin(Omega).*sin(u).*cos(inc));

    y = R*(sin(Omega).*cos(u) ...
        + cos(Omega).*sin(u).*cos(inc));

    z = R*(sin(u).*sin(inc));

    positions = [x y z];
end
