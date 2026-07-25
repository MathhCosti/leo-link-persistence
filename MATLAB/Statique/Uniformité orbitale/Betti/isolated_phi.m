clear; clc; close all; rng(2);

R = 6371 + 550;
inc_deg = 53;
inc_rad = deg2rad(inc_deg);
N = 250;
dmax = 1500;
n_realizations = 600;
nbins = 20;

lat_edges = linspace(-inc_rad,inc_rad,nbins+1);
lat_centers = 0.5*(lat_edges(1:end-1)+lat_edges(2:end));

nodes_per_bin = zeros(nbins,1);
degree_sum_bin = zeros(nbins,1);
isolated_sum_bin = zeros(nbins,1);
isolated_emp = zeros(n_realizations,1);

for r = 1:n_realizations
    [positions,~,~,lat] = walker_delta_static_sample(N,R,inc_rad);
    D = squareform(pdist(positions));
    A = (D <= dmax) & (D > 0);
    deg = sum(A,2);
    isolated_emp(r) = sum(deg==0);

    bin_id = discretize(lat,lat_edges);
    for b = 1:nbins
        mask = (bin_id==b);
        nodes_per_bin(b) = nodes_per_bin(b) + sum(mask);
        degree_sum_bin(b) = degree_sum_bin(b) + sum(deg(mask));
        isolated_sum_bin(b) = isolated_sum_bin(b) + sum(deg(mask)==0);
    end
end

mean_degree_lat = degree_sum_bin ./ max(nodes_per_bin,1);
p_link_lat = mean_degree_lat/(N-1);
p_isolated_lat_emp = isolated_sum_bin ./ max(nodes_per_bin,1);
p_isolated_lat_th = (1-p_link_lat).^(N-1);

weights_lat = nodes_per_bin/sum(nodes_per_bin);
N1_th_local = N*sum(weights_lat.*p_isolated_lat_th);
p_link_global = sum(weights_lat.*p_link_lat);
N1_th_global = N*(1-p_link_global)^(N-1);

fprintf('\n=== Correction locale en latitude ===\n');
fprintf('Isoles simulation      : %.3f +/- %.3f\n',mean(isolated_emp),std(isolated_emp));
fprintf('Prediction locale      : %.3f\n',N1_th_local);
fprintf('Prediction globale     : %.3f\n',N1_th_global);

figure;
yyaxis left
plot(rad2deg(lat_centers),p_link_lat,'o-','LineWidth',1.6);
ylabel('p_{link}(\varphi)');
yyaxis right
plot(rad2deg(lat_centers),mean_degree_lat,'s-','LineWidth',1.6);
ylabel('E[k|\varphi]');
grid on; xlabel('Latitude (deg)');
title('Probabilite de lien et degre locaux');

figure;
plot(rad2deg(lat_centers),p_isolated_lat_emp,'o-','LineWidth',1.6); hold on;
plot(rad2deg(lat_centers),p_isolated_lat_th,'--','LineWidth',1.8);
grid on; xlabel('Latitude (deg)');
ylabel('Probabilite d''etre isole');
title('Satellites isoles');
legend('Simulation','Prediction locale','Location','best');

figure;
histogram(isolated_emp,30,'Normalization','pdf'); hold on;
xline(N1_th_local,'--','Prediction locale','LineWidth',1.8);
xline(N1_th_global,':','Prediction globale','LineWidth',1.8);
grid on; xlabel('Nombre de satellites isoles'); ylabel('Densite');
title('Nombre de satellites isoles');

save('verify_local_plink_isolated_results.mat');
