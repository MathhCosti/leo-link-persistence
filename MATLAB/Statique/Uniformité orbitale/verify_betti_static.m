clear; clc; close all; rng(3);

R = 6371 + 550;
inc_deg = 53;
inc_rad = deg2rad(inc_deg);
N = 250;
dmax = 1500;
n_realizations = 500;
nbins = 20;
C_macro = 1;

lat_edges = linspace(-inc_rad,inc_rad,nbins+1);
nodes_per_bin = zeros(nbins,1);
degree_sum_bin = zeros(nbins,1);

beta0_emp = zeros(n_realizations,1);
beta1_emp = zeros(n_realizations,1);
E_emp = zeros(n_realizations,1);
N1_emp = zeros(n_realizations,1);
N2_emp = zeros(n_realizations,1);
N3_emp = zeros(n_realizations,1);

for r = 1:n_realizations
    [positions,~,~,lat] = walker_delta_static_sample(N,R,inc_rad);
    D = squareform(pdist(positions));
    A = sparse((D <= dmax) & (D > 0));
    G = graph(A);

    comp = conncomp(G);
    comp_sizes = accumarray(comp',1);

    beta0_emp(r) = numel(comp_sizes);
    E_emp(r) = nnz(triu(A,1));
    beta1_emp(r) = E_emp(r)-N+beta0_emp(r);

    N1_emp(r) = sum(comp_sizes==1);
    N2_emp(r) = sum(comp_sizes==2);
    N3_emp(r) = sum(comp_sizes==3);

    deg = full(sum(A,2));
    bin_id = discretize(lat,lat_edges);
    for b = 1:nbins
        mask = (bin_id==b);
        nodes_per_bin(b) = nodes_per_bin(b)+sum(mask);
        degree_sum_bin(b) = degree_sum_bin(b)+sum(deg(mask));
    end
end

weights_lat = nodes_per_bin/sum(nodes_per_bin);
mean_degree_lat = degree_sum_bin./max(nodes_per_bin,1);
p_link_lat = mean_degree_lat/(N-1);
p_link_global = sum(weights_lat.*p_link_lat);

E_th = nchoosek(N,2)*p_link_global;
beta0_th_forest = max(N-E_th,1);

N1_th_global = N*(1-p_link_global)^(N-1);
N1_th_local = N*sum(weights_lat.*(1-p_link_lat).^(N-1));

beta0_th_global = C_macro + N1_th_global;
beta0_th_local = C_macro + N1_th_local;

beta1_th_forest = E_th-N+beta0_th_forest;
beta1_th_global = E_th-N+beta0_th_global;
beta1_th_local = E_th-N+beta0_th_local;

fprintf('\n=== Nombres de Betti ===\n');
fprintf('beta0 simulation : %.3f +/- %.3f\n',mean(beta0_emp),std(beta0_emp));
fprintf('beta0 foret      : %.3f\n',beta0_th_forest);
fprintf('beta0 global     : %.3f\n',beta0_th_global);
fprintf('beta0 local      : %.3f\n',beta0_th_local);
fprintf('beta1 simulation : %.3f +/- %.3f\n',mean(beta1_emp),std(beta1_emp));
fprintf('beta1 foret      : %.3f\n',beta1_th_forest);
fprintf('beta1 global     : %.3f\n',beta1_th_global);
fprintf('beta1 local      : %.3f\n',beta1_th_local);
fprintf('N1 moyen         : %.3f\n',mean(N1_emp));
fprintf('N2 moyen         : %.3f\n',mean(N2_emp));
fprintf('N3 moyen         : %.3f\n',mean(N3_emp));

figure;
histogram(beta0_emp,30,'Normalization','pdf'); hold on;
xline(beta0_th_forest,'--','Foret','LineWidth',1.6);
xline(beta0_th_global,':','Isoles global','LineWidth',1.8);
xline(beta0_th_local,'-.','Isoles local','LineWidth',1.8);
grid on; xlabel('\beta_0'); ylabel('Densite');
title('\beta_0 statique');

figure;
histogram(beta1_emp,30,'Normalization','pdf'); hold on;
xline(beta1_th_forest,'--','Foret','LineWidth',1.6);
xline(beta1_th_global,':','Isoles global','LineWidth',1.8);
xline(beta1_th_local,'-.','Isoles local','LineWidth',1.8);
grid on; xlabel('\beta_1'); ylabel('Densite');
title('\beta_1 statique');

figure;
bar([mean(N1_emp),mean(N2_emp),mean(N3_emp)]);
set(gca,'XTickLabel',{'Isoles','Dimeres','Trimeres'});
grid on; ylabel('Nombre moyen de composantes');
title('Petites composantes observees');

save('verify_betti_static_results.mat');
