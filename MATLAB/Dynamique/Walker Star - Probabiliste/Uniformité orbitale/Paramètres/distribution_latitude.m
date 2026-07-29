clear; clc; close all;

%% TEST DE LA DISTRIBUTION SPATIALE
R = 1;
inc_deg = 90;
inc = deg2rad(inc_deg);
M = 1e6;
nbins = 120;
rng('default');

Omega = 2*pi*rand(M,1);
u = 2*pi*rand(M,1);
positions = walker_delta_positions(R, inc, Omega, u);

x = positions(:,1);
y = positions(:,2);
z = positions(:,3);
longitude = atan2(y,x);
latitude = asin(z/R);
sin_latitude = z/R;

%% Latitude
lat_edges = linspace(-inc, inc, nbins+1);
lat_centers = 0.5*(lat_edges(1:end-1)+lat_edges(2:end));
[counts_lat,~] = histcounts(latitude, lat_edges, 'Normalization','pdf');

f_lat_theory = cos(lat_centers) ./ ...
    (pi*sqrt(max(eps, sin(inc)^2 - sin(lat_centers).^2)));
f_lat_uniform_band = cos(lat_centers)/(2*sin(inc));

figure; hold on; grid on;
bar(rad2deg(lat_centers), counts_lat, 1, 'FaceAlpha',0.35, ...
    'EdgeColor','none', 'DisplayName','Simulation');
plot(rad2deg(lat_centers), f_lat_theory, 'LineWidth',2, ...
    'DisplayName','Theorie induite par u uniforme');
plot(rad2deg(lat_centers), f_lat_uniform_band, '--', 'LineWidth',2, ...
    'DisplayName','Uniforme dans la bande');
xlabel('Latitude (deg)'); ylabel('Densite de probabilite');
title(sprintf('Distribution des latitudes, i = %.1f deg', inc_deg));
legend('Location','best'); hold off;

%% sin(latitude)
smax = sin(inc);
s_edges = linspace(-smax, smax, nbins+1);
s_centers = 0.5*(s_edges(1:end-1)+s_edges(2:end));
[counts_s,~] = histcounts(sin_latitude, s_edges, 'Normalization','pdf');

f_s_theory = 1 ./ (pi*sqrt(max(eps, smax^2 - s_centers.^2)));
f_s_uniform_band = ones(size(s_centers))/(2*smax);

figure; hold on; grid on;
bar(s_centers, counts_s, 1, 'FaceAlpha',0.35, ...
    'EdgeColor','none', 'DisplayName','Simulation');
plot(s_centers, f_s_theory, 'LineWidth',2, ...
    'DisplayName','Theorie induite');
plot(s_centers, f_s_uniform_band, '--', 'LineWidth',2, ...
    'DisplayName','Uniforme dans la bande');
xlabel('sin(latitude)'); ylabel('Densite de probabilite');
title('Distribution de sin(latitude)');
legend('Location','best'); hold off;

%% Longitude
lon_edges = linspace(-pi, pi, nbins+1);
lon_centers = 0.5*(lon_edges(1:end-1)+lon_edges(2:end));
[counts_lon,~] = histcounts(longitude, lon_edges, 'Normalization','pdf');

figure; hold on; grid on;
bar(rad2deg(lon_centers), counts_lon, 1, 'FaceAlpha',0.35, ...
    'EdgeColor','none', 'DisplayName','Simulation');
yline(1/(2*pi), '--', 'LineWidth',2, ...
    'DisplayName','Uniforme sur [-pi,pi]');
xlabel('Longitude (deg)'); ylabel('Densite de probabilite');
title('Distribution des longitudes');
legend('Location','best'); hold off;

%% Nuage 3D
Mplot = min(10000, M);
idx = randperm(M, Mplot);
figure;
scatter3(x(idx), y(idx), z(idx), 8, 'filled');
axis equal; grid on; view(3);
xlabel('x'); ylabel('y'); zlabel('z');
title(sprintf('%d positions, inclinaison fixe %.1f deg', Mplot, inc_deg));

%% Test KS contre uniforme dans la bande en surface
U_band = (sin_latitude + smax)/(2*smax);
[h_band,p_band] = kstest(U_band);

fprintf('\n=== Test de distribution spatiale ===\n');
fprintf('Nombre de points            : %d\n', M);
fprintf('Inclinaison                 : %.2f deg\n', inc_deg);
fprintf('Latitude min / max          : %.3f / %.3f deg\n', ...
    rad2deg(min(latitude)), rad2deg(max(latitude)));
fprintf('Moyenne latitude            : %.6f deg\n', rad2deg(mean(latitude)));
fprintf('Moyenne sin(latitude)       : %.6e\n', mean(sin_latitude));
fprintf('Test KS uniforme bande      : h = %d, p = %.6e\n', h_band, p_band);

if h_band == 1
    fprintf('Conclusion : distribution non uniforme dans la bande.\n');
else
    fprintf('Conclusion : uniformite dans la bande non rejetee.\n');
end

function positions = walker_delta_positions(R, inc, Omega, u)
    x = R * (cos(Omega).*cos(u) - sin(Omega).*sin(u).*cos(inc));
    y = R * (sin(Omega).*cos(u) + cos(Omega).*sin(u).*cos(inc));
    z = R * (sin(u).*sin(inc));
    positions = [x y z];
end
