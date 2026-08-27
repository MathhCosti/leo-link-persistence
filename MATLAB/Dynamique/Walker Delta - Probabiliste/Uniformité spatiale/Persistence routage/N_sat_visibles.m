clear; clc; close all;

%% NOMBRE MOYEN DE SATELLITES VISIBLES - DELTA UNIFORMITE SPATIALE
rng(1);

R_earth = 6371;
h = 550;
R = R_earth+h;
mu = 398600;
omega_sat = sqrt(mu/R^3);
omega_earth = 2*pi/86164;
T_orbit = 2*pi/omega_sat;

lambda = 4e-7;
N_mean = lambda*4*pi*R^2;

inc_deg = 58;
inc = deg2rad(inc_deg);

elevation_min_deg = 20;
elevation_min = deg2rad(elevation_min_deg);
psi_max = acos((R_earth/R)*cos(elevation_min))-elevation_min;
A_visible = 2*pi*R_earth^2*(1-cos(psi_max));

user_lat_deg = [-50 -35 -20 0 20 35 50].';
user_lat = deg2rad(user_lat_deg);
N_users = numel(user_lat);
user_lon0 = 2*pi*rand(N_users,1)-pi;

dt = 60;
time_values = (0:dt:T_orbit).';
Nt = numel(time_values);

N_realizations = 100;
sum_visible_t = zeros(Nt,N_users);

for r = 1:N_realizations
    N_sat = poissrnd(N_mean);
    Omega = 2*pi*rand(N_sat,1);
    u0 = sample_u_spatial(N_sat);

    for k = 1:Nt
        t = time_values(k);
        u_t = mod(u0+omega_sat*t,2*pi);

        x_sat = R*(cos(Omega).*cos(u_t)-sin(Omega).*sin(u_t).*cos(inc));
        y_sat = R*(sin(Omega).*cos(u_t)+cos(Omega).*sin(u_t).*cos(inc));
        z_sat = R*(sin(u_t).*sin(inc));
        sat_pos = [x_sat y_sat z_sat];

        lon_t = mod(user_lon0+omega_earth*t+pi,2*pi)-pi;
        user_pos = [R_earth*cos(user_lat).*cos(lon_t), ...
                    R_earth*cos(user_lat).*sin(lon_t), ...
                    R_earth*sin(user_lat)];

        for q = 1:N_users
            rho = sat_pos-user_pos(q,:);
            dist = sqrt(sum(rho.^2,2));
            zenith = user_pos(q,:)/R_earth;
            sin_el = (rho*zenith.') ./ dist;
            el = asin(max(-1,min(1,sin_el)));
            sum_visible_t(k,q) = sum_visible_t(k,q)+nnz(el >= elevation_min);
        end
    end
end

MeanVisible_emp = sum_visible_t/N_realizations;

MeanVisible_exact = zeros(Nt,N_users);
MeanVisible_local = zeros(Nt,N_users);
p_vis_exact = zeros(Nt,N_users);
rho_sat_local = zeros(Nt,N_users);

Nu = 20000;
u_grid = linspace(0,2*pi,Nu);
phi_s = asin(sin(inc).*sin(u_grid));

for k = 1:Nt
    t = time_values(k);
    fU_t = abs(cos(u_grid-omega_sat*t))/4;

    for q = 1:N_users
        phi_u = user_lat(q);

        denominator = cos(phi_u).*cos(phi_s);
        numerator = cos(psi_max)-sin(phi_u).*sin(phi_s);
        qq = numerator ./ denominator;

        fraction = zeros(size(u_grid));
        full = qq <= -1;
        none = qq >= 1;
        partial = ~(full | none);
        fraction(full)=1;
        fraction(none)=0;
        fraction(partial)=acos(qq(partial))/pi;

        singular = abs(denominator)<1e-14;
        if any(singular)
            angular_distance = abs(phi_s(singular)-phi_u);
            fraction(singular)=double(angular_distance<=psi_max);
        end

        p_vis_exact(k,q) = trapz(u_grid,fU_t.*fraction);
        MeanVisible_exact(k,q) = N_mean*p_vis_exact(k,q);

        fphi_local = latitude_pdf_spatial(phi_u,t,inc,omega_sat);
        rho_sat_local(k,q) = N_mean*fphi_local/(2*pi*R_earth^2*max(cos(phi_u),eps));
        MeanVisible_local(k,q) = rho_sat_local(k,q)*A_visible;
    end
end

err_exact = abs(MeanVisible_emp-MeanVisible_exact)./max(MeanVisible_emp,eps);

fprintf('\nErreur relative moyenne theorie exacte : %.2f %%\n', ...
    100*mean(err_exact,'all','omitnan'));

figure;
imagesc(user_lat_deg,time_values,MeanVisible_emp);
set(gca,'YDir','normal'); colorbar;
xlabel('Latitude utilisateur (deg)'); ylabel('Temps (s)');
title('N_{vis} moyen empirique - uniformite spatiale');

figure;
imagesc(user_lat_deg,time_values,MeanVisible_exact);
set(gca,'YDir','normal'); colorbar;
xlabel('Latitude utilisateur (deg)'); ylabel('Temps (s)');
title('N_{vis} moyen theorique - uniformite spatiale');

save('N_sat_visibles_results.mat');

function u = sample_u_spatial(N)
    s = 2*rand(N,1)-1;
    a = asin(s);
    branch = rand(N,1)<0.5;
    u = a;
    u(~branch)=pi-a(~branch);
    u = mod(u,2*pi);
end

function fphi = latitude_pdf_spatial(phi,t,inc,omega)
    if abs(phi)>=inc
        fphi=0; return;
    end
    x=max(-1,min(1,sin(phi)/sin(inc)));
    up=asin(x);
    um=pi-up;
    fUp=abs(cos(up-omega*t))/4;
    fUm=abs(cos(um-omega*t))/4;
    fphi=cos(phi)/sqrt(max(sin(inc)^2-sin(phi)^2,eps))*(fUp+fUm);
end
