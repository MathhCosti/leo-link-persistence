clear; clc; close all;

%% PROBABILITE D'OUTAGE - DELTA UNIFORMITE SPATIALE
rng(1);

R_earth=6371; h=550; R=R_earth+h; mu=398600;
omega_sat=sqrt(mu/R^3);
omega_earth=2*pi/86164;
T_orbit=2*pi/omega_sat;

lambda=4e-7;
N_mean=lambda*4*pi*R^2;

inc_deg=58;
inc=deg2rad(inc_deg);

elevation_min_deg=20;
elevation_min=deg2rad(elevation_min_deg);
psi_max=acos((R_earth/R)*cos(elevation_min))-elevation_min;
A_visible=2*pi*R_earth^2*(1-cos(psi_max));

user_lat_deg=[-50 -35 -20 0 20 35 50].';
user_lat=deg2rad(user_lat_deg);
N_users=numel(user_lat);
user_lon0=2*pi*rand(N_users,1)-pi;

dt=60;
time_values=(0:dt:T_orbit).';
Nt=numel(time_values);

N_realizations=1000;
sum_outage_t=zeros(Nt,N_users);

for r=1:N_realizations
    N_sat=poissrnd(N_mean);
    Omega=2*pi*rand(N_sat,1);
    u0=sample_u_spatial(N_sat);

    for k=1:Nt
        t=time_values(k);
        u_t=mod(u0+omega_sat*t,2*pi);

        x_sat=R*(cos(Omega).*cos(u_t)-sin(Omega).*sin(u_t).*cos(inc));
        y_sat=R*(sin(Omega).*cos(u_t)+cos(Omega).*sin(u_t).*cos(inc));
        z_sat=R*(sin(u_t).*sin(inc));
        sat_pos=[x_sat y_sat z_sat];

        lon_t=mod(user_lon0+omega_earth*t+pi,2*pi)-pi;
        user_pos=[R_earth*cos(user_lat).*cos(lon_t), ...
                  R_earth*cos(user_lat).*sin(lon_t), ...
                  R_earth*sin(user_lat)];

        for q=1:N_users
            rho=sat_pos-user_pos(q,:);
            dist=sqrt(sum(rho.^2,2));
            zenith=user_pos(q,:)/R_earth;
            sin_el=(rho*zenith.')./dist;
            el=asin(max(-1,min(1,sin_el)));
            nvis=nnz(el>=elevation_min);
            sum_outage_t(k,q)=sum_outage_t(k,q)+(nvis==0);
        end
    end
end

Outage_emp=sum_outage_t/N_realizations;

MeanVisible_exact=zeros(Nt,N_users);
MeanVisible_local=zeros(Nt,N_users);
Outage_exact=zeros(Nt,N_users);
Outage_local=zeros(Nt,N_users);

Nu=20000;
u_grid=linspace(0,2*pi,Nu);
phi_s=asin(sin(inc).*sin(u_grid));

for k=1:Nt
    t=time_values(k);
    fU_t=abs(cos(u_grid-omega_sat*t))/4;

    for q=1:N_users
        phi_u=user_lat(q);

        denominator=cos(phi_u).*cos(phi_s);
        numerator=cos(psi_max)-sin(phi_u).*sin(phi_s);
        qq=numerator./denominator;

        fraction=zeros(size(u_grid));
        full=qq<=-1; none=qq>=1; partial=~(full|none);
        fraction(full)=1;
        fraction(none)=0;
        fraction(partial)=acos(qq(partial))/pi;

        singular=abs(denominator)<1e-14;
        if any(singular)
            angular_distance=abs(phi_s(singular)-phi_u);
            fraction(singular)=double(angular_distance<=psi_max);
        end

        p_vis_exact=trapz(u_grid,fU_t.*fraction);
        MeanVisible_exact(k,q)=N_mean*p_vis_exact;

        fphi_local=latitude_pdf_spatial(phi_u,t,inc,omega_sat);
        rho_sat_local=N_mean*fphi_local/(2*pi*R_earth^2*max(cos(phi_u),eps));
        MeanVisible_local(k,q)=rho_sat_local*A_visible;

        Outage_exact(k,q)=exp(-MeanVisible_exact(k,q));
        Outage_local(k,q)=exp(-MeanVisible_local(k,q));
    end
end

%% ============================================================
% DIAGNOSTIC DES ERREURS
%% ============================================================

% Erreur absolue : toujours bien definie, meme lorsque
% l'outage empirique vaut exactement zero.
err_abs = abs(Outage_emp-Outage_exact);

% Avec N_realizations realisations, une probabilite inferieure a
% environ 1/N_realizations est tres mal resolue empiriquement.
% On ne calcule donc l'erreur relative que dans les cases pour
% lesquelles l'outage empirique est suffisamment non nul.
threshold_outage = 1/N_realizations;

valid_rel = Outage_emp > threshold_outage;

err_rel = NaN(size(Outage_emp));
err_rel(valid_rel) = ...
    abs(Outage_emp(valid_rel)-Outage_exact(valid_rel)) ...
    ./ Outage_emp(valid_rel);

fprintf('\n============================================================\n');
fprintf('DIAGNOSTIC OUTAGE - UNIFORMITE SPATIALE\n');
fprintf('============================================================\n');

fprintf('Nombre de realisations                    : %d\n', ...
    N_realizations);

fprintf('Seuil pour erreur relative               : %.6f\n', ...
    threshold_outage);

fprintf('Erreur absolue moyenne                   : %.6f\n', ...
    mean(err_abs,'all','omitnan'));

fprintf('Erreur absolue maximale                  : %.6f\n', ...
    max(err_abs,[],'all'));

fprintf('Erreur relative moyenne (zones valides)  : %.2f %%\n', ...
    100*mean(err_rel,'all','omitnan'));

fprintf('Nombre de cases utilisees pour err. rel. : %d / %d\n', ...
    nnz(valid_rel),numel(valid_rel));

fprintf('============================================================\n');

figure;
imagesc(user_lat_deg,time_values,Outage_emp);
set(gca,'YDir','normal'); colorbar;
xlabel('Latitude utilisateur (deg)'); ylabel('Temps (s)');
title('P_{out} empirique - uniformite spatiale');

figure;
imagesc(user_lat_deg,time_values,Outage_exact);
set(gca,'YDir','normal'); colorbar;
xlabel('Latitude utilisateur (deg)'); ylabel('Temps (s)');
title('P_{out} theorique - uniformite spatiale');

figure;
imagesc(user_lat_deg,time_values,err_abs);
set(gca,'YDir','normal');
colorbar;
xlabel('Latitude utilisateur (deg)');
ylabel('Temps (s)');
title('|P_{out}^{emp}-P_{out}^{th}|');

save('p_outage_results.mat');

function u=sample_u_spatial(N)
    s=2*rand(N,1)-1;
    a=asin(s);
    branch=rand(N,1)<0.5;
    u=a;
    u(~branch)=pi-a(~branch);
    u=mod(u,2*pi);
end

function fphi=latitude_pdf_spatial(phi,t,inc,omega)
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
