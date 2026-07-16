clear; clc; close all;

%% p_link(phi) - Walker Delta a uniformite orbitale

R = 6371 + 550;       % km
inc_deg = 70;
inc = deg2rad(inc_deg);
dmax = 1500;          % km

n_phi = 200;
nQuad = 300;

eps_phi = 1e-6;
phi_vals = linspace(-inc+eps_phi, inc-eps_phi, n_phi);
p_link_phi = zeros(size(phi_vals));

alpha_max = 2*asin(min(dmax/(2*R),1));
cmax = cos(alpha_max);

[u2,w2] = gauss_legendre_interval(nQuad,0,2*pi);
u2 = u2(:); w2 = w2(:);

c2 = cos(u2);
s2 = sin(u2);
ci = cos(inc);
si = sin(inc);

for k = 1:numel(phi_vals)

    phi = phi_vals(k);

    s1 = sin(phi)/si;
    s1 = max(min(s1,1),-1);

    u1a = asin(s1);
    u1b = pi-u1a;

    p_branch = zeros(2,1);
    u1_list = [u1a,u1b];

    for b = 1:2
        u1 = u1_list(b);
        c1 = cos(u1);
        s1b = sin(u1);

        A = c1.*c2 + ci^2.*s1b.*s2;
        B = ci.*(s1b.*c2 - c1.*s2);
        C = si^2.*s1b.*s2;

        rho = sqrt(A.^2+B.^2);

        g = zeros(size(rho));
        mask = rho > 1e-14;

        q = zeros(size(rho));
        q(mask) = (cmax-C(mask))./rho(mask);

        g(mask & q <= -1) = 1;

        middle = mask & q > -1 & q < 1;
        g(middle) = acos(q(middle))/pi;

        g(~mask) = double(C(~mask) >= cmax);

        p_branch(b) = sum(w2.*g)/(2*pi);
    end

    p_link_phi(k) = mean(p_branch);
end

f_phi = cos(phi_vals) ./ ...
    (pi*sqrt(max(si^2-sin(phi_vals).^2,eps)));

f_phi = f_phi/trapz(phi_vals,f_phi);

p_link_global_from_phi = trapz(phi_vals,p_link_phi.*f_phi);
p_link_sphere = (1-cos(alpha_max))/2;

figure;
yyaxis left
plot(rad2deg(phi_vals),p_link_phi,'LineWidth',2);
ylabel('p_{link}(\phi)');

yyaxis right
plot(rad2deg(phi_vals),f_phi,'--','LineWidth',1.8);
ylabel('f_\phi(\phi)');

grid on;
xlabel('Latitude \phi (deg)');
title('Probabilite locale de lien et densite de latitude');
legend('p_{link}(\phi)','f_\phi(\phi)','Location','best');

fprintf('\n=== p_link(phi) Walker Delta ===\n');
fprintf('Inclinaison               : %.2f deg\n',inc_deg);
fprintf('dmax                      : %.2f km\n',dmax);
fprintf('p_link global reconstruit : %.10f\n',p_link_global_from_phi);
fprintf('p_link uniforme sphere    : %.10f\n',p_link_sphere);

save('plink_phi_walker_delta_results.mat', ...
    'R','inc_deg','inc','dmax','alpha_max','phi_vals', ...
    'p_link_phi','f_phi','p_link_global_from_phi','p_link_sphere','nQuad');

function [x,w] = gauss_legendre_interval(n,a,b)
k = (1:n-1)';
beta = k./sqrt(4*k.^2-1);
J = diag(beta,1)+diag(beta,-1);
[V,D] = eig(J);
x0 = diag(D);
[x0,idx] = sort(x0);
V = V(:,idx);
w0 = 2*(V(1,:).^2)';
x = (b-a)/2*x0+(a+b)/2;
w = (b-a)/2*w0;
end
