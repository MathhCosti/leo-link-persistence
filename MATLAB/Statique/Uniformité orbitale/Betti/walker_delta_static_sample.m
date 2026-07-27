function [positions, Omega, u0, latitude] = walker_delta_static_sample(N, R, inc_rad)
Omega = 2*pi*rand(N,1);
u0 = 2*pi*rand(N,1);

cO = cos(Omega); sO = sin(Omega);
cu = cos(u0);    su = sin(u0);
ci = cos(inc_rad); si = sin(inc_rad);

x = R .* (cO.*cu - sO.*su.*ci);
y = R .* (sO.*cu + cO.*su.*ci);
z = R .* (su.*si);

positions = [x y z];
latitude = asin(max(min(z./R,1),-1));
end
