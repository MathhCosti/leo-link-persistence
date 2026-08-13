%% betti_t_phi.m
% Approximation locale et globale de beta_0(t,phi) a partir de :
%   N1_t_phi_results.mat, N2_t_phi_results.mat, N3_t_phi_results.mat.
%
% Approximation globale :
%   beta_0(t) ~= N_1(t) + N_2(t) + N_3(t) + 2.
%
% Les deux composantes macroscopiques sont reparties localement selon
% la masse residuelle de satellites n'appartenant ni aux isoles, ni aux
% dimeres, ni aux trimeres.
%
% La theorie reste approximee par :
%
%   beta_0^th(t) ~= N_1^th(t)+N_2^th(t)+N_3^th(t)+2.
%
% En revanche, la valeur empirique est calculee directement sur les
% graphes simules en comptant TOUTES les composantes connexes.
% Localement, une composante de taille s contribue pour 1/s dans la
% tranche de chacun de ses membres. Sa contribution totale vaut donc 1.

clear; clc; close all;

%% 1. Dossier et fichiers
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir = pwd; end

file_N1 = fullfile(script_dir,'N1_t_phi_results.mat');
file_N2 = fullfile(script_dir,'N2_t_phi_results.mat');
file_N3 = fullfile(script_dir,'N3_t_phi_results.mat');

for f = {file_N1,file_N2,file_N3}
    if ~isfile(f{1})
        error('Fichier introuvable : %s',f{1});
    end
end

S1 = load(file_N1);
S2 = load(file_N2);
S3 = load(file_N3);

%% 2. Verification des variables
check_fields(S1,{ ...
    'N','R','dmax','inc','omega', ...
    'time_values','phi_vals_emp','phi_edges_emp', ...
    'satellite_count_th','isolated_count_th', ...
    'isolated_count_emp','isolated_count_emp_sem'},file_N1);

check_fields(S2,{ ...
    'time_values','phi_vals_emp','component_count_th', ...
    'component_count_emp','component_count_emp_sem'},file_N2);

check_fields(S3,{ ...
    'time_values','phi_vals_emp','component_count_th', ...
    'component_count_emp','component_count_emp_sem'},file_N3);

N = double(S1.N);
R = double(S1.R);
dmax = double(S1.dmax);
inc = double(S1.inc);
omega = double(S1.omega);

time_values = double(S1.time_values(:));
phi_vals = double(S1.phi_vals_emp(:).');
phi_edges = double(S1.phi_edges_emp(:).');

Nt = numel(time_values);
Nb = numel(phi_vals);

assert_same_grid(time_values,phi_vals,S2,file_N2);
assert_same_grid(time_values,phi_vals,S3,file_N3);


%% 2.b Parametres de la mesure empirique complete

% Nombre de constellations independantes utilisees pour mesurer le vrai
% beta_0 empirique. Cette mesure recompte toutes les composantes.
n_iterations_beta0_emp = 100;

rng_seed_beta0 = 31;
rng(rng_seed_beta0);

%% 3. Petites composantes
N1_th = double(S1.isolated_count_th);
N1_emp = double(S1.isolated_count_emp);
N1_emp_sem = double(S1.isolated_count_emp_sem);

N2_th = double(S2.component_count_th);
N2_emp = double(S2.component_count_emp);
N2_emp_sem = double(S2.component_count_emp_sem);

N3_th = double(S3.component_count_th);
N3_emp = double(S3.component_count_emp);
N3_emp_sem = double(S3.component_count_emp_sem);

check_size(N1_th,Nt,Nb,'N1_th');
check_size(N2_th,Nt,Nb,'N2_th');
check_size(N3_th,Nt,Nb,'N3_th');

%% 4. Nombre de satellites par tranche
satellite_count_th = double(S1.satellite_count_th);

if isfield(S1,'satellite_count_emp_iterations')
    satellite_count_emp = squeeze(mean( ...
        double(S1.satellite_count_emp_iterations),1,'omitnan'));
elseif isfield(S1,'satellite_count_emp')
    satellite_count_emp = double(S1.satellite_count_emp);
elseif isfield(S1,'satellite_count_emp_total') && isfield(S1,'n_iterations')
    satellite_count_emp = double(S1.satellite_count_emp_total) ...
        / double(S1.n_iterations);
else
    warning(['Nombre empirique de satellites par tranche absent : ', ...
        'la repartition empirique des composantes macroscopiques ', ...
        'utilise satellite_count_th.']);
    satellite_count_emp = satellite_count_th;
end

check_size(satellite_count_th,Nt,Nb,'satellite_count_th');
check_size(satellite_count_emp,Nt,Nb,'satellite_count_emp');

%% 5. Satellites appartenant aux petites composantes
small_member_count_th = N1_th + 2*N2_th + 3*N3_th;
small_member_count_emp = N1_emp + 2*N2_emp + 3*N3_emp;

macro_member_count_th = max(satellite_count_th-small_member_count_th,0);
macro_member_count_emp = max(satellite_count_emp-small_member_count_emp,0);

%% 6. Repartition des deux composantes macroscopiques
n_macro_components = 2;
beta0_macro_th = zeros(Nt,Nb);
beta0_macro_emp = zeros(Nt,Nb);

for it = 1:Nt
    mass_th = sum(macro_member_count_th(it,:),'omitnan');
    if mass_th > 0
        beta0_macro_th(it,:) = n_macro_components ...
            * macro_member_count_th(it,:) / mass_th;
    else
        sat_mass = sum(satellite_count_th(it,:),'omitnan');
        if sat_mass > 0
            beta0_macro_th(it,:) = n_macro_components ...
                * satellite_count_th(it,:) / sat_mass;
        end
    end

    mass_emp = sum(macro_member_count_emp(it,:),'omitnan');
    if mass_emp > 0
        beta0_macro_emp(it,:) = n_macro_components ...
            * macro_member_count_emp(it,:) / mass_emp;
    else
        sat_mass = sum(satellite_count_emp(it,:),'omitnan');
        if sat_mass > 0
            beta0_macro_emp(it,:) = n_macro_components ...
                * satellite_count_emp(it,:) / sat_mass;
        end
    end
end

%% 6.b Mesure empirique complete de beta_0

beta0_emp_iterations = zeros( ...
    n_iterations_beta0_emp,Nt,Nb);

beta0_total_emp_iterations = zeros( ...
    n_iterations_beta0_emp,Nt);

for r = 1:n_iterations_beta0_emp

    [u0,Omega] = sample_initial_orbits(N,inc);

    for it = 1:Nt

        X = positions_from_orbits( ...
            u0,Omega,time_values(it),R,inc,omega);

        latitude = asin(max(min(X(:,3)/R,1),-1));
        bin_id = discretize(latitude,phi_edges);

        A = adjacency_from_positions(X,dmax);
        G = graph(sparse(A));

        component_id = conncomp(G);
        component_sizes = accumarray(component_id(:),1);

        beta0_total_emp_iterations(r,it) = ...
            numel(component_sizes);

        % Attribution locale : chaque composante de taille s contribue
        % pour 1/s dans la tranche de chacun de ses s membres.
        for c = 1:numel(component_sizes)

            members = find(component_id == c);
            s = component_sizes(c);

            for m = members(:).'
                b = bin_id(m);

                if ~isnan(b)
                    beta0_emp_iterations(r,it,b) = ...
                        beta0_emp_iterations(r,it,b)+1/s;
                end
            end
        end
    end

    fprintf('Beta0 empirique complet : realisation %d/%d\n', ...
        r,n_iterations_beta0_emp);
end

beta0_emp_true = squeeze(mean( ...
    beta0_emp_iterations,1,'omitnan'));

beta0_emp_true_std = squeeze(std( ...
    beta0_emp_iterations,0,1,'omitnan'));

beta0_emp_true_sem = ...
    beta0_emp_true_std/sqrt(n_iterations_beta0_emp);

beta0_total_emp_true = mean( ...
    beta0_total_emp_iterations,1,'omitnan').';

beta0_total_emp_true_std = std( ...
    beta0_total_emp_iterations,0,1,'omitnan').';

beta0_total_emp_true_sem = ...
    beta0_total_emp_true_std/sqrt(n_iterations_beta0_emp);

% Verification : la somme locale doit redonner exactement le nombre
% global de composantes dans chaque realisation.
beta0_total_emp_from_local_iterations = squeeze(sum( ...
    beta0_emp_iterations,3,'omitnan'));

max_local_global_error_emp = max(abs( ...
    beta0_total_emp_from_local_iterations(:) ...
    - beta0_total_emp_iterations(:)));

%% 7. beta_0 local et global
beta0_th = N1_th + N2_th + N3_th + beta0_macro_th;

% Ancienne approximation empirique conservee uniquement comme diagnostic.
beta0_emp_approx = N1_emp + N2_emp + N3_emp + beta0_macro_emp;

beta0_total_th = sum(beta0_th,2,'omitnan');
beta0_total_emp_approx = sum(beta0_emp_approx,2,'omitnan');

beta0_total_th_formula = sum(N1_th,2,'omitnan') ...
    + sum(N2_th,2,'omitnan') ...
    + sum(N3_th,2,'omitnan') ...
    + n_macro_components;

beta0_total_emp_formula = sum(N1_emp,2,'omitnan') ...
    + sum(N2_emp,2,'omitnan') ...
    + sum(N3_emp,2,'omitnan') ...
    + n_macro_components;

%% 8. Incertitudes empiriques

% Incertitude de l'ancienne approximation N1+N2+N3+2.
beta0_emp_approx_sem = sqrt( ...
    N1_emp_sem.^2+N2_emp_sem.^2+N3_emp_sem.^2);

beta0_total_emp_approx_sem = sqrt(sum( ...
    beta0_emp_approx_sem.^2,2,'omitnan'));

% Pour la mesure empirique complete, les SEM ont deja ete calculees
% directement a partir des graphes simules.

%% 9. Diagnostics
macro_sum_th = sum(beta0_macro_th,2,'omitnan');
macro_sum_emp = sum(beta0_macro_emp,2,'omitnan');

max_macro_error_th = max(abs(macro_sum_th-n_macro_components));
max_macro_error_emp = max(abs(macro_sum_emp-n_macro_components));
max_formula_error_th = max(abs(beta0_total_th-beta0_total_th_formula));
max_formula_error_emp = max(abs( ...
    beta0_total_emp_approx-beta0_total_emp_formula));

difference = beta0_emp_true-beta0_th;
valid = isfinite(difference);

rmse_local = sqrt(mean(difference(valid).^2));
mae_local = mean(abs(difference(valid)));
bias_local = mean(difference(valid));

%% 10. Figures
figure;
imagesc(time_values,rad2deg(phi_vals),beta0_th.');
axis xy; colorbar;
xlabel('Temps (s)');
ylabel('Latitude \phi (deg)');
title('\beta_0^{th}(t,\phi) : isoles + dimeres + trimeres + 2 macro');

figure;
imagesc(time_values,rad2deg(phi_vals),beta0_emp_true.');
axis xy; colorbar;
xlabel('Temps (s)');
ylabel('Latitude \phi (deg)');
title('\beta_0^{emp}(t,\phi) : toutes les composantes');

figure;
imagesc(time_values,rad2deg(phi_vals),difference.');
axis xy; colorbar;
xlabel('Temps (s)');
ylabel('Latitude \phi (deg)');
title('\beta_0^{emp}-\beta_0^{th}');

n_selected_times = 5;
selected_indices = unique(round(linspace(1,Nt,n_selected_times)));

figure;
tiledlayout(numel(selected_indices),1, ...
    'TileSpacing','compact','Padding','compact');

for k = 1:numel(selected_indices)
    it = selected_indices(k);
    nexttile; hold on;

    plot(rad2deg(phi_vals),beta0_th(it,:), ...
        's-','LineWidth',1.8,'DisplayName','Theorie');

    errorbar(rad2deg(phi_vals),beta0_emp_true(it,:), ...
        beta0_emp_true_sem(it,:), ...
        'o-','LineWidth',1.2, ...
        'DisplayName','Empirique complet \pm SEM');

    plot(rad2deg(phi_vals),beta0_macro_th(it,:), ...
        '--','LineWidth',1.2, ...
        'DisplayName','Part macroscopique th');

    grid on;
    ylabel('\beta_0 par tranche');
    title(sprintf('t = %.1f s',time_values(it)));

    if k == 1, legend('Location','best'); end
    if k == numel(selected_indices), xlabel('Latitude \phi (deg)'); end
end

figure; hold on;
plot(time_values,beta0_total_th, ...
    'LineWidth',2,'DisplayName','Theorie N_1+N_2+N_3+2');
plot(time_values,beta0_total_emp_true, ...
    'LineWidth',1.5, ...
    'DisplayName','Empirique complet');

plot(time_values, ...
    beta0_total_emp_true+beta0_total_emp_true_sem, ...
    ':','LineWidth',1,'DisplayName','Empirique + SEM');

plot(time_values, ...
    beta0_total_emp_true-beta0_total_emp_true_sem, ...
    ':','LineWidth',1,'DisplayName','Empirique - SEM');
grid on;
xlabel('Temps (s)');
ylabel('\beta_0(t)');
title('Comparaison de \beta_0 theorique et empirique complet');
legend('Location','best');
hold off;

figure; hold on;
plot(time_values,sum(N1_th,2,'omitnan'), ...
    'LineWidth',1.6,'DisplayName','N_1^{th}');
plot(time_values,sum(N2_th,2,'omitnan'), ...
    'LineWidth',1.6,'DisplayName','N_2^{th}');
plot(time_values,sum(N3_th,2,'omitnan'), ...
    'LineWidth',1.6,'DisplayName','N_3^{th}');
yline(n_macro_components,'--','2 composantes macroscopiques', ...
    'DisplayName','Composantes macroscopiques');
grid on;
xlabel('Temps (s)');
ylabel('Nombre de composantes');
title('Contributions a \beta_0^{th}');
legend('Location','best');
hold off;

%% 11. Affichage console
fprintf('\n');
fprintf('============================================================\n');
fprintf(' BETTI_0(t,phi) - THEORIE APPROX. / EMPIRIQUE COMPLET\n');
fprintf('============================================================\n');
fprintf('Fichier N1                         : %s\n',file_N1);
fprintf('Fichier N2                         : %s\n',file_N2);
fprintf('Fichier N3                         : %s\n',file_N3);
fprintf('N                                  : %d\n',round(N));
fprintf('Nombre d''instants                 : %d\n',Nt);
fprintf('Nombre de tranches                 : %d\n',Nb);
fprintf('Composantes macroscopiques         : %d\n',n_macro_components);
fprintf('------------------------------------------------------------\n');
fprintf('Moyenne beta0 theorique            : %.10f\n', ...
    mean(beta0_total_th,'omitnan'));
fprintf('Moyenne beta0 empirique approx.    : %.10f\n', ...
    mean(beta0_total_emp_approx,'omitnan'));
fprintf('RMSE locale                        : %.10e\n',rmse_local);
fprintf('MAE locale                         : %.10e\n',mae_local);
fprintf('Biais local emp-th                 : %.10e\n',bias_local);
fprintf('------------------------------------------------------------\n');
fprintf('Erreur max somme macro th          : %.10e\n',max_macro_error_th);
fprintf('Erreur max somme macro emp         : %.10e\n',max_macro_error_emp);
fprintf('Erreur max formule globale th      : %.10e\n',max_formula_error_th);
fprintf('Erreur max formule globale emp     : %.10e\n',max_formula_error_emp);
fprintf('============================================================\n');

%% 12. Sauvegarde
output_file = fullfile(script_dir,'betti_t_phi_results.mat');

save(output_file, ...
    'file_N1','file_N2','file_N3', ...
    'N','R','dmax','inc','omega', ...
    'Nt','Nb','time_values','phi_vals','phi_edges', ...
    'n_macro_components', ...
    'N1_th','N1_emp','N1_emp_sem', ...
    'N2_th','N2_emp','N2_emp_sem', ...
    'N3_th','N3_emp','N3_emp_sem', ...
    'satellite_count_th','satellite_count_emp', ...
    'small_member_count_th','small_member_count_emp', ...
    'macro_member_count_th','macro_member_count_emp', ...
    'beta0_macro_th','beta0_macro_emp', ...
    'beta0_th','beta0_emp_approx', ...
    'beta0_emp_approx_sem', ...
    'beta0_total_th','beta0_total_emp_approx', ...
    'beta0_total_emp_approx_sem', ...
    'n_iterations_beta0_emp','rng_seed_beta0', ...
    'beta0_emp_iterations', ...
    'beta0_emp_true','beta0_emp_true_std','beta0_emp_true_sem', ...
    'beta0_total_emp_iterations', ...
    'beta0_total_emp_true', ...
    'beta0_total_emp_true_std', ...
    'beta0_total_emp_true_sem', ...
    'beta0_total_emp_from_local_iterations', ...
    'max_local_global_error_emp', ...
    'beta0_total_th_formula','beta0_total_emp_formula', ...
    'macro_sum_th','macro_sum_emp', ...
    'max_macro_error_th','max_macro_error_emp', ...
    'max_formula_error_th','max_formula_error_emp', ...
    'difference','rmse_local','mae_local','bias_local', ...
    '-v7.3');

fprintf('Resultats sauvegardes dans %s\n',output_file);

%% Fonctions locales

function [u0,Omega] = sample_initial_orbits(N,inc)

    sin_phi0 = -sin(inc)+2*sin(inc)*rand(N,1);

    sin_u0 = sin_phi0/sin(inc);
    sin_u0 = min(max(sin_u0,-1),1);

    u_base = asin(sin_u0);
    ascending = rand(N,1)<0.5;

    u0 = zeros(N,1);
    u0(ascending) = u_base(ascending);
    u0(~ascending) = pi-u_base(~ascending);
    u0 = mod(u0,2*pi);

    Omega = 2*pi*rand(N,1);
end

function X = positions_from_orbits(u0,Omega,t,R,inc,omega)

    u = mod(u0+omega*t,2*pi);

    cO = cos(Omega);
    sO = sin(Omega);
    cu = cos(u);
    su = sin(u);

    X = R*[ ...
        cO.*cu-sO.*su*cos(inc), ...
        sO.*cu+cO.*su*cos(inc), ...
        su*sin(inc)];
end

function A = adjacency_from_positions(X,dmax)

    gram = X*X.';
    norm2 = sum(X.^2,2);

    D2 = max(norm2+norm2.'-2*gram,0);

    A = D2<=dmax^2;
    A(1:size(A,1)+1:end) = false;

    A = triu(A,1);
    A = A|A.';
end

function check_fields(S,fields,file_name)
    for k = 1:numel(fields)
        if ~isfield(S,fields{k})
            error('Le fichier %s doit contenir %s.',file_name,fields{k});
        end
    end
end

function assert_same_grid(time_ref,phi_ref,S,file_name)
    time_other = double(S.time_values(:));
    phi_other = double(S.phi_vals_emp(:).');

    if numel(time_other) ~= numel(time_ref) ...
            || max(abs(time_other-time_ref)) > 1e-9
        error('Grille temporelle incompatible dans %s.',file_name);
    end

    if numel(phi_other) ~= numel(phi_ref) ...
            || max(abs(phi_other-phi_ref)) > 1e-12
        error('Grille de latitude incompatible dans %s.',file_name);
    end
end

function check_size(X,Nt,Nb,name)
    if ~isequal(size(X),[Nt,Nb])
        error('%s doit etre de taille %d x %d.',name,Nt,Nb);
    end
end
