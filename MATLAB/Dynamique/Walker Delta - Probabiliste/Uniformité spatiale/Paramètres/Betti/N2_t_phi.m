%% dimeres_t_phi_compare_geometrique.m
% Comparaison theorie / empirique du nombre local de dimeres.
%
% La fermeture Erdos-Renyi est remplacee par la formule geometrique :
%
% C_2,b^th(t)
%   = N(N-1)/2 * P(Phi_1 dans b)
%     * E[ 1_{X1~X2} (1-q_2(X1,X2))^(N-2)
%          | Phi_1 dans b ],
%
% avec
%
% q_2(X1,X2)
%   = P(X3~X1 ou X3~X2 | X1,X2).
%
% L'esperance est estimee par Monte-Carlo imbrique, en conditionnant
% X2 a etre lie a X1 afin d'eviter les evenements rares.
%
% Entree :
%   plink_t_phi_results.mat
%
% Sortie :
%   dimeres_t_phi_results.mat

clear; clc; close all;

script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir = pwd; end

input_file = fullfile(script_dir,'..', 'plink_t_phi_results.mat');
if ~isfile(input_file)
    error('Fichier introuvable : %s',input_file);
end

S = load(input_file);

required = {'N','R','dmax','inc','omega','time_theory', ...
    'phi_edges_emp','phi_vals_emp','dphi_emp'};

for k = 1:numel(required)
    if ~isfield(S,required{k})
        error('Le fichier doit contenir %s.',required{k});
    end
end

N = double(S.N);
R = double(S.R);
dmax = double(S.dmax);
inc = double(S.inc);
omega = double(S.omega);

time_values = double(S.time_theory(:).');
phi_edges_emp = double(S.phi_edges_emp(:).');
phi_vals_emp = double(S.phi_vals_emp(:).');
dphi_emp = double(S.dphi_emp(:).');

Nt = numel(time_values);
Nb = numel(phi_vals_emp);

%% Parametres
n_iterations_emp = 100;

% Theorie geometrique : reglages precision / temps de calcul.
n_time_theory = 31;
n_outer = 500;
% Pools independants :
%   - estimation de p1 ;
%   - tirage de X2 conditionnellement au lien ;
%   - estimation non biaisee de l'exclusion.
n_probe_link_est = 1200;
n_probe_link_draw = 1200;
n_probe_exclusion = 1200;

n_selected_times = 5;
rng_seed = 14;
rng(rng_seed);

theory_indices = unique(round(linspace(1,Nt,n_time_theory)));
time_theory_eval = time_values(theory_indices);
Nte = numel(theory_indices);

%% Masse theorique de chaque tranche
% Estimee par Monte-Carlo independamment de la formule de composante.
n_mass_samples = 200000;
mass_bin_th = zeros(Nte,Nb);

for it = 1:Nte
    Xmass = sample_global_positions( ...
        n_mass_samples,time_theory_eval(it),R,inc,omega);

    phi_mass = asin(max(min(Xmass(:,3)/R,1),-1));
    counts = histcounts(phi_mass,phi_edges_emp);
    mass_bin_th(it,:) = counts/n_mass_samples;
end

%% Theorie geometrique des dimeres
component_count_th_eval = zeros(Nte,Nb);
p_member_th_eval = zeros(Nte,Nb);
mean_q2_eval = nan(Nte,Nb);
valid_outer_eval = zeros(Nte,Nb);

for it = 1:Nte
    t = time_theory_eval(it);

    fprintf('Theorie dimeres : instant %d/%d\n',it,Nte);

    for b = 1:Nb
        X1_all = sample_positions_in_bin( ...
            n_outer,t,R,inc,omega, ...
            phi_edges_emp(b),phi_edges_emp(b+1));

        weights = nan(n_outer,1);
        q2_values = nan(n_outer,1);

        for s = 1:n_outer
            X1 = X1_all(s,:);

            % Pool A : estimation de p1 = P(X2~X1).
            P_est = sample_global_positions( ...
                n_probe_link_est,t,R,inc,omega);

            link_est = ...
                squared_distance_to_point(P_est,X1) <= dmax^2;

            p1 = mean(link_est);

            % Pool B independant : tirage de X2 | X2~X1.
            P_draw = sample_global_positions( ...
                n_probe_link_draw,t,R,inc,omega);

            link_draw = ...
                squared_distance_to_point(P_draw,X1) <= dmax^2;

            idx_link = find(link_draw);

            if isempty(idx_link)
                continue;
            end

            X2 = P_draw(idx_link(randi(numel(idx_link))),:);

            % Pool C independant : estimation non biaisee de
            %
            %   (1-q2)^(N-2),
            %
            % ou q2 est la probabilite qu'un satellite exterieur
            % touche X1 ou X2.
            P_excl = sample_global_positions( ...
                n_probe_exclusion,t,R,inc,omega);

            union_link = ...
                squared_distance_to_point(P_excl,X1) <= dmax^2 ...
                | squared_distance_to_point(P_excl,X2) <= dmax^2;

            K2 = nnz(union_link);
            M2 = numel(union_link);

            p_exclusion_2 = ...
                exclusion_probability_unbiased(K2,M2,N-2);

            q2 = K2/M2;

            weights(s) = p1*p_exclusion_2;
            q2_values(s) = q2;
        end

        valid = isfinite(weights);
        valid_outer_eval(it,b) = nnz(valid);

        if ~any(valid)
            component_count_th_eval(it,b) = NaN;
            continue;
        end

        mean_weight = mean(weights(valid));

        % Probabilite qu'un satellite de la tranche appartienne
        % a un dimere.
        p_member = (N-1)*mean_weight;
        p_member_th_eval(it,b) = p_member;

        % Une composante est partagee entre ses deux membres.
        component_count_th_eval(it,b) = ...
            N/2 * mass_bin_th(it,b) * p_member;

        mean_q2_eval(it,b) = mean(q2_values(valid),'omitnan');
    end
end

%% Interpolation temporelle de la theorie
component_count_th = nan(Nt,Nb);
p_member_th = nan(Nt,Nb);

for b = 1:Nb
    component_count_th(:,b) = interp1( ...
        time_theory_eval,component_count_th_eval(:,b), ...
        time_values,'pchip','extrap');

    p_member_th(:,b) = interp1( ...
        time_theory_eval,p_member_th_eval(:,b), ...
        time_values,'pchip','extrap');
end

component_count_th = max(component_count_th,0);
p_member_th = max(p_member_th,0);

%% Resimulation empirique
component_count_emp_iterations = ...
    zeros(n_iterations_emp,Nt,Nb);

for r = 1:n_iterations_emp
    [u0,Omega] = sample_initial_orbits(N,inc);

    for it = 1:Nt
        X = positions_from_orbits( ...
            u0,Omega,time_values(it),R,inc,omega);

        latitude = asin(max(min(X(:,3)/R,1),-1));
        A = adjacency_from_positions(X,dmax);

        G = graph(sparse(A));
        cid = conncomp(G);
        sizes = accumarray(cid(:),1);
        targets = find(sizes == 2);

        bin_id = discretize(latitude,phi_edges_emp);

        for c = targets(:).'
            members = find(cid == c);

            for m = members(:).'
                b = bin_id(m);
                if ~isnan(b)
                    component_count_emp_iterations(r,it,b) = ...
                        component_count_emp_iterations(r,it,b)+1/2;
                end
            end
        end
    end

    fprintf('Empirique dimeres : realisation %d/%d\n', ...
        r,n_iterations_emp);
end

component_count_emp = squeeze(mean( ...
    component_count_emp_iterations,1,'omitnan'));

component_count_emp_std = squeeze(std( ...
    component_count_emp_iterations,0,1,'omitnan'));

component_count_emp_sem = ...
    component_count_emp_std/sqrt(n_iterations_emp);

component_total_th = sum(component_count_th,2,'omitnan');

component_total_emp_iterations = squeeze(sum( ...
    component_count_emp_iterations,3,'omitnan'));

component_total_emp = mean( ...
    component_total_emp_iterations,1,'omitnan').';

component_total_emp_std = std( ...
    component_total_emp_iterations,0,1,'omitnan').';

component_total_emp_sem = ...
    component_total_emp_std/sqrt(n_iterations_emp);

%% Diagnostics
difference = component_count_emp-component_count_th;
valid = isfinite(difference);

rmse_grid = sqrt(mean(difference(valid).^2));
mae_grid = mean(abs(difference(valid)));
bias_grid = mean(difference(valid));

%% Figures
selected = unique(round(linspace(1,Nt,n_selected_times)));

figure;
tiledlayout(numel(selected),1, ...
    'TileSpacing','compact','Padding','compact');

for k = 1:numel(selected)
    it = selected(k);
    nexttile; hold on;

    plot(rad2deg(phi_vals_emp),component_count_th(it,:), ...
        's-','LineWidth',1.8,'DisplayName','Theorie geometrique');

    errorbar(rad2deg(phi_vals_emp),component_count_emp(it,:), ...
        component_count_emp_sem(it,:), ...
        'o-','LineWidth',1.2, ...
        'DisplayName','Empirique moyen \pm SEM');

    grid on;
    ylabel('Nombre de dimeres');
    title(sprintf('t = %.1f s',time_values(it)));

    if k == 1, legend('Location','best'); end
    if k == numel(selected), xlabel('Latitude \phi (deg)'); end
end

figure; hold on;
plot(time_values,component_total_th, ...
    'LineWidth',2,'DisplayName','Theorie geometrique');
plot(time_values,component_total_emp, ...
    'LineWidth',1.5,'DisplayName','Empirique moyen');
plot(time_values,component_total_emp+component_total_emp_sem, ...
    ':','DisplayName','Empirique + SEM');
plot(time_values,component_total_emp-component_total_emp_sem, ...
    ':','DisplayName','Empirique - SEM');
grid on;
xlabel('Temps (s)');
ylabel('Nombre total de dimeres');
title('Nombre total de dimeres');
legend('Location','best');
hold off;

figure;
imagesc(time_values,rad2deg(phi_vals_emp),difference.');
axis xy; colorbar;
xlabel('Temps (s)');
ylabel('Latitude \phi (deg)');
title('Ecart dimeres empiriques - theorie geometrique');

fprintf('\nDIMÈRES - MODELE GEOMETRIQUE\n');
fprintf('RMSE locale : %.6e composantes\n',rmse_grid);
fprintf('MAE locale  : %.6e composantes\n',mae_grid);
fprintf('Biais       : %.6e composantes\n',bias_grid);
fprintf('Moyenne totale empirique  : %.6f\n', ...
    mean(component_total_emp,'omitnan'));
fprintf('Moyenne totale theorique  : %.6f\n', ...
    mean(component_total_th,'omitnan'));

%% Sauvegarde
output_file = fullfile(script_dir,'N2_t_phi_results.mat');

save(output_file, ...
    'input_file','N','R','dmax','inc','omega', ...
    'time_values','Nt','Nb', ...
    'phi_edges_emp','phi_vals_emp','dphi_emp', ...
    'n_iterations_emp','n_time_theory', ...
    'n_outer','n_probe_link_est','n_probe_link_draw', ...
    'n_probe_exclusion', ...
    'theory_indices','time_theory_eval', ...
    'mass_bin_th','mean_q2_eval','valid_outer_eval', ...
    'p_member_th_eval','p_member_th', ...
    'component_count_th_eval','component_count_th', ...
    'component_count_emp_iterations', ...
    'component_count_emp','component_count_emp_std', ...
    'component_count_emp_sem', ...
    'component_total_th', ...
    'component_total_emp_iterations', ...
    'component_total_emp','component_total_emp_std', ...
    'component_total_emp_sem', ...
    'difference','rmse_grid','mae_grid','bias_grid', ...
    '-v7.3');

fprintf('Resultats sauvegardes dans %s\n',output_file);

%% Fonctions locales
function X = sample_global_positions(M,t,R,inc,omega)
    [u0,Omega] = sample_initial_orbits(M,inc);
    X = positions_from_orbits(u0,Omega,t,R,inc,omega);
end

function X = sample_positions_in_bin(M,t,R,inc,omega,phi_min,phi_max)
    X = zeros(M,3);
    n = 0;

    while n < M
        batch = max(1000,4*(M-n));
        C = sample_global_positions(batch,t,R,inc,omega);
        phi = asin(max(min(C(:,3)/R,1),-1));

        mask = phi >= phi_min & phi < phi_max;
        accepted = C(mask,:);

        take = min(size(accepted,1),M-n);
        if take > 0
            X(n+(1:take),:) = accepted(1:take,:);
            n = n+take;
        end
    end
end

function [u0,Omega] = sample_initial_orbits(M,inc)
    sin_phi0 = -sin(inc)+2*sin(inc)*rand(M,1);
    sin_u0 = sin_phi0/sin(inc);
    sin_u0 = min(max(sin_u0,-1),1);

    u_base = asin(sin_u0);
    ascending = rand(M,1) < 0.5;

    u0 = zeros(M,1);
    u0(ascending) = u_base(ascending);
    u0(~ascending) = pi-u_base(~ascending);
    u0 = mod(u0,2*pi);

    Omega = 2*pi*rand(M,1);
end

function X = positions_from_orbits(u0,Omega,t,R,inc,omega)
    u = mod(u0+omega*t,2*pi);

    cO = cos(Omega); sO = sin(Omega);
    cu = cos(u); su = sin(u);

    X = R*[ ...
        cO.*cu-sO.*su*cos(inc), ...
        sO.*cu+cO.*su*cos(inc), ...
        su*sin(inc)];
end

function D2 = squared_distance_to_point(X,x)
    D2 = sum((X-x).^2,2);
end

function p = exclusion_probability_unbiased(K,M,m)
% Estimateur non biaise de (1-q)^m lorsque K~Bin(M,q) :
%
%   E[ C(M-K,m)/C(M,m) ] = (1-q)^m.
%
% Le calcul est effectue sous forme logarithmique pour eviter les
% debordements et sous-flux numeriques.

    n_free = M-K;

    if m < 0 || M <= 0
        error('Parametres invalides pour la probabilite d''exclusion.');
    end

    if m == 0
        p = 1;
        return;
    end

    if n_free < m || M < m
        p = 0;
        return;
    end

    j = 0:(m-1);

    log_p = sum( ...
        log(n_free-j) ...
        - log(M-j));

    p = exp(log_p);
end

function A = adjacency_from_positions(X,dmax)
    g = X*X.';
    n2 = sum(X.^2,2);
    D2 = max(n2+n2.'-2*g,0);

    A = D2 <= dmax^2;
    A(1:size(A,1)+1:end) = false;
    A = triu(A,1);
    A = A|A.';
end
