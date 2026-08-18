%% beta0_th_cluster_ppp_sphere_importance_adaptive.m
% Approximation theorique de E[beta_0] pour un PPP homogene sur une sphere
% par expansion en tailles de composantes connexes.
%
% VERSION CORRIGEE POUR LES GRANDES TAILLES s
% ============================================
%
% Probleme de la version precedente :
% ------------------------------------
% Le sampler par arbres de Cayley est exact, mais il sur-echantillonne
% fortement les configurations tres denses. En effet, une configuration G
% est generee autant de fois qu'elle possede d'arbres couvrants tau(G).
% Le poids contient donc 1/tau(G), ce qui provoque une chute de l'ESS.
%
% Correction utilisee ici :
% --------------------------
% On conserve le principe exact des arbres couvrants, mais on modifie la
% loi geometrique le long des aretes de l'arbre afin de favoriser un peu
% les distances proches du bord de la calotte de connexion.
%
% Pour une arete parent-enfant, on definit
%
%   u = (1 - cos(theta)) / (1 - cos(alpha_max))  dans [0,1].
%
% Sous la loi uniforme dans une calotte, u ~ Uniforme[0,1].
%
% On utilise la densite surfacique relative
%
%   h(u) = (1-rho) + rho * a * u^(a-1),
%
% qui est normalisee car integral_0^1 h(u) du = 1.
%
% - rho = 0 redonne EXACTEMENT l'ancien sampler uniforme.
% - rho > 0 favorise les aretes plus longues, donc des composantes moins
%   compactes et en general moins riches en arbres couvrants.
%
% Pour un graphe geometrique G, la densite marginale de proposition devient
%
%   q(X) =
%      tau_h(G)
%      -----------------------------------------
%      s^(s-2) * A_cap^(s-1)
%
% avec
%
%   tau_h(G) = somme_{T arbre couvrant de G} produit_{e in T} h(u_e).
%
% tau_h(G) est calcule EXACTEMENT par le theoreme de Matrix-Tree pondere.
%
% Le poids d'importance devient donc
%
%   exp(-lambda*A_union) / tau_h(G)
%
% et
%
%   E[C_s] =
%       lambda^s/s! * A_sphere * s^(s-2) * A_cap^(s-1)
%       * E_q[ exp(-lambda*A_union) / tau_h(G) ].
%
% Le code travaille en LOG pour eviter les overflow/underflow et augmente
% automatiquement le nombre de tirages jusqu'a atteindre un ESS cible
% ou un nombre maximal de tirages.
%
% IMPORTANT :
% -----------
% Le Monte-Carlo sert uniquement a calculer l'integrale theorique.
% Il ne simule pas le PPP lui-meme.

clear; clc; close all;
rng(1);

%% ===================== PARAMETRES =====================

R_earth = 6371;          % km
h       = 550;           % km
R       = R_earth + h;   % rayon orbital [km]

lambda  = 4e-7;          % intensite PPP [sat/km^2]
d_max   = 1500;          % distance max de lien chordale [km]

Smax    = 25;            % pousser la serie au-dela de 15

% ----- quadrature de A_union -----
Nprobe = 5000;

% ----- echantillonnage adaptatif -----
Nsamp_min = 6000;        % minimum par taille s
batch_size = 3000;       % ajouts successifs
Nsamp_max = 60000;       % plafond par taille s
ESS_target = 500;        % objectif absolu
ESS_ratio_target = 0.03; % ou 3 % des tirages, si plus exigeant

% ----- biais geometrique DEFENSIF -----
% h(u) = (1-rho) + rho*a*u^(a-1)
%
% Valeurs conseillees :
% rho = 0.60 a 0.85
% a   = 2 a 4
%
% Le terme (1-rho)>0 garantit que toutes les configurations restent
% accessibles et evite des poids pathologiques.
rho_tilt = 0.75;
beta_shape = 3.0;

% Pas de +1 heuristique par defaut.
add_giant_component = false;

% ----- comparaison empirique -----
% Fichier produit par analysis_temp.m / analysis_temp_E_sur_beta0.m.
% Il doit contenir au minimum :
%   mean_component_count_by_size
% et idealement :
%   beta0, component_size_counts, lambda, R, dmax
empirical_mat_file = 'analysis_temp_results.mat';

% Nombre maximal de tailles affichees pour la comparaison.
% On utilisera min(Smax, longueur disponible dans le .mat).
compare_with_empirical = true;


%% ===================== DONNEES EMPIRIQUES =====================

emp = struct();
has_empirical = false;

if compare_with_empirical
    if isfile(empirical_mat_file)
        emp = load(empirical_mat_file);

        if isfield(emp,'mean_component_count_by_size')
            has_empirical = true;
            fprintf('Fichier empirique charge : %s\n', empirical_mat_file);

            % Diagnostics simples de coherence des parametres.
            if isfield(emp,'lambda')
                rel = abs(emp.lambda-lambda)/max(abs(lambda),eps);
                if rel > 1e-10
                    warning('lambda theorique = %.6e, lambda empirique = %.6e.', ...
                        lambda, emp.lambda);
                end
            end

            if isfield(emp,'R')
                rel = abs(emp.R-R)/max(abs(R),eps);
                if rel > 1e-10
                    warning('R theorique = %.6f, R empirique = %.6f.', R, emp.R);
                end
            end

            if isfield(emp,'dmax')
                rel = abs(emp.dmax-d_max)/max(abs(d_max),eps);
                if rel > 1e-10
                    warning('d_max theorique = %.6f, dmax empirique = %.6f.', ...
                        d_max, emp.dmax);
                end
            end
        else
            warning(['Le fichier %s ne contient pas ' ...
                     'mean_component_count_by_size. Comparaison ignoree.'], ...
                     empirical_mat_file);
        end
    else
        warning(['Fichier empirique %s introuvable. ' ...
                 'Le calcul theorique sera execute sans comparaison.'], ...
                 empirical_mat_file);
    end
end

%% ===================== GEOMETRIE =====================

A_sphere = 4*pi*R^2;

alpha_max = 2*asin(min(1, d_max/(2*R)));
cos_alpha = cos(alpha_max);
one_minus_cos_alpha = 1-cos_alpha;

A_cap = 2*pi*R^2*one_minus_cos_alpha;

N_mean = lambda*A_sphere;
k_mean = lambda*A_cap;

fprintf('============================================================\n');
fprintf('PPP homogene sur sphere - IS pondere + adaptatif\n');
fprintf('R                 = %.3f km\n', R);
fprintf('A_sphere          = %.6e km^2\n', A_sphere);
fprintf('lambda            = %.6e sat/km^2\n', lambda);
fprintf('E[N]              = %.3f\n', N_mean);
fprintf('d_max             = %.3f km\n', d_max);
fprintf('alpha_max         = %.6f rad = %.3f deg\n', ...
    alpha_max, rad2deg(alpha_max));
fprintf('A_cap             = %.6e km^2\n', A_cap);
fprintf('degre moyen PPP   = %.4f\n', k_mean);
fprintf('Smax              = %d\n', Smax);
fprintf('Nprobe            = %d\n', Nprobe);
fprintf('rho_tilt          = %.3f\n', rho_tilt);
fprintf('beta_shape        = %.3f\n', beta_shape);
fprintf('Nsamp min/max     = %d / %d\n', Nsamp_min, Nsamp_max);
fprintf('ESS cible         = max(%g, %.1f%% des tirages)\n', ...
    ESS_target, 100*ESS_ratio_target);
fprintf('============================================================\n\n');

%% ===== POINTS QUASI-UNIFORMES POUR A_union =====

probe_unit = fibonacci_sphere(Nprobe);

%% ===================== TABLEAUX =====================

EC          = zeros(Smax,1);
stderr_EC   = zeros(Smax,1);
EN_by_size  = zeros(Smax,1);

ess         = nan(Smax,1);
Nsamp_used  = zeros(Smax,1);
mean_Aunion = nan(Smax,1);
mean_log10_tauh = nan(Smax,1);

%% ===================== s = 1 EXACT =====================

EC(1) = lambda*A_sphere*exp(-lambda*A_cap);
stderr_EC(1) = 0;
EN_by_size(1) = EC(1);

ess(1) = Inf;
Nsamp_used(1) = 0;
mean_Aunion(1) = A_cap;
mean_log10_tauh(1) = 0;

fprintf('s =  1 : E[C_1] = %.6f (exact)\n', EC(1));

%% ===================== s >= 2 =====================

x1 = [0 0 1];

for s = 2:Smax

    % Prefacteur theorique, en logarithme :
    %
    % lambda^s/s! * A_sphere * s^(s-2) * A_cap^(s-1)
    log_prefactor = ...
          s*log(lambda) ...
        - gammaln(s+1) ...
        + log(A_sphere) ...
        + (s-2)*log(s) ...
        + (s-1)*log(A_cap);

    % Stockage dynamique des poids LOG.
    logw_all = zeros(Nsamp_max,1);
    Aunion_all = zeros(Nsamp_max,1);
    logtau_all = zeros(Nsamp_max,1);

    n_done = 0;
    current_ess = 0;

    while n_done < Nsamp_max

        if n_done == 0
            n_batch = Nsamp_min;
        else
            n_batch = min(batch_size, Nsamp_max-n_done);
        end

        if n_batch <= 0
            break;
        end

        idx = n_done + (1:n_batch);

        for mm = 1:n_batch

            %% 1) arbre etiquete uniforme de Cayley
            edges = random_labeled_tree_prufer(s);

            %% 2) orientation depuis la racine 1
            [parent, order] = root_tree(edges, s, 1);

            %% 3) generation geometrique avec densite h(u)/A_cap
            X = zeros(s,3);
            X(1,:) = x1;

            for kk = 2:s
                v = order(kk);
                p = parent(v);

                X(v,:) = sample_cap_around_tilted( ...
                    X(p,:), alpha_max, rho_tilt, beta_shape);
            end

            %% 4) graphe geometrique induit
            Ddot = max(-1,min(1, X*X.'));
            Adj = (Ddot >= cos_alpha - 1e-12);
            Adj(1:s+1:end) = false;
            Adj = Adj | Adj.';

            if ~is_connected_adj(Adj)
                error('Le sampler par arbre a produit un graphe non connexe.');
            end

            %% 5) poids h_ij pour TOUTES les aretes du graphe
            %
            % u_ij est uniforme [0,1] sous la mesure surfacique uniforme
            % conditionnee au fait que ij est une arete.
            U = (1-Ddot) / one_minus_cos_alpha;
            U = max(0,min(1,U));

            H = zeros(s,s);

            edge_mask = Adj;
            H(edge_mask) = ...
                (1-rho_tilt) ...
                + rho_tilt*beta_shape .* ...
                  U(edge_mask).^(beta_shape-1);

            H = (H+H.')/2;
            H(1:s+1:end) = 0;

            %% 6) log tau_h(G) par Matrix-Tree pondere
            log_tau_h = log_weighted_spanning_tree_sum(H);

            %% 7) aire union des calottes
            maxdot = max(probe_unit * X.', [], 2);
            frac_union = mean(maxdot >= cos_alpha);
            A_union = A_sphere * frac_union;

            %% 8) poids LOG
            % log w = -lambda*A_union - log(tau_h)
            j = idx(mm);
            logw_all(j) = -lambda*A_union - log_tau_h;
            Aunion_all(j) = A_union;
            logtau_all(j) = log_tau_h;
        end

        n_done = n_done + n_batch;

        lw = logw_all(1:n_done);

        % ESS calcule proprement en log :
        %
        % ESS = (sum w)^2 / sum(w^2)
        log_sum_w  = logsumexp_vec(lw);
        log_sum_w2 = logsumexp_vec(2*lw);
        current_ess = exp(2*log_sum_w - log_sum_w2);

        target_now = max(ESS_target, ESS_ratio_target*n_done);

        if current_ess >= target_now && n_done >= Nsamp_min
            break;
        end
    end

    %% ----- moyenne et erreur standard en log/scale stable -----

    lw = logw_all(1:n_done);

    shift = max(lw);
    ws = exp(lw-shift);

    mean_scaled = mean(ws);
    std_scaled  = std(ws);

    log_mean_w = shift + log(mean_scaled);

    EC(s) = exp(log_prefactor + log_mean_w);

    % stderr = prefactor * exp(shift) * std(ws)/sqrt(n)
    stderr_EC(s) = ...
        exp(log_prefactor + shift) * std_scaled/sqrt(n_done);

    EN_by_size(s) = s*EC(s);

    ess(s) = current_ess;
    Nsamp_used(s) = n_done;
    mean_Aunion(s) = mean(Aunion_all(1:n_done));
    mean_log10_tauh(s) = mean(logtau_all(1:n_done))/log(10);

    fprintf(['s = %2d : E[C_%d] = %10.6f +/- %-9.3g | ' ...
             'sE[C_s] = %9.5f | ESS = %7.1f/%-6d | ' ...
             '<log10 tau_h> = %7.2f\n'], ...
        s, s, EC(s), stderr_EC(s), EN_by_size(s), ...
        ess(s), n_done, mean_log10_tauh(s));

    % Avertissement si le plafond est atteint avec ESS encore faible.
    if n_done >= Nsamp_max
        target_final = max(ESS_target, ESS_ratio_target*n_done);
        if current_ess < target_final
            fprintf(['         ATTENTION : ESS cible non atteinte. ' ...
                     'Augmenter Nsamp_max ou ajuster rho_tilt/beta_shape.\n']);
        end
    end
end

%% ===================== BETA_0 =====================

beta0_partial = cumsum(EC);
N_partial = cumsum(EN_by_size);

beta0_th_trunc = beta0_partial(end);
N_small_th = N_partial(end);
N_remaining = max(0,N_mean-N_small_th);

giant_correction = 0;
if add_giant_component && N_remaining > 0.5
    giant_correction = 1;
end

beta0_th_plus_giant = beta0_th_trunc + giant_correction;

fprintf('\n============================================================\n');
fprintf('RESULTATS\n');
fprintf('Somme des composantes s <= %d :\n', Smax);
fprintf('  beta0_th_trunc       = %.6f\n', beta0_th_trunc);
fprintf('  sommets expliques    = %.6f / %.6f (%.2f %%)\n', ...
    N_small_th, N_mean, 100*N_small_th/N_mean);
fprintf('  sommets restants     = %.6f\n', N_remaining);

if add_giant_component
    fprintf('  correction geante    = %d\n', giant_correction);
    fprintf('  beta0_th_plus_giant  = %.6f\n', beta0_th_plus_giant);
else
    fprintf('  correction geante DESACTIVEE\n');
end

fprintf('============================================================\n');

%% ===================== GRAPHIQUES =====================

svec = (1:Smax).';

figure;
errorbar(svec, EC, stderr_EC, 'o-', 'LineWidth', 1.4);
grid on;
xlabel('Taille s de la composante');
ylabel('E[C_s]');
title('Composantes theoriques par taille');

figure;
bar(svec, EN_by_size);
grid on;
xlabel('Taille s de la composante');
ylabel('s E[C_s]');
title('Sommets portes par les composantes de taille s');

figure;
plot(svec, beta0_partial, 'o-', 'LineWidth', 1.4);
grid on;
xlabel('s_{max}');
ylabel('\Sigma_{s=1}^{s_{max}} E[C_s]');
title('Convergence de l''expansion theorique de \beta_0');

figure;
plot(svec, N_partial, 'o-', 'LineWidth', 1.4);
hold on;
yline(N_mean, '--', 'E[N]');
grid on;
xlabel('s_{max}');
ylabel('\Sigma_{s=1}^{s_{max}} sE[C_s]');
title('Masse de sommets expliquee');

figure;
yyaxis left
plot(svec(2:end), ess(2:end)./Nsamp_used(2:end), ...
    'o-', 'LineWidth', 1.4);
ylabel('ESS / N_{samp}');
ylim([0 1]);

yyaxis right
plot(svec(2:end), Nsamp_used(2:end), ...
    's-', 'LineWidth', 1.2);
ylabel('N_{samp} utilise');

grid on;
xlabel('Taille s');
title('Qualite et cout de l''importance sampling adaptatif');


%% ===================== COMPARAISON THEORIE / EMPIRIQUE =====================

if has_empirical

    emp_counts_all = emp.mean_component_count_by_size(:);

    Scompare = min([Smax, numel(emp_counts_all)]);
    ss = (1:Scompare).';

    EC_emp = emp_counts_all(1:Scompare);
    EC_th  = EC(1:Scompare);
    SE_th  = stderr_EC(1:Scompare);

    % Erreur empirique temporelle, si component_size_counts est disponible.
    SE_emp = nan(Scompare,1);
    if isfield(emp,'component_size_counts') && ~isempty(emp.component_size_counts)
        CSC = emp.component_size_counts;

        % Convention attendue : lignes = temps, colonnes = tailles.
        if size(CSC,2) >= Scompare
            SE_emp = std(CSC(:,1:Scompare),0,1).' / sqrt(size(CSC,1));
        elseif size(CSC,1) >= Scompare
            % Secours si la matrice a ete sauvegardee transposee.
            SE_emp = std(CSC(1:Scompare,:),0,2) / sqrt(size(CSC,2));
        end
    end

    abs_err = EC_th - EC_emp;
    rel_err = abs_err ./ max(EC_emp, eps);

    beta0_emp_partial = cumsum(EC_emp);
    beta0_th_compare  = cumsum(EC_th);

    if isfield(emp,'beta0')
        beta0_emp_mean = mean(emp.beta0(:));
    else
        beta0_emp_mean = sum(emp_counts_all);
    end

    fprintf('\n============================================================\n');
    fprintf('COMPARAISON THEORIE / EMPIRIQUE PAR TAILLE\n');
    fprintf('beta0 empirique moyen total = %.6f\n', beta0_emp_mean);
    fprintf('Comparaison jusqu''a s = %d\n', Scompare);
    fprintf('------------------------------------------------------------\n');
    fprintf(' s       C_s emp       C_s th       erreur rel.      ESS\n');
    fprintf('------------------------------------------------------------\n');

    for s = 1:Scompare
        if s == 1
            ess_print = Inf;
        else
            ess_print = ess(s);
        end

        fprintf('%2d   %11.6f   %11.6f   %+10.2f %%   %8.1f\n', ...
            s, EC_emp(s), EC_th(s), 100*rel_err(s), ess_print);
    end

    fprintf('------------------------------------------------------------\n');
    fprintf('Somme emp s<=%d = %.6f\n', Scompare, sum(EC_emp));
    fprintf('Somme th  s<=%d = %.6f\n', Scompare, sum(EC_th));
    fprintf('Ecart cumulatif  = %.6f\n', sum(EC_th)-sum(EC_emp));
    fprintf('============================================================\n');

    % ---- Figure 1 : C_s emp vs theorique ----
    figure;
    hold on;

    if all(isnan(SE_emp))
        plot(ss, EC_emp, 'o-', 'LineWidth', 1.5, ...
            'DisplayName', 'Empirique');
    else
        errorbar(ss, EC_emp, SE_emp, 'o-', 'LineWidth', 1.5, ...
            'DisplayName', 'Empirique');
    end

    errorbar(ss, EC_th, SE_th, 's-', 'LineWidth', 1.5, ...
        'DisplayName', 'Theorie PPP');

    grid on;
    xlabel('Taille s de la composante');
    ylabel('E[C_s]');
    title('Nombre moyen de composantes par taille : theorie vs empirique');
    legend('Location','best');

    % ---- Figure 2 : echelle logarithmique ----
    figure;
    hold on;
    semilogy(ss, max(EC_emp,eps), 'o-', 'LineWidth', 1.5, ...
        'DisplayName', 'Empirique');
    semilogy(ss, max(EC_th,eps), 's-', 'LineWidth', 1.5, ...
        'DisplayName', 'Theorie PPP');
    grid on;
    xlabel('Taille s de la composante');
    ylabel('E[C_s] (echelle log)');
    title('Queue de la distribution des tailles de composantes');
    legend('Location','best');

    % ---- Figure 3 : erreur relative ----
    figure;
    bar(ss, 100*rel_err);
    yline(0,'--');
    grid on;
    xlabel('Taille s de la composante');
    ylabel('Erreur relative (%)');
    title('(Theorie - empirique) / empirique');

    % ---- Figure 4 : beta0 partiel ----
    figure;
    hold on;
    plot(ss, beta0_emp_partial, 'o-', 'LineWidth', 1.5, ...
        'DisplayName', '\Sigma_{j\leq s} C_j empirique');
    plot(ss, beta0_th_compare, 's-', 'LineWidth', 1.5, ...
        'DisplayName', '\Sigma_{j\leq s} E[C_j] theorique');
    yline(beta0_emp_mean, '--', ...
        sprintf('\\beta_0 emp moyen = %.2f', beta0_emp_mean), ...
        'DisplayName', '\beta_0 empirique total');
    grid on;
    xlabel('Taille maximale s');
    ylabel('Contribution cumulee a \beta_0');
    title('Construction progressive de \beta_0');
    legend('Location','best');

    % ---- Figure 5 : nombre de satellites par classe ----
    figure;
    hold on;
    plot(ss, ss.*EC_emp, 'o-', 'LineWidth', 1.5, ...
        'DisplayName', 's C_s empirique');
    plot(ss, ss.*EC_th, 's-', 'LineWidth', 1.5, ...
        'DisplayName', 's E[C_s] theorique');
    grid on;
    xlabel('Taille s de la composante');
    ylabel('Nombre moyen de satellites');
    title('Masse de sommets par taille de composante');
    legend('Location','best');

else
    Scompare = 0;
    EC_emp = [];
    SE_emp = [];
    rel_err = [];
    beta0_emp_mean = NaN;
end

%% ===================== SAUVEGARDE =====================

results = struct();

results.lambda = lambda;
results.R = R;
results.d_max = d_max;
results.alpha_max = alpha_max;
results.A_sphere = A_sphere;
results.A_cap = A_cap;
results.N_mean = N_mean;
results.k_mean = k_mean;

results.Smax = Smax;
results.Nprobe = Nprobe;

results.Nsamp_min = Nsamp_min;
results.Nsamp_max = Nsamp_max;
results.batch_size = batch_size;
results.ESS_target = ESS_target;
results.ESS_ratio_target = ESS_ratio_target;

results.rho_tilt = rho_tilt;
results.beta_shape = beta_shape;

results.EC = EC;
results.stderr_EC = stderr_EC;
results.EN_by_size = EN_by_size;

results.ess = ess;
results.Nsamp_used = Nsamp_used;
results.mean_Aunion = mean_Aunion;
results.mean_log10_tauh = mean_log10_tauh;

results.beta0_partial = beta0_partial;
results.N_partial = N_partial;

results.beta0_th_trunc = beta0_th_trunc;
results.N_small_th = N_small_th;
results.N_remaining = N_remaining;
results.giant_correction = giant_correction;
results.beta0_th_plus_giant = beta0_th_plus_giant;

results.has_empirical = has_empirical;
results.empirical_mat_file = empirical_mat_file;
results.Scompare = Scompare;
results.EC_emp = EC_emp;
results.SE_emp = SE_emp;
results.relative_error_EC = rel_err;
results.beta0_emp_mean = beta0_emp_mean;

save('beta0_th_results.mat', ...
    '-struct', 'results');

%% ===================== FONCTIONS =====================

function edges = random_labeled_tree_prufer(n)
% Arbre etiquete uniforme parmi les n^(n-2) arbres de Cayley.

    if n == 2
        edges = [1 2];
        return;
    end

    P = randi(n,n-2,1);

    deg = ones(n,1);
    for k = 1:numel(P)
        deg(P(k)) = deg(P(k))+1;
    end

    edges = zeros(n-1,2);

    for k = 1:n-2
        leaf = find(deg==1,1,'first');
        v = P(k);

        edges(k,:) = [leaf v];

        deg(leaf) = deg(leaf)-1;
        deg(v) = deg(v)-1;
    end

    last = find(deg==1);
    edges(n-1,:) = last(1:2).';
end

function [parent,order] = root_tree(edges,n,root)
% Oriente l'arbre depuis root.

    A = false(n,n);

    for k = 1:size(edges,1)
        i = edges(k,1);
        j = edges(k,2);
        A(i,j)=true;
        A(j,i)=true;
    end

    parent = zeros(n,1);
    visited = false(n,1);
    order = zeros(n,1);

    queue = zeros(n,1);
    head=1;
    tail=1;
    queue(1)=root;
    visited(root)=true;

    no=0;

    while head<=tail
        v=queue(head);
        head=head+1;

        no=no+1;
        order(no)=v;

        neigh=find(A(v,:));

        for u=neigh
            if ~visited(u)
                visited(u)=true;
                parent(u)=v;

                tail=tail+1;
                queue(tail)=u;
            end
        end
    end

    if no~=n
        error('Arbre non connexe.');
    end
end

function x = sample_cap_around_tilted(c,alpha,rho,a)
% Tire un point dans la calotte autour de c suivant
%
%   h(u) = (1-rho) + rho*a*u^(a-1)
%
% par melange :
%   avec proba 1-rho : u ~ Uniforme(0,1)
%   avec proba rho   : u ~ Beta(a,1), donc u = rand^(1/a).

    c = c(:).';
    c = c/norm(c);

    if rand < rho
        u = rand^(1/a);
    else
        u = rand;
    end

    cos_theta = 1 - u*(1-cos(alpha));
    sin_theta = sqrt(max(0,1-cos_theta^2));

    phi = 2*pi*rand;

    if abs(c(3)) < 0.9
        ref=[0 0 1];
    else
        ref=[1 0 0];
    end

    e1=cross(ref,c);
    e1=e1/norm(e1);

    e2=cross(c,e1);
    e2=e2/norm(e2);

    x = cos_theta*c ...
        + sin_theta*cos(phi)*e1 ...
        + sin_theta*sin(phi)*e2;

    x=x/norm(x);
end

function log_tau = log_weighted_spanning_tree_sum(W)
% Calcule
%
%   tau_h(G) = somme_T produit_{e in T} W_e
%
% via le Matrix-Tree theorem pondere.
%
% Le resultat est retourne en logarithme afin d'eviter les overflow.

    n=size(W,1);

    if n==1
        log_tau=0;
        return;
    end

    d=sum(W,2);
    L=diag(d)-W;

    M=L(2:end,2:end);
    M=(M+M.')/2;

    % Pour un graphe connexe avec poids strictement positifs sur les
    % aretes, le cofacteur laplacien est defini positif.
    [Rchol,p]=chol(M);

    if p==0
        log_tau = 2*sum(log(diag(Rchol)));
    else
        % Secours numerique si Chol echoue legerement.
        ev=eig(M);
        ev=real(ev);

        tol=max(1e-14,1e-12*max(abs(ev)));
        ev(ev<tol)=tol;

        log_tau=sum(log(ev));
    end
end

function tf=is_connected_adj(A)

    n=size(A,1);
    seen=false(n,1);
    queue=zeros(n,1);

    head=1;
    tail=1;

    queue(1)=1;
    seen(1)=true;

    while head<=tail
        v=queue(head);
        head=head+1;

        neigh=find(A(v,:) & ~seen.');

        for u=neigh
            seen(u)=true;
            tail=tail+1;
            queue(tail)=u;
        end
    end

    tf=all(seen);
end

function X=fibonacci_sphere(n)

    k=(0:n-1).';
    z=1-2*(k+0.5)/n;

    golden_angle=pi*(3-sqrt(5));
    phi=golden_angle*k;

    rxy=sqrt(max(0,1-z.^2));

    X=[rxy.*cos(phi),rxy.*sin(phi),z];
end

function y=logsumexp_vec(x)
% log(sum(exp(x))) stable numeriquement.

    xmax=max(x);
    y=xmax+log(sum(exp(x-xmax)));
end
