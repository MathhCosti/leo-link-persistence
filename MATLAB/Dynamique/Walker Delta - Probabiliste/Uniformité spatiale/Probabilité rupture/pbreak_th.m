%% pbreak_th.m
% Probabilite theorique de rupture pour le Walker-Delta spatial.
%
% Ce script applique directement la formule du LaTeX :
%
%   p_break(t) =
%       E[|E(t)|]/E[beta0(t)]
%       * 2 v_rel(t) Delta_t/(pi dmax)
%       * exp[-lambda_eff(t) A_inter(dmax)].
%
% Fichiers charges :
%   - ../Vitesse relative/vitesse_rel_temp.mat
%   - ../Nombre liens/liens_quadrature_results.mat
%   - ../lambda_eff_th_results.mat
%
% Sortie :
%   - pbreak_th_results.mat

clear; clc; close all;

%% ============================================================
%  LOCALISATION DES FICHIERS
%% ============================================================

script_dir = fileparts(mfilename('fullpath'));
vrel_file = fullfile(script_dir, '..', 'Paramètres', 'Vitesse relative', 'vitesse_rel_temp.mat');
links_file = fullfile(script_dir, '..', 'Paramètres', 'Nombre liens', 'liens_quadrature_results.mat');
lambda_file = fullfile(script_dir, '..', 'Paramètres', 'lambda_eff_th_results.mat');

if isempty(vrel_file)
    error('Fichier vitesse_rel_temp.mat introuvable.');
end
if isempty(links_file)
    error('Fichier liens_quadrature_results.mat introuvable.');
end
if isempty(lambda_file)
    error('Fichier lambda_eff_th_results.mat introuvable.');
end

fprintf('Fichier v_rel      : %s\n',vrel_file);
fprintf('Fichier liens      : %s\n',links_file);
fprintf('Fichier lambda_eff : %s\n',lambda_file);

%% ============================================================
%  CHARGEMENT DE v_rel(t)
%% ============================================================

Sv = load(vrel_file);

required_vrel = {'time_values','vrel_theory'};
for q = 1:numel(required_vrel)
    if ~isfield(Sv,required_vrel{q})
        error('Variable %s absente de %s.', ...
            required_vrel{q},vrel_file);
    end
end

time_vrel = double(Sv.time_values(:));
vrel_theory_t = double(Sv.vrel_theory(:));

if isfield(Sv,'vrel_emp')
    vrel_emp_t = double(Sv.vrel_emp(:));
else
    vrel_emp_t = NaN(size(vrel_theory_t));
end

if isfield(Sv,'dt')
    dt = double(Sv.dt);
elseif numel(time_vrel) >= 2
    dt = median(diff(time_vrel));
else
    error('Impossible de determiner dt.');
end

%% ============================================================
%  CHARGEMENT DE lambda_eff(t)
%% ============================================================

Slambda = load(lambda_file);

required_lambda = {'time_values','lambda_eff_t','dmax'};
for q = 1:numel(required_lambda)
    if ~isfield(Slambda,required_lambda{q})
        error('Variable %s absente de %s.', ...
            required_lambda{q},lambda_file);
    end
end

time_lambda = double(Slambda.time_values(:));
lambda_eff_t = double(Slambda.lambda_eff_t(:));
dmax = double(Slambda.dmax);

if isfield(Slambda,'inc_deg')
    inc_deg = double(Slambda.inc_deg);
else
    inc_deg = NaN;
end

%% ============================================================
%  CHARGEMENT DU NOMBRE THEORIQUE DE LIENS
%% ============================================================

Slinks = load(links_file);

if ~isfield(Slinks,'time_values')
    error('time_values absent de %s.',links_file);
end

time_links = double(Slinks.time_values(:));

if isfield(Slinks,'L_theory_time')
    E_edges_t = double(Slinks.L_theory_time(:));
    edges_source = 'L_theory_time';
elseif isfield(Slinks,'mean_edges')
    E_edges_t = double(Slinks.mean_edges(:));
    edges_source = 'mean_edges';
else
    error(['Aucune serie temporelle de liens trouvee dans %s. ', ...
           'Variables attendues : L_theory_time ou mean_edges.'], ...
           links_file);
end

%% ============================================================
%  CONSTRUCTION DE beta0 THEORIQUE INITIAL
%
%  L'approximation par isoles, dimeres et trimeres est utilisee
%  uniquement a l'instant initial :
%
%      beta0_initial ~= 1 + E[N1(0)] + E[N2(0)] + E[N3(0)].
%
%  Cette valeur est ensuite conservee constante sur toute la simulation :
%
%      beta0(t) = beta0_initial.
%
%  Ce choix evite la forte sous-estimation de beta0 lorsque les satellites
%  se concentrent aux limites de la bande et forment de grandes composantes.
%% ============================================================

if isfield(Slinks,'N_mean_theory')
    N = double(Slinks.N_mean_theory);
elseif isfield(Slinks,'N_all')
    N = mean(double(Slinks.N_all(:)),'omitnan');
else
    error('Impossible de determiner N depuis %s.',links_file);
end

if isfield(Slinks,'p_link_theory_time')
    p_link_t = double(Slinks.p_link_theory_time(:));
elseif isfield(Slinks,'L_theory_time')
    p_link_t = 2*double(Slinks.L_theory_time(:)) / (N*(N-1));
else
    p_link_t = 2*E_edges_t/(N*(N-1));
end

p_link_t = min(max(p_link_t,0),1);

c2_union = 1 + 3*sqrt(3)/(4*pi);
c3_conn = 1 + 3*sqrt(3)/(2*pi);
c3_union = 1.80;

N1_t = N .* max(1-p_link_t,0).^(N-1);

if N >= 2
    C2 = exp(gammaln(N+1)-gammaln(3)-gammaln(N-1));
    N2_t = C2 .* p_link_t .* ...
        max(1-c2_union*p_link_t,0).^(N-2);
else
    N2_t = zeros(size(p_link_t));
end

if N >= 3
    C3 = exp(gammaln(N+1)-gammaln(4)-gammaln(N-2));
    p_conn_3_t = min(max(c3_conn*p_link_t.^2,0),1);
    N3_t = C3 .* p_conn_3_t .* ...
        max(1-c3_union*p_link_t,0).^(N-3);
else
    N3_t = zeros(size(p_link_t));
end

% Approximation N1-N2-N3 conservee pour diagnostic.
beta0_N123_t = 1 + N1_t + N2_t + N3_t;
beta0_N123_t = min(max(beta0_N123_t,1),N);

% Valeur initiale retenue dans la formule de p_break.
beta0_initial = beta0_N123_t(1);

% beta0 est ensuite suppose constant.
beta0_t = beta0_initial * ones(size(beta0_N123_t));

beta0_source = ...
    'beta0 constant egal a 1 + E[N1(0)] + E[N2(0)] + E[N3(0)]';

%% ============================================================
%  ALIGNEMENT TEMPOREL
%% ============================================================

time_values = time_lambda;

vrel_theory_aligned = interp1( ...
    time_vrel,vrel_theory_t,time_values,'linear','extrap');

vrel_emp_aligned = interp1( ...
    time_vrel,vrel_emp_t,time_values,'linear','extrap');

E_edges_aligned = interp1( ...
    time_links,E_edges_t,time_values,'linear','extrap');

beta0_aligned = interp1( ...
    time_links,beta0_t,time_values,'linear','extrap');

p_link_aligned = interp1( ...
    time_links,p_link_t,time_values,'linear','extrap');

if numel(time_values) < 2
    error('La grille temporelle doit contenir au moins deux points.');
end

% Probabilites sur les transitions [t_k,t_{k+1}].
t_transition = time_values(1:end-1);
lambda_transition = lambda_eff_t(1:end-1);
vrel_transition = vrel_theory_aligned(1:end-1);
vrel_emp_transition = vrel_emp_aligned(1:end-1);
E_edges_transition = E_edges_aligned(1:end-1);
beta0_transition = beta0_aligned(1:end-1);
p_link_transition = p_link_aligned(1:end-1);

%% ============================================================
%  CALCUL DE p_break SELON LA FORMULE DU LATEX
%% ============================================================

[p_break_t,details_theory] = calc_p_break_th( ...
    vrel_transition, ...
    lambda_transition, ...
    E_edges_transition, ...
    beta0_transition, ...
    dmax,dt);

%% Diagnostic avec la vitesse empirique uniquement
[p_break_vrel_emp_t,details_emp] = calc_p_break_th( ...
    vrel_emp_transition, ...
    lambda_transition, ...
    E_edges_transition, ...
    beta0_transition, ...
    dmax,dt);

%% Reference constante construite avec les moyennes temporelles
vrel_mean = mean(vrel_transition,'omitnan');
lambda_eff_mean = mean(lambda_transition,'omitnan');
E_edges_mean = mean(E_edges_transition,'omitnan');
beta0_mean = mean(beta0_transition,'omitnan');

[p_break_const_t,details_const] = calc_p_break_th( ...
    vrel_mean*ones(size(vrel_transition)), ...
    lambda_eff_mean*ones(size(lambda_transition)), ...
    E_edges_mean*ones(size(E_edges_transition)), ...
    beta0_mean*ones(size(beta0_transition)), ...
    dmax,dt);

p_break_const = p_break_const_t(1);

%% ============================================================
%  AFFICHAGE
%% ============================================================

fprintf('\n=== p_break theorique spatial : formule du LaTeX ===\n');
fprintf('Inclinaison                        : %.2f deg\n',inc_deg);
fprintf('N utilise                          : %.6f\n',N);
fprintf('dmax                               : %.2f km\n',dmax);
fprintf('dt                                 : %.2f s\n',dt);
fprintf('Source de E[|E(t)|]                : %s\n',edges_source);
fprintf('Source de E[beta0(t)]              : %s\n',beta0_source);
fprintf('lambda_eff moyen                   : %.6e sat/km^2\n', ...
    lambda_eff_mean);
fprintf('E[|E(t)|] moyen                    : %.6f\n', ...
    E_edges_mean);
fprintf('beta0 initial constant             : %.6f\n', ...
    beta0_initial);
fprintf('E[beta0(t)] moyen                  : %.6f\n', ...
    beta0_mean);
fprintf('v_rel theorique moyen              : %.6f km/s\n', ...
    vrel_mean);
fprintf('p_break moyen / min / max          : %.6f / %.6f / %.6f\n', ...
    mean(p_break_t,'omitnan'), ...
    min(p_break_t,[],'omitnan'), ...
    max(p_break_t,[],'omitnan'));
fprintf('Reference constante                : %.6f\n',p_break_const);

figure;
hold on; grid on;
plot(t_transition,p_break_t,'LineWidth',1.7, ...
    'DisplayName','Formule théorique du LaTeX');
plot(t_transition,p_break_vrel_emp_t,'--','LineWidth',1.2, ...
    'DisplayName','Diagnostic avec v_{rel} empirique');
yline(mean(p_break_t,'omitnan'),':', ...
    sprintf('Moyenne = %.4f',mean(p_break_t,'omitnan')));
yline(p_break_const,'--', ...
    sprintf('Référence constante = %.4f',p_break_const));
xlabel('Temps (s)');
ylabel('p_{break}^{th}(t)');
title('Walker-Delta spatial : formule théorique de rupture');
legend('Location','best');
hold off;

figure;
hold on; grid on;
plot(t_transition,lambda_transition,'LineWidth',1.5);
xlabel('Temps (s)');
ylabel('\lambda_{eff}(t) (satellites/km^2)');
title('Intensité effective utilisée');

figure;
hold on; grid on;
plot(t_transition,E_edges_transition,'LineWidth',1.5, ...
    'DisplayName','E[|E(t)|]');
plot(t_transition,beta0_transition,'LineWidth',1.5, ...
    'DisplayName','E[\beta_0(t)]');
xlabel('Temps (s)');
ylabel('Nombre moyen');
title('Grandeurs topologiques théoriques');
legend('Location','best');


%% Comparaison de l'ancienne approximation N1-N2-N3 et de beta0 constant
beta0_N123_aligned = interp1( ...
    time_links,beta0_N123_t,time_values,'linear','extrap');

figure;
hold on; grid on;
plot(time_values,beta0_N123_aligned,'LineWidth',1.5, ...
    'DisplayName','Approximation N_1+N_2+N_3 variable');
yline(beta0_initial,'--','LineWidth',1.8, ...
    'DisplayName',sprintf('\beta_0 constant = %.3f',beta0_initial));
xlabel('Temps (s)');
ylabel('\beta_0');
title('Choix de \beta_0 constant égal à sa valeur initiale');
legend('Location','best');
hold off;


%% ============================================================
%  DIAGNOSTIC DES FACTEURS DE LA FORMULE
%
%  p_break_raw(t) =
%      [E(t)/beta0(t)]
%      * [2 v_rel(t) dt/(pi dmax)]
%      * exp[-lambda_eff(t) A_inter(dmax)]
%% ============================================================

ratio_E_beta0_t = details_theory.mean_links_per_component_t;
exp_factor_t = details_theory.p_bridge_bord_t;
q_break_link_t = details_theory.q_break_link_t;
p_break_raw_t = details_theory.p_break_raw_t;

%% 1. Terme exponentiel et lambda_eff brute
figure;
yyaxis left;
plot(t_transition,lambda_transition,'LineWidth',1.5, ...
    'DisplayName','\lambda_{eff}(t)');
ylabel('\lambda_{eff}(t) (satellites/km^2)');

yyaxis right;
plot(t_transition,exp_factor_t,'--','LineWidth',1.8, ...
    'DisplayName','exp[-\lambda_{eff}(t)A_{\cap}]');
ylabel('Facteur exponentiel');

grid on;
xlabel('Temps (s)');
title('Intensité effective et facteur exponentiel');
legend('Location','best');

%% 2. E(t), beta0(t) et ratio E/beta0
figure;
yyaxis left;
plot(t_transition,E_edges_transition,'LineWidth',1.5, ...
    'DisplayName','E[|\mathcal{E}(t)|]');
plot(t_transition,beta0_transition,'LineWidth',1.5, ...
    'DisplayName','E[\beta_0(t)]');
ylabel('Valeurs brutes');

yyaxis right;
plot(t_transition,ratio_E_beta0_t,'--','LineWidth',1.8, ...
    'DisplayName','E[|\mathcal{E}(t)|]/E[\beta_0(t)]');
ylabel('Nombre moyen de liens par composante');

grid on;
xlabel('Temps (s)');
title('Grandeurs topologiques brutes et ratio E/\beta_0');
legend('Location','best');

%% 3. Les trois facteurs multiplicatifs de p_break
figure;
hold on; grid on;
plot(t_transition,ratio_E_beta0_t,'LineWidth',1.5, ...
    'DisplayName','E/\beta_0');
plot(t_transition,q_break_link_t,'LineWidth',1.5, ...
    'DisplayName','2v_{rel}\Delta t/(\pi d_{max})');
plot(t_transition,exp_factor_t,'LineWidth',1.5, ...
    'DisplayName','exp[-\lambda_{eff}A_{\cap}]');
xlabel('Temps (s)');
ylabel('Valeur du facteur');
title('Facteurs multiplicatifs de la formule de p_{break}');
legend('Location','best');
hold off;

%% 4. Valeur brute avant saturation et valeur finale
figure;
hold on; grid on;
plot(t_transition,p_break_raw_t,'LineWidth',1.7, ...
    'DisplayName','p_{break}^{brut}(t)');
plot(t_transition,p_break_t,'--','LineWidth',1.7, ...
    'DisplayName','p_{break}(t) après saturation');
yline(1,':','Seuil de saturation');
xlabel('Temps (s)');
ylabel('Probabilité / valeur brute');
title('Effet de la saturation à 1');
legend('Location','best');
hold off;

fprintf('\n--- Diagnostic des facteurs ---\n');
fprintf('A_inter(dmax)                          : %.6e km^2\n', ...
    details_theory.A_inter_at_dmax);
fprintf('Facteur exp min / moyen / max          : %.6f / %.6f / %.6f\n', ...
    min(exp_factor_t,[],'omitnan'), ...
    mean(exp_factor_t,'omitnan'), ...
    max(exp_factor_t,[],'omitnan'));
fprintf('Ratio E/beta0 min / moyen / max        : %.6f / %.6f / %.6f\n', ...
    min(ratio_E_beta0_t,[],'omitnan'), ...
    mean(ratio_E_beta0_t,'omitnan'), ...
    max(ratio_E_beta0_t,[],'omitnan'));
fprintf('q_break lien min / moyen / max         : %.6f / %.6f / %.6f\n', ...
    min(q_break_link_t,[],'omitnan'), ...
    mean(q_break_link_t,'omitnan'), ...
    max(q_break_link_t,[],'omitnan'));
fprintf('p_break brut min / moyen / max         : %.6f / %.6f / %.6f\n', ...
    min(p_break_raw_t,[],'omitnan'), ...
    mean(p_break_raw_t,'omitnan'), ...
    max(p_break_raw_t,[],'omitnan'));
fprintf('Fraction de points saturés à 1         : %.4f\n', ...
    mean(p_break_raw_t >= 1,'omitnan'));

%% ============================================================
%  SAUVEGARDE
%% ============================================================

save('pbreak_th_results.mat', ...
    'time_values','t_transition', ...
    'vrel_theory_aligned','vrel_emp_aligned', ...
    'vrel_transition','vrel_emp_transition', ...
    'lambda_eff_t','lambda_transition','lambda_eff_mean', ...
    'E_edges_t','E_edges_aligned','E_edges_transition','E_edges_mean', ...
    'p_link_t','p_link_aligned','p_link_transition', ...
    'N1_t','N2_t','N3_t', ...
    'beta0_N123_t','beta0_N123_aligned','beta0_initial', ...
    'beta0_t','beta0_aligned','beta0_transition','beta0_mean', ...
    'N','dmax','dt','inc_deg', ...
    'p_break_t','p_break_vrel_emp_t','p_break_const', ...
    'details_theory','details_emp','details_const', ...
    'ratio_E_beta0_t','exp_factor_t','q_break_link_t','p_break_raw_t', ...
    'edges_source','beta0_source', ...
    'vrel_file','links_file','lambda_file');

fprintf('\nResultats sauvegardes dans pbreak_th_results.mat\n');

%% ============================================================
%  FONCTION LOCALE
%% ============================================================

function file = first_existing_file(candidates)
    file = '';
    for k = 1:numel(candidates)
        if isfile(candidates{k})
            file = candidates{k};
            return;
        end
    end
end
