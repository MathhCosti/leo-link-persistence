clear; clc; close all;

%% ============================================================
%  p_break THEORIQUE - WALKER DELTA
%
%  Etape 1 : probabilite conditionnelle de rupture d'un pont
%
%  Pour un lien quelconque :
%
%      q_break^lien = 2 v_rad^Delta Delta_t / d_max.
%
%  Mais un pont n'a pas la meme distribution de longueur qu'un lien
%  quelconque. Sous l'approximation locale "absence de voisin commun" :
%
%      P(pont | D=d, lien) ~= exp[-lambda A_inter(d)].
%
%  Par conditionnement :
%
%      f_{D|pont}(d_max)
%        = (2/d_max) *
%          P(pont | D=d_max,lien) / chi_bridge.
%
%  Ainsi :
%
%      q_break^pont
%        = q_break^lien *
%          P(pont | D=d_max,lien) / chi_bridge.
%
%  Etape 2 : passage au niveau d'une composante
%
%      E[B_C] = E[L_C] chi_bridge,
%
%      p_break^Delta
%        = 1 - (1-q_break^pont)^(E[B_C]).
%
%  La version lineaire se simplifie en :
%
%      p_break^Delta
%        ~= E[L_C] q_break^lien
%            P(pont | D=d_max,lien).
%% ============================================================

%% Parametres physiques
R_earth = 6371;          % km
h = 550;                 % km
R = R_earth + h;         % km

mu_earth = 398600;  % km^3/s^2

v_orb = sqrt(mu_earth/R);    % km/s
v_rel = v_orb/sqrt(2);       % km/s
v_rad = v_rel/pi;            % km/s

%% Parametres du reseau
N = 204;
d_max = 1500;            % km
Delta_t = 20;            % s

%% ============================================================
%  1. Recuperation des facteurs topologiques orbitaux
%% ============================================================

script_dir = fileparts(mfilename('fullpath'));
bridge_results_file = fullfile(script_dir, '..', 'chi_bridge_uniformite_orbitale_results.mat');

if ~isfile(bridge_results_file)
    error(['Fichier introuvable : %s\n' ...
           'Executer d''abord chi_bridge_uniformite_orbitale.'], ...
           bridge_results_file);
end

bridge_data = load(bridge_results_file);

if ~isfield(bridge_data,'results')
    error('Le fichier %s ne contient pas la structure results.', ...
        bridge_results_file);
end

bridge_results = bridge_data.results;

required_fields = { ...
    'chi_bridge', ...
    'mean_links_per_component', ...
    'mean_bridges_per_component', ...
    'beta0_used', ...
    'E_edges', ...
    'phi', ...
    'lambda_orb'};

for k = 1:numel(required_fields)
    if ~isfield(bridge_results,required_fields{k})
        error('Champ manquant dans results : %s', required_fields{k});
    end
end

chi_bridge_delta = double(bridge_results.chi_bridge);
mean_links_per_component = ...
    double(bridge_results.mean_links_per_component);
mean_bridges_per_component = ...
    double(bridge_results.mean_bridges_per_component);
beta0_delta = double(bridge_results.beta0_used);
E_edges_delta = double(bridge_results.E_edges);
phi_orb = double(bridge_results.phi(:));
lambda_orb = double(bridge_results.lambda_orb(:));

factor_source = ...
    'chi_bridge et nombre moyen de ponts charges depuis le fichier orbital';

% Vérifications de cohérence
if isfield(bridge_results,'N') && double(bridge_results.N) ~= N
    warning(['N differe entre ce code (%d) et le fichier chi_bridge ' ...
             '(%.0f).'], N, double(bridge_results.N));
end

if isfield(bridge_results,'R') && ...
        abs(double(bridge_results.R)-R) > 1e-9
    warning(['R differe entre ce code (%.6f km) et le fichier ' ...
             'chi_bridge (%.6f km).'], ...
             R,double(bridge_results.R));
end

if isfield(bridge_results,'dmax') && ...
        abs(double(bridge_results.dmax)-d_max) > 1e-9
    warning(['d_max differe entre ce code (%.6f km) et le fichier ' ...
             'chi_bridge (%.6f km).'], ...
             d_max,double(bridge_results.dmax));
end

validateattributes(chi_bridge_delta,{'numeric'}, ...
    {'scalar','real','finite','>=',0,'<=',1});

validateattributes(mean_bridges_per_component,{'numeric'}, ...
    {'scalar','real','finite','nonnegative'});

%% ============================================================
%  2. Probabilite de rupture d'un lien quelconque
%% ============================================================

q_break_link_raw = 2*v_rad*Delta_t/d_max;
q_break_link = min(max(q_break_link_raw,0),1);

if abs(q_break_link-q_break_link_raw) > 1e-12
    warning(['q_break lien brut = %.6f hors de [0,1]. ' ...
             'Valeur tronquee a %.6f.'], ...
             q_break_link_raw,q_break_link);
end

%% ============================================================
%  3. Conditionnement par le fait que le lien soit un pont
%% ============================================================

% Aire d'intersection plane de deux disques de rayon d_max
% dont les centres sont separes par d_max.
A_inter_at_dmax = ...
    2*d_max^2*acos(1/2) ...
    - 0.5*d_max*sqrt(3*d_max^2);

% Element de surface integre sur la longitude.
surface_weight = 2*pi*R^2.*cos(phi_orb);

% Les liens sont vus avec un poids proportionnel a lambda_orb^2.
link_weight = lambda_orb.^2 .* surface_weight;

% Probabilite locale qu'un lien de longueur d_max n'ait aucun
% voisin commun.
p_bridge_at_boundary_local = ...
    exp(-lambda_orb .* A_inter_at_dmax);

% Moyenne conditionnee par l'existence d'un lien.
p_bridge_at_boundary = ...
    trapz(phi_orb, p_bridge_at_boundary_local .* link_weight) ...
    / trapz(phi_orb, link_weight);

p_bridge_at_boundary = min(max(p_bridge_at_boundary,0),1);

% Facteur de Bayes :
%
% f_{D|pont}(d_max) / f_{D|lien}(d_max)
%   = P(pont | D=d_max,lien) / P(pont | lien).
bridge_boundary_factor = ...
    p_bridge_at_boundary / chi_bridge_delta;

% Probabilite conditionnelle qu'un pont se rompe.
q_break_bridge_raw = ...
    q_break_link * bridge_boundary_factor;

q_break_bridge = min(max(q_break_bridge_raw,0),1);

if abs(q_break_bridge-q_break_bridge_raw) > 1e-12
    warning(['q_break pont brut = %.6f hors de [0,1]. ' ...
             'Valeur tronquee a %.6f.'], ...
             q_break_bridge_raw,q_break_bridge);
end

%% ============================================================
%  4. Probabilite theorique de rupture d'une composante
%% ============================================================

bridge_factor_delta = mean_bridges_per_component;

% Forme probabiliste.
p_break_delta = ...
    1 - (1-q_break_bridge)^bridge_factor_delta;

p_break_delta = min(max(p_break_delta,0),1);

% Approximation lineaire conditionnelle :
% E[B_C] q_break^pont.
p_break_delta_linear = ...
    min(max(bridge_factor_delta*q_break_bridge,0),1);

% Forme lineaire equivalente apres simplification de chi_bridge :
% E[L_C] q_break^lien P(pont | D=d_max,lien).
p_break_delta_linear_simplified = ...
    min(max(mean_links_per_component * ...
            q_break_link * p_bridge_at_boundary,0),1);

%% ============================================================
%  5. Evolution en fonction du nombre moyen de ponts
%% ============================================================

bridge_factor_values = linspace(0, ...
    max(5,2*bridge_factor_delta),500);

p_break_vs_bridge_factor = ...
    1-(1-q_break_bridge).^bridge_factor_values;

figure;
plot(bridge_factor_values,p_break_vs_bridge_factor,'LineWidth',2);
hold on;
xline(bridge_factor_delta,':', ...
    sprintf('\\overline{B}_{\\mathcal C}=%.3f',bridge_factor_delta), ...
    'LineWidth',1.5);
yline(p_break_delta,':', ...
    sprintf('p_{break}^{\\Delta}=%.4f',p_break_delta), ...
    'LineWidth',1.5);
grid on;
xlabel('Nombre moyen de ponts par composante');
ylabel('p_{break}^{\Delta}');
title('Probabilite theorique de rupture corrigee');
legend('1-(1-q_{break}^{pont})^{\overline{B}_{\mathcal C}}', ...
       'Facteur utilise', ...
       'p_{break}^{\Delta}', ...
       'Location','best');
hold off;

%% ============================================================
%  6. Evolution en fonction du pas temporel
%% ============================================================

Delta_t_values = linspace(0,120,500);

q_break_link_vs_dt = ...
    2*v_rad.*Delta_t_values/d_max;

q_break_link_vs_dt = min(max(q_break_link_vs_dt,0),1);

q_break_bridge_vs_dt = ...
    q_break_link_vs_dt .* bridge_boundary_factor;

q_break_bridge_vs_dt = min(max(q_break_bridge_vs_dt,0),1);

p_break_vs_dt = ...
    1-(1-q_break_bridge_vs_dt).^bridge_factor_delta;

figure;
plot(Delta_t_values,p_break_vs_dt,'LineWidth',2);
xline(Delta_t,'--', ...
    sprintf('\\Delta t = %.1f s',Delta_t), ...
    'LineWidth',1.5);
yline(p_break_delta,':', ...
    sprintf('p_{break}^{Delta}=%.4f',p_break_delta), ...
    'LineWidth',1.5);
grid on;
xlabel('Pas temporel \Delta t (s)');
ylabel('p_{break}^{\Delta}');
title('Evolution de p_{break}^{\Delta} avec le pas temporel');
ylim([0,1]);

%% ============================================================
%  7. Affichage console
%% ============================================================

fprintf('\n=== p_break theorique Walker Delta ===\n');
fprintf('N                                  : %d\n',N);
fprintf('Rayon orbital                      : %.3f km\n',R);
fprintf('Vitesse orbitale                   : %.6f km/s\n',v_orb);
fprintf('Vitesse relative approximee        : %.6f km/s\n',v_rel);
fprintf('Vitesse radiale approximee         : %.6f km/s\n',v_rad);
fprintf('d_max                              : %.3f km\n',d_max);
fprintf('Delta_t                            : %.3f s\n',Delta_t);
fprintf('beta0_Delta utilise                : %.6f\n',beta0_delta);
fprintf('E[|E_Delta|]                       : %.6f\n',E_edges_delta);
fprintf('chi_bridge_Delta                   : %.10f\n',chi_bridge_delta);
fprintf('P(pont | D=d_max,lien)             : %.10f\n',p_bridge_at_boundary);
fprintf('Facteur conditionnel frontière     : %.10f\n',bridge_boundary_factor);
fprintf('Liens moyens par composante        : %.10f\n',mean_links_per_component);
fprintf('Ponts moyens par composante        : %.10f\n',mean_bridges_per_component);
fprintf('Source des facteurs                : %s\n',factor_source);
fprintf('q_break par lien quelconque        : %.10f\n',q_break_link);
fprintf('q_break conditionnel par pont      : %.10f\n',q_break_bridge);
fprintf('p_break_Delta probabiliste         : %.10f\n',p_break_delta);
fprintf('p_break_Delta lineaire             : %.10f\n',p_break_delta_linear);
fprintf('p_break lineaire simplifie         : %.10f\n', ...
    p_break_delta_linear_simplified);

%% ============================================================
%  8. Sauvegarde
%% ============================================================

save('pbreak_theorique_walker_delta_results.mat', ...
    'R_earth','h','R','mu_earth', ...
    'v_orb','v_rel','v_rad', ...
    'N','d_max','Delta_t', ...
    'beta0_delta','E_edges_delta', ...
    'chi_bridge_delta', ...
    'p_bridge_at_boundary','bridge_boundary_factor', ...
    'A_inter_at_dmax', ...
    'mean_links_per_component', ...
    'mean_bridges_per_component','bridge_factor_delta', ...
    'factor_source', ...
    'q_break_link_raw','q_break_link', ...
    'q_break_bridge_raw','q_break_bridge', ...
    'p_break_delta','p_break_delta_linear', ...
    'p_break_delta_linear_simplified', ...
    'bridge_factor_values','p_break_vs_bridge_factor', ...
    'Delta_t_values','q_break_link_vs_dt', ...
    'q_break_bridge_vs_dt','p_break_vs_dt');

fprintf('\nResultats sauvegardes dans :\n');
fprintf('pbreak_theorique_walker_delta_results.mat\n');
