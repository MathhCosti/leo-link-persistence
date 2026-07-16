clear; clc; close all;

%% ============================================================
%  p_break THEORIQUE - WALKER DELTA
%
%  Etape 1 : probabilite de rupture d'un lien moyen
%
%      q_break^Delta
%        = 2 v_rad^Delta Delta_t / d_max
%
%  avec les approximations retenues :
%
%      v_rel^Delta = v_orb/sqrt(2)
%      v_rad^Delta = v_rel^Delta/pi
%
%  donc :
%
%      q_break^Delta
%        = (sqrt(2)/pi) v_orb Delta_t / d_max.
%
%  Etape 2 : correction topologique de la rupture d'un lien
%
%      chi^Delta = (N-beta0^Delta) / E[|E^Delta|]
%
%  La probabilite de rupture topologiquement efficace est :
%
%      p_break^Delta = q_break^Delta * chi^Delta.
%% ============================================================

%% Parametres physiques
R_earth = 6371;          % km
h = 550;                 % km
R = R_earth + h;         % km

mu_earth = 398600.4418;  % km^3/s^2

v_orb = sqrt(mu_earth/R);    % km/s
v_rel = v_orb/sqrt(2);       % km/s
v_rad = v_rel/pi;            % km/s

%% Parametres du reseau
N = 250;
d_max = 1500;            % km
Delta_t = 20;            % s

%% ============================================================
%  1. Recuperation de beta0 theorique
%% ============================================================

script_dir = fileparts(mfilename('fullpath'));
betti_results_file = fullfile(script_dir, '..', 'N1_N2_N3_walker_delta_results.mat');

use_betti_file = true;

if use_betti_file && isfile(betti_results_file)

    betti_data = load(betti_results_file);

    if isfield(betti_data,'beta0_th_123')
        beta0_delta = betti_data.beta0_th_123;
        beta0_source = 'beta0_th_123 charge depuis le fichier';
    elseif isfield(betti_data,'beta0_delta')
        beta0_delta = betti_data.beta0_delta;
        beta0_source = 'beta0_delta charge depuis le fichier';
    else
        error(['Le fichier %s ne contient ni beta0_th_123 ' ...
               'ni beta0_delta.'],betti_results_file);
    end

    if isfield(betti_data,'N') && betti_data.N ~= N
        warning(['N differe entre ce code (%d) et le fichier Betti (%d).'], ...
            N,betti_data.N);
    end

else
    % Valeur manuelle de secours.
    beta0_delta = 30;
    beta0_source = 'valeur manuelle';
end

%% ============================================================
%  1.b Recuperation du nombre moyen d'aretes
%% ============================================================

plink_results_file = fullfile(script_dir, '..', ...
    'plink_walker_delta_semi_analytique_results.mat');

if isfile(plink_results_file)
    plink_data = load(plink_results_file);

    if isfield(plink_data,'E_delta')
        E_edges_delta = plink_data.E_delta;
        edges_source = 'E_delta charge depuis le fichier';
    elseif isfield(plink_data,'p_link_delta')
        p_link_delta = plink_data.p_link_delta;
        E_edges_delta = nchoosek(N,2)*p_link_delta;
        edges_source = 'recalcule depuis p_link_delta';
    else
        error(['Le fichier %s ne contient ni E_delta ' ...
               'ni p_link_delta.'],plink_results_file);
    end
else
    error(['Fichier introuvable : %s\n' ...
           'Executer d''abord le code de quadrature de p_link.'], ...
           plink_results_file);
end

%% ============================================================
%  2. Probabilite theorique de rupture d'un lien
%% ============================================================

q_break_raw = 2*v_rad*Delta_t/d_max;

% Une probabilite doit rester dans [0,1].
q_break_delta = min(max(q_break_raw,0),1);

if abs(q_break_delta-q_break_raw) > 1e-12
    warning(['q_break brut = %.6f hors de [0,1]. ' ...
             'Valeur tronquee a %.6f.'], ...
             q_break_raw,q_break_delta);
end

%% ============================================================
%  3. Facteur correctif topologique
%% ============================================================

chi_delta_raw = (N-beta0_delta)/E_edges_delta;
chi_delta = min(max(chi_delta_raw,0),1);

if abs(chi_delta-chi_delta_raw) > 1e-12
    warning(['chi_Delta brut = %.6f hors de [0,1]. ' ...
             'Valeur tronquee a %.6f.'], ...
             chi_delta_raw,chi_delta);
end

%% ============================================================
%  4. Probabilite theorique de rupture topologiquement efficace
%% ============================================================

p_break_delta = q_break_delta*chi_delta;

%% ============================================================
%  5. Evolution en fonction du facteur correctif
%% ============================================================

chi_values = linspace(0,1,500);
p_break_vs_chi = q_break_delta.*chi_values;

figure;
plot(chi_values,p_break_vs_chi,'LineWidth',2); hold on;
xline(chi_delta,':',sprintf('\chi_{\Delta}=%.3f',chi_delta),'LineWidth',1.5);
yline(p_break_delta,':',sprintf('p_{break}^{\Delta}=%.4f',p_break_delta),'LineWidth',1.5);
grid on;
xlabel('Facteur correctif \chi_{\Delta}');
ylabel('p_{break}^{\Delta}');
title('Probabilite theorique de rupture corrigee');
legend('q_{break}^{\Delta}\chi_{\Delta}', ...
       '\chi_{\Delta} utilise', ...
       'p_{break}^{\Delta}', ...
       'Location','best');
hold off;

%% ============================================================
%  6. Evolution en fonction du pas temporel
%% ============================================================

Delta_t_values = linspace(0,120,500);

q_break_vs_dt = ...
    2*v_rad.*Delta_t_values/d_max;

q_break_vs_dt = min(max(q_break_vs_dt,0),1);

p_break_vs_dt = ...
    q_break_vs_dt.*chi_delta;

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
fprintf('beta0_Delta                        : %.6f\n',beta0_delta);
fprintf('Source de beta0                    : %s\n',beta0_source);
fprintf('E[|E_Delta|]                       : %.6f\n',E_edges_delta);
fprintf('Source des aretes                  : %s\n',edges_source);
fprintf('chi_Delta brut                     : %.10f\n',chi_delta_raw);
fprintf('chi_Delta utilise                  : %.10f\n',chi_delta);
fprintf('q_break_Delta par lien             : %.10f\n',q_break_delta);
fprintf('p_break_Delta = q_break*chi        : %.10f\n',p_break_delta);

%% ============================================================
%  8. Sauvegarde
%% ============================================================

save('pbreak_theorique_walker_delta_results.mat', ...
    'R_earth','h','R','mu_earth', ...
    'v_orb','v_rel','v_rad', ...
    'N','d_max','Delta_t', ...
    'beta0_delta','beta0_source', ...
    'E_edges_delta','edges_source', ...
    'chi_delta_raw','chi_delta', ...
    'q_break_raw','q_break_delta','p_break_delta', ...
    'chi_values','p_break_vs_chi', ...
    'Delta_t_values','p_break_vs_dt');

fprintf('\nResultats sauvegardes dans :\n');
fprintf('pbreak_theorique_walker_delta_results.mat\n');
