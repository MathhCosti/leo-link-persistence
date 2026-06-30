clear; clc; close all;

%% ============================================================
%  ÉTUDE TOPOLOGIQUE TEMPORELLE D'UN RÉSEAU LEO
%  Sans animation
%
%  Sorties :
%  - beta0(t) : nombre de composantes connexes
%  - beta1_graph(t) : nombre de cycles du graphe non rempli
%  - beta0 et beta1 sur la suite zigzag par unions
%
%  Zigzag construit :
%  G1 -> G1 union G2 <- G2 -> G2 union G3 <- G3 ...
%% ============================================================

%% Paramètres physiques
R_earth = 6371;      % km
h = 550;             % km
R = R_earth + h;     % rayon orbital

mu = 398600;              % km^3/s^2
omega = sqrt(mu / R^3);   % vitesse angulaire orbitale rad/s

%% Paramètres du processus de Poisson
lambda = 4e-7;       % satellites / km^2
surface_sphere = 4*pi*R^2;

N = poissrnd(lambda * surface_sphere);

fprintf('Nombre de satellites générés : N = %d\n', N);

%% Génération uniforme des positions initiales sur la sphère
u = rand(N,1);
phi = 2*pi*rand(N,1);
theta = acos(1 - 2*u);

x = R * sin(theta) .* cos(phi);
y = R * sin(theta) .* sin(phi);
z = R * cos(theta);

positions0 = [x y z];

%% Sens de rotation défini par un plan séparateur passant par les pôles
rotation_sign = ones(N, 1);
rotation_sign(y >= 0) = 1;
rotation_sign(y < 0) = -1;

%% Paramètres des liens et du temps
dmax = 1500;     % km
dt = 60;         % pas temporel en secondes
Tmax = 12000;     % durée totale de simulation

time_values = 0:dt:Tmax;
Nt = length(time_values);

%% Stockage
Positions = cell(Nt,1);
Adjacency = cell(Nt,1);

num_edges = zeros(Nt,1);

% Nombre de liens théorique issu directement de la géométrie Walker Star
% L_W(t) = sum_{i<j} 1_{d_ij(t) <= dmax}
% calculé via l'angle central, sans utiliser pdist.
num_edges_walker_theory = zeros(Nt,1);

beta0 = zeros(Nt,1);
beta1_graph = zeros(Nt,1);
largest_component = zeros(Nt,1);

% Seuil angulaire équivalent à la distance euclidienne dmax
% On tient aussi compte de la limite de visibilité géométrique LOS.
d_LOS = 2*sqrt(R^2 - R_earth^2);
d_eff = min(dmax, d_LOS);

alpha_max = 2*asin(d_eff/(2*R));
cos_alpha_max = cos(alpha_max);

% Minimum théorique approximatif : configuration uniforme sur la sphère
% L_min_th ≈ C(N,2) p_link, avec p_link = d_eff^2/(4R^2).
% Dans le Walker Star, cette valeur sert d'approximation du creux
% lorsque les satellites sont les moins concentrés.
p_link_uniform = d_eff^2/(4*R^2);
L_min_theory_uniform = N*(N-1)/2 * p_link_uniform;

% Période théorique du nombre de liens.
% La période orbitale complète vaut T_orb = 2*pi/omega.
% Pour le nombre total de liens dans un Walker Star polaire, la configuration
% de connectivité se répète approximativement toutes les demi-orbites.
T_orb = 2*pi/omega;
T_links_theory = T_orb/2;     % = pi/omega

%% ============================================================
%  1. CONSTRUCTION DES GRAPHES TEMPORELS G(t)
%% ============================================================

for k = 1:Nt

    t = time_values(k);

    %% Mouvement orbital
    phi_t = phi;
    theta_t = theta + rotation_sign * omega * t;

    x_t = R * sin(theta_t) .* cos(phi_t);
    y_t = R * sin(theta_t) .* sin(phi_t);
    z_t = R * cos(theta_t);

    positions_t = [x_t y_t z_t];

    %% Graphe de lien
    D = squareform(pdist(positions_t));
    A = (D <= dmax) & (D > 0);
    A = sparse(A);

    %% Nombre de liens théorique Walker Star
    % Angle central gamma_ij(t) entre deux satellites :
    % cos(gamma_ij) = cos(theta_i(t)) cos(theta_j(t))
    %               + sin(theta_i(t)) sin(theta_j(t)) cos(phi_i - phi_j)
    % Le lien existe si gamma_ij <= alpha_max,
    % soit cos(gamma_ij) >= cos(alpha_max).
    cos_gamma = cos(theta_t)*cos(theta_t)' + ...
                (sin(theta_t)*sin(theta_t)') .* cos(phi - phi');
    A_theory = (cos_gamma >= cos_alpha_max) & ~eye(N);
    num_edges_walker_theory(k) = nnz(triu(A_theory,1));

    %% Stockage
    Positions{k} = positions_t;
    Adjacency{k} = A;

    %% Mesures topologiques sur le graphe
    G = graph(A);

    comp = conncomp(G);
    beta0(k) = max(comp);

    comp_sizes = accumarray(comp', 1);
    largest_component(k) = max(comp_sizes);

    E = nnz(triu(A,1));
    num_edges(k) = E;

    % Nombre cyclomatique du graphe :
    % beta1 = E - V + C
    beta1_graph(k) = E - N + beta0(k);
end

% Moyennes temporelles du nombre de liens
% - moyenne théorique Walker Star : calcul géométrique déterministe
% - moyenne empirique : moyenne des liens réellement mesurés dans la simulation
L_walker_mean_theory = mean(num_edges_walker_theory);
L_empirical_mean = mean(num_edges);

% Sinusoïde théorique du nombre de liens.
% On impose :
%   - la moyenne théorique Walker Star,
%   - le minimum théorique uniforme,
%   - la période théorique T_links_theory.
%
% Comme le minimum est supposé atteint à t = 0 :
% L_sin(0) = L_min_theory_uniform.
A_links_theory = L_walker_mean_theory - L_min_theory_uniform;
L_links_sinus_theory = L_walker_mean_theory ...
    - A_links_theory * cos(2*pi*time_values(:)/T_links_theory);

fprintf('Nombre moyen théorique de liens Walker Star : %.2f\n', L_walker_mean_theory);
fprintf('Nombre moyen empirique de liens par simulation : %.2f\n', L_empirical_mean);
fprintf('Minimum théorique uniforme approximatif : %.2f\n', L_min_theory_uniform);
fprintf('Période théorique des liens : %.2f s\n', T_links_theory);
fprintf('Amplitude de la sinusoïde théorique : %.2f\n', A_links_theory);

%% ============================================================
%  MODELE PERIODIQUE DU TAUX DE DISPARITION p_disp(t)
%% ============================================================
% On utilise la chaîne de modélisation :
%   L(t) périodique  ->  p_disp(t) périodique
%
% Modèle retenu : opposition de phase avec le nombre de liens.
%   p_disp(t) = p0_theory * ((2*L_mean - L_sin(t))/L_mean)^beta_pdisp
%
% p0_theory est estimé par un modèle cinétique de rencontre entre
% composantes :
%   p0 ~= 2 * r_eff * v_rel * lambda_comp
% avec :
%   r_eff = R * alpha_max
%   v_rel ~= 4/pi * v_orb
%   lambda_comp ~= mean(beta0)/(4*pi*R^2)

r_eff = R * alpha_max;          % rayon de connexion sur la sphère, en km
v_orb = sqrt(mu/R);             % vitesse orbitale, en km/s
v_rel_mean = (4/pi) * v_orb;    % vitesse relative moyenne approximative, en km/s

beta0_mean = mean(beta0);
lambda_comp = beta0_mean / (4*pi*R^2);

p0_theory = 2 * r_eff * v_rel_mean * lambda_comp;   % en s^-1

% Exposant de sensibilité au nombre de liens.
% beta_pdisp = 1 : modèle linéaire.
% beta_pdisp > 1 : amplification non linéaire.
beta_pdisp = 1;

% Ancien modèle, en phase avec le nombre de liens :
% pdisp_inphase(t) = p0 * (L(t)/L_moy)^beta
pdisp_theory_inphase = p0_theory * ...
    (L_links_sinus_theory / L_walker_mean_theory).^beta_pdisp;
pdisp_theory_inphase = max(pdisp_theory_inphase, 0);

% Nouveau modèle théorique en opposition de phase :
% on remplace L(t) par son symétrique par rapport à la moyenne,
% L_anti(t) = 2*L_moy - L(t).
% Ainsi, quand le nombre de liens est maximal, pdisp est minimal,
% et inversement.
L_links_sinus_antiphase = 2*L_walker_mean_theory - L_links_sinus_theory;
L_links_sinus_antiphase = max(L_links_sinus_antiphase, 0);

pdisp_theory_antiphase = p0_theory * ...
    (L_links_sinus_antiphase / L_walker_mean_theory).^beta_pdisp;
pdisp_theory_antiphase = max(pdisp_theory_antiphase, 0);

% Par défaut, la courbe théorique affichée est maintenant l'opposition de phase.
pdisp_theory = pdisp_theory_antiphase;

fprintf('beta0 moyen : %.2f\n', beta0_mean);
fprintf('p0 théorique moyen : %.4e s^-1\n', p0_theory);
fprintf('Durée caractéristique 1/p0 : %.2f s\n', 1/p0_theory);
fprintf('pdisp anti-phase min / moyen / max : %.4e / %.4e / %.4e s^-1\n', ...
    min(pdisp_theory), mean(pdisp_theory), max(pdisp_theory));

%% ============================================================
%  2. GRAPHES TEMPORELS CLASSIQUES
%% ============================================================

figure;
plot(time_values, beta0, 'LineWidth', 2);
grid on;
xlabel('Temps (s)');
ylabel('\beta_0');
title('\beta_0(t) : nombre de composantes connexes');

figure;
plot(time_values, beta1_graph, 'LineWidth', 2);
grid on;
xlabel('Temps (s)');
ylabel('\beta_1 graphe');
title('\beta_1(t) du graphe non rempli');

figure;
plot(time_values, largest_component / N, 'LineWidth', 2);
grid on;
xlabel('Temps (s)');
ylabel('|C_{max}| / N');
title('Fraction de satellites dans la plus grande composante');

figure;
plot(time_values, num_edges, 'LineWidth', 2); hold on;
plot(time_values, L_links_sinus_theory, 'm--', 'LineWidth', 2);
yline(L_walker_mean_theory, 'r--', ...
    sprintf('Moyenne th. Walker = %.1f', L_walker_mean_theory), ...
    'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
yline(L_empirical_mean, 'k:', ...
    sprintf('Moyenne empirique = %.1f', L_empirical_mean), ...
    'LineWidth', 1.8, 'LabelHorizontalAlignment', 'left');
yline(L_min_theory_uniform, 'g-.', ...
    sprintf('Minimum th. uniforme = %.1f', L_min_theory_uniform), ...
    'LineWidth', 1.8, 'LabelHorizontalAlignment', 'left');
grid on;
xlabel('Temps (s)');
ylabel('Nombre de liens');
title('Nombre de liens inter-satellites');
legend('Simulation', ...
       'Sinusoïde théorique', ...
       'Moyenne théorique Walker Star', ...
       'Moyenne empirique simulation', ...
       'Minimum théorique uniforme', ...
       'Location', 'best');

% Proxy empirique basé uniquement sur le nombre de liens simulé.
% Attention : ce n'est PAS encore le vrai taux empirique de disparition.
% Le vrai taux empirique sera calculé plus bas à partir des morts des barres H0.
pdisp_proxy_links = p0_theory * (num_edges(:) / L_empirical_mean).^beta_pdisp;
pdisp_proxy_links = max(pdisp_proxy_links, 0);

%% ============================================================
%  3. CONSTRUCTION DU ZIGZAG PAR UNIONS
%
%  G1 -> G1 U G2 <- G2 -> G2 U G3 <- G3 ...
%% ============================================================

Nz = 2*Nt - 1;

ZigzagAdjacency = cell(Nz,1);
ZigzagLabels = zeros(Nz,1);

idx = 1;

for k = 1:Nt

    % Graphe réel G_k
    ZigzagAdjacency{idx} = Adjacency{k};
    ZigzagLabels(idx) = k;
    idx = idx + 1;

    % Graphe union G_k U G_{k+1}
    if k < Nt
        ZigzagAdjacency{idx} = Adjacency{k} | Adjacency{k+1};
        ZigzagLabels(idx) = k + 0.5;
        idx = idx + 1;
    end
end

%% ============================================================
%  4. BETTI SUR LA SUITE ZIGZAG
%% ============================================================

beta0_zigzag = zeros(Nz,1);
beta1_zigzag_graph = zeros(Nz,1);
num_edges_zigzag = zeros(Nz,1);
largest_component_zigzag = zeros(Nz,1);

for k = 1:Nz

    A = ZigzagAdjacency{k};
    G = graph(A);

    comp = conncomp(G);
    beta0_zigzag(k) = max(comp);

    comp_sizes = accumarray(comp', 1);
    largest_component_zigzag(k) = max(comp_sizes);

    E = nnz(triu(A,1));
    num_edges_zigzag(k) = E;

    beta1_zigzag_graph(k) = E - N + beta0_zigzag(k);
end

%% ============================================================
%  5. GRAPHES SUR LA SUITE ZIGZAG
%% ============================================================

figure;
plot(ZigzagLabels, beta0_zigzag, '-o', 'LineWidth', 1.5);
grid on;
xlabel('Indice temporel / demi-indice');
ylabel('\beta_0');
title('\beta_0 sur le zigzag par unions');

figure;
plot(ZigzagLabels, beta1_zigzag_graph, '-o', 'LineWidth', 1.5);
grid on;
xlabel('Indice temporel / demi-indice');
ylabel('\beta_1 graphe');
title('\beta_1 du graphe sur le zigzag par unions');

figure;
plot(ZigzagLabels, largest_component_zigzag / N, '-o', 'LineWidth', 1.5);
grid on;
xlabel('Indice temporel / demi-indice');
ylabel('|C_{max}| / N');
title('Composante géante sur le zigzag par unions');

% Moyenne empirique du nombre de liens sur le zigzag par unions
L_empirical_mean_zigzag = mean(num_edges_zigzag);
fprintf('Nombre moyen empirique de liens sur le zigzag : %.2f\n', L_empirical_mean_zigzag);

figure;
plot(ZigzagLabels, num_edges_zigzag, '-o', 'LineWidth', 1.5); hold on;
yline(L_empirical_mean_zigzag, 'k:', ...
    sprintf('Moyenne empirique zigzag = %.1f', L_empirical_mean_zigzag), ...
    'LineWidth', 1.8, 'LabelHorizontalAlignment', 'left');
grid on;
xlabel('Indice temporel / demi-indice');
ylabel('Nombre de liens');
title('Nombre de liens sur le zigzag par unions');
legend('Zigzag par unions', 'Moyenne empirique zigzag', 'Location', 'best');

%% ============================================================
%  6. TAUX EMPIRIQUE REEL p_disp(t) DEPUIS LES MORTS H0
%% ============================================================
% Définition empirique utilisée :
% p_disp_emp(t_k) ~= (# barres H0 vivantes à t_k qui meurent avant t_{k+1})
%                   / (# barres H0 vivantes à t_k * dt)
%
% Contrairement à pdisp_proxy_links, cette courbe vient directement des
% naissances/morts des composantes H0 dans le barcode zigzag.

% Conversion des labels zigzag en temps physiques
ZigzagTime = zeros(Nz,1);
for k = 1:Nz
    lab = ZigzagLabels(k);

    if abs(lab - round(lab)) < 1e-12
        idx_time = round(lab);
        ZigzagTime(k) = time_values(idx_time);
    else
        idx_time = floor(lab);
        ZigzagTime(k) = 0.5 * (time_values(idx_time) + time_values(idx_time+1));
    end
end

% Espaces H0 du zigzag : une base est donnée par les composantes connexes
component_labels_zigzag = cell(Nz,1);
h0_dims = zeros(Nz,1);

for k = 1:Nz
    A = ZigzagAdjacency{k};
    G = graph(A);

    comp = conncomp(G);
    comp = comp(:);

    component_labels_zigzag{k} = comp;
    h0_dims(k) = max(comp);
end

% Applications H0 entre objets successifs du zigzag
maps = cell(Nz-1,1);

for k = 1:Nz-1
    if mod(k,2) == 1
        % G_i -> G_i U G_{i+1}
        maps{k}.type = 'f';
        maps{k}.mat = build_H0_map( ...
            component_labels_zigzag{k}, ...
            component_labels_zigzag{k+1}, ...
            h0_dims(k), ...
            h0_dims(k+1));
    else
        % G_i U G_{i+1} <- G_{i+1}
        maps{k}.type = 'g';
        maps{k}.mat = build_H0_map( ...
            component_labels_zigzag{k+1}, ...
            component_labels_zigzag{k}, ...
            h0_dims(k+1), ...
            h0_dims(k));
    end
end

% Barcode zigzag H0
intervals_H0 = zigzag_barcode_from_module_mod2(h0_dims, maps);

birth_index_H0 = intervals_H0(:,1);
death_index_H0 = intervals_H0(:,2);

birth_time_H0 = ZigzagTime(birth_index_H0);
death_time_H0 = ZigzagTime(death_index_H0);
lifetimes_H0 = death_time_H0 - birth_time_H0;

% Taux empirique réel sur les intervalles [t_k, t_{k+1}]
pdisp_emp_H0 = nan(Nt,1);
num_alive_H0 = zeros(Nt,1);
num_deaths_H0 = zeros(Nt,1);

for k = 1:Nt-1
    t0 = time_values(k);
    t1 = time_values(k+1);

    alive_k = (birth_time_H0 <= t0) & (death_time_H0 > t0);
    deaths_k = alive_k & (death_time_H0 <= t1);

    num_alive_H0(k) = sum(alive_k);
    num_deaths_H0(k) = sum(deaths_k);

    if num_alive_H0(k) > 0
        pdisp_emp_H0(k) = num_deaths_H0(k) / (num_alive_H0(k) * dt);
    end
end

pdisp_emp_H0(end) = pdisp_emp_H0(end-1);
num_alive_H0(end) = num_alive_H0(end-1);
num_deaths_H0(end) = num_deaths_H0(end-1);

% Lissage léger pour rendre le taux lisible
smoothing_window = 5;
pdisp_emp_H0_smooth = movmean(pdisp_emp_H0, smoothing_window, 'omitnan');

fprintf('Nombre total de barres H0 : %d\n', size(intervals_H0,1));
fprintf('pdisp empirique H0 min / moyen / max : %.4e / %.4e / %.4e s^-1\n', ...
    min(pdisp_emp_H0_smooth, [], 'omitnan'), ...
    mean(pdisp_emp_H0_smooth, 'omitnan'), ...
    max(pdisp_emp_H0_smooth, [], 'omitnan'));

% Graphe comparatif de p_disp(t)
figure;
plot(time_values, pdisp_theory, 'm--', 'LineWidth', 2); hold on;
plot(time_values, pdisp_proxy_links, 'b-', 'LineWidth', 1.2);
plot(time_values, pdisp_emp_H0_smooth, 'k-', 'LineWidth', 2);
yline(p0_theory, 'r--', ...
    sprintf('p_0 théorique moyen = %.2e s^{-1}', p0_theory), ...
    'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
grid on;
xlabel('Temps (s)');
ylabel('p_{disp}(t) (s^{-1})');
title('Taux de disparition p_{disp}(t)');
legend('p_{disp}^{th}(t) anti-phase', ...
       'Proxy depuis nombre de liens simulés', ...
       'p_{disp}^{emp}(t) depuis morts H_0', ...
       'p_0 moyen', ...
       'Location', 'best');


%% ============================================================
%  7. ETUDE SPECTRALE DE p_disp^{emp}(t)
%% ============================================================
% Objectif : vérifier si la courbe empirique noire peut être comprise
% comme une constante plus quelques harmoniques périodiques.
%
% On analyse la série pdisp_emp_H0_smooth(t), issue des morts H0.

p_spectral = pdisp_emp_H0_smooth(:);
t_spectral = time_values(:);

% Remplacement des NaN éventuels pour pouvoir appliquer la FFT
p_spectral = fillmissing(p_spectral, 'linear', 'EndValues', 'nearest');

% Pas de temps et fréquence d'échantillonnage
Delta_t_spectral = mean(diff(t_spectral));
Fs_spectral = 1/Delta_t_spectral;
M_spectral = length(p_spectral);

% Partie constante
pdisp_DC = mean(p_spectral, 'omitnan');

% Signal centré pour visualiser les harmoniques non nulles
x_spectral = p_spectral - pdisp_DC;

% Fenêtre de Hann manuelle pour limiter les fuites spectrales
n_spectral = (0:M_spectral-1)';
w_spectral = 0.5 * (1 - cos(2*pi*n_spectral/(M_spectral-1)));
coherent_gain = mean(w_spectral);

% FFT et spectre d'amplitude mono-côté
Y_spectral = fft(x_spectral .* w_spectral);
P2_spectral = abs(Y_spectral) / (M_spectral * coherent_gain);
P1_spectral = P2_spectral(1:floor(M_spectral/2)+1);

if length(P1_spectral) > 2
    P1_spectral(2:end-1) = 2*P1_spectral(2:end-1);
end

f_spectral = Fs_spectral*(0:floor(M_spectral/2))'/M_spectral;
period_spectral = 1 ./ f_spectral;
period_spectral(1) = Inf;

% Recherche des pics dominants hors fréquence nulle
P_for_peaks = P1_spectral;
P_for_peaks(1) = 0;

% On ignore les périodes plus longues que toute la simulation,
% sinon la tendance très basse fréquence peut dominer artificiellement.
f_min_useful = 1/(t_spectral(end)-t_spectral(1));
P_for_peaks(f_spectral < f_min_useful) = 0;

[~, peak_order] = sort(P_for_peaks, 'descend');
num_peaks_to_print = min(5, length(peak_order));
peak_indices = peak_order(1:num_peaks_to_print);

fprintf('\nEtude spectrale de pdisp empirique H0 :\n');
fprintf('Composante constante DC = %.4e s^-1\n', pdisp_DC);
for q = 1:num_peaks_to_print
    idxp = peak_indices(q);
    if f_spectral(idxp) > 0
        fprintf('Pic %d : f = %.4e Hz, T = %.1f s, amplitude = %.4e s^-1\n', ...
            q, f_spectral(idxp), 1/f_spectral(idxp), P1_spectral(idxp));
    end
end

% Fréquences théoriques attendues à partir de la période des liens
f1_links_theory = 1/T_links_theory;
f2_links_theory = 2/T_links_theory;

% Graphe spectral
figure;
plot(f_spectral, P1_spectral, 'k-', 'LineWidth', 1.5); hold on;
plot(f_spectral(peak_indices), P1_spectral(peak_indices), 'ro', ...
    'MarkerFaceColor', 'r', 'MarkerSize', 5);
xline(f1_links_theory, 'm--', ...
    sprintf('f_L = 1/T_L = %.2e Hz', f1_links_theory), ...
    'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
xline(f2_links_theory, 'b--', ...
    sprintf('2f_L = %.2e Hz', f2_links_theory), ...
    'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
grid on;
xlabel('Fréquence (Hz)');
ylabel('Amplitude spectrale de p_{disp}^{emp} (s^{-1})');
title('Analyse spectrale du taux empirique p_{disp}^{emp}(t)');
legend('Spectre FFT', 'Pics dominants', 'Fondamentale théorique', ...
       'Deuxième harmonique théorique', 'Location', 'best');
xlim([0, min(5*f1_links_theory, max(f_spectral))]);


% Graphe diagnostic : nombre de barres vivantes et mortes par intervalle
figure;
yyaxis left;
plot(time_values, num_alive_H0, 'LineWidth', 1.5);
ylabel('Nombre de barres H_0 vivantes');
yyaxis right;
stem(time_values, num_deaths_H0, 'filled');
ylabel('Morts H_0 par intervalle');
grid on;
xlabel('Temps (s)');
title('Diagnostic du taux empirique H_0');
legend('Barres vivantes', 'Morts par intervalle', 'Location', 'best');

%% ============================================================
%  6. SAUVEGARDE DES DONNÉES
%% ============================================================

save('leo_zigzag_analysis_results.mat', ...
    'N', 'R', 'h', 'lambda', 'dmax', 'dt', 'Tmax', ...
    'time_values', ...
    'Positions', 'Adjacency', ...
    'beta0', 'beta1_graph', 'largest_component', 'num_edges', ...
    'num_edges_walker_theory', 'L_walker_mean_theory', ...
    'L_empirical_mean', 'L_min_theory_uniform', ...
    'L_links_sinus_theory', 'A_links_theory', 'T_orb', 'T_links_theory', ...
    'p_link_uniform', 'd_LOS', 'd_eff', 'alpha_max', ...
    'r_eff', 'v_orb', 'v_rel_mean', 'beta0_mean', 'lambda_comp', ...
    'p0_theory', 'beta_pdisp', 'pdisp_theory', 'pdisp_theory_inphase', ...
    'pdisp_theory_antiphase', 'L_links_sinus_antiphase', 'pdisp_proxy_links', ...
    'ZigzagAdjacency', 'ZigzagLabels', 'ZigzagTime', ...
    'intervals_H0', 'birth_time_H0', 'death_time_H0', 'lifetimes_H0', ...
    'pdisp_emp_H0', 'pdisp_emp_H0_smooth', 'num_alive_H0', 'num_deaths_H0', ...
    'pdisp_DC', 'f_spectral', 'P1_spectral', 'period_spectral', ...
    'peak_indices', 'f1_links_theory', 'f2_links_theory', ...
    'beta0_zigzag', 'beta1_zigzag_graph', ...
    'largest_component_zigzag', 'num_edges_zigzag', ...
    'L_empirical_mean_zigzag');

fprintf('\nAnalyse terminée.\n');
fprintf('Résultats sauvegardés dans leo_zigzag_analysis_results.mat\n');
%% ============================================================
%  FONCTIONS LOCALES POUR LE BARCODE ZIGZAG H0
%% ============================================================

function M = build_H0_map(labels_source, labels_target, dim_source, dim_target)
    % Construit la matrice induite en H0 par une inclusion de graphes.
    % Chaque composante source est envoyée vers la composante cible
    % qui la contient. M est de taille dim_target x dim_source.

    M = zeros(dim_target, dim_source);

    for c = 1:dim_source
        vertices = find(labels_source == c);
        target_comps = unique(labels_target(vertices));

        if length(target_comps) ~= 1
            error(['Inclusion invalide pour H0 : une composante source ', ...
                   'est envoyée dans plusieurs composantes cibles.']);
        end

        target_c = target_comps(1);
        M(target_c, c) = 1;
    end

    M = mod(M,2);
end

function intervals = zigzag_barcode_from_module_mod2(dims, maps)
    % Calcule le barcode zigzag d'un module sur F2 à partir des dimensions
    % et des matrices d'applications.

    n = length(dims);

    R = cell(2,1);
    R{1} = zeros(dims(1),0);
    R{2} = eye(dims(1));

    b = 1;
    r = filtration_quotient_dims(R);
    intervals = [];

    for k = 1:n-1
        current_type = maps{k}.type;

        if current_type == 'f'
            M = maps{k}.mat;
            Rnext = cell(length(R)+1,1);

            for i = 1:length(R)
                Rnext{i} = gf2_col_basis(M * R{i});
            end

            Rnext{end} = eye(dims(k+1));
            bnext = [b, k+1];
            rnext = filtration_quotient_dims(Rnext);

            for i = 1:length(r)
                c = r(i) - rnext(i);
                if c > 0
                    intervals = [intervals; repmat([b(i), k], c, 1)]; %#ok<AGROW>
                end
            end

        elseif current_type == 'g'
            Nmat = maps{k}.mat;
            Rnext = cell(length(R)+1,1);
            Rnext{1} = zeros(dims(k+1),0);

            for i = 1:length(R)
                Rnext{i+1} = gf2_preimage(Nmat, R{i});
            end

            bnext = [k+1, b];
            rnext = filtration_quotient_dims(Rnext);

            for i = 1:length(r)
                c = r(i) - rnext(i+1);
                if c > 0
                    intervals = [intervals; repmat([b(i), k], c, 1)]; %#ok<AGROW>
                end
            end
        else
            error('Type de flèche inconnu.');
        end

        R = Rnext;
        b = bnext;
        r = rnext;
    end

    for i = 1:length(r)
        c = r(i);
        if c > 0
            intervals = [intervals; repmat([b(i), n], c, 1)]; %#ok<AGROW>
        end
    end
end

function dims = filtration_quotient_dims(R)
    m = length(R) - 1;
    dims = zeros(1,m);

    for i = 1:m
        dims(i) = gf2_rank(R{i+1}) - gf2_rank(R{i});
    end
end

function P = gf2_preimage(A, S)
    A = mod(full(A),2);
    S = mod(full(S),2);

    n = size(A,2);
    Big = [A S];
    Z = gf2_null(Big);
    X = Z(1:n,:);
    P = gf2_col_basis(X);
end

function B = gf2_col_basis(A)
    A = mod(full(A),2);

    if isempty(A)
        B = zeros(size(A,1),0);
        return;
    end

    [~, pivots] = gf2_rref(A);

    if isempty(pivots)
        B = zeros(size(A,1),0);
    else
        B = mod(A(:,pivots),2);
    end
end

function r = gf2_rank(A)
    A = mod(full(A),2);

    if isempty(A)
        r = 0;
        return;
    end

    [~, pivots] = gf2_rref(A);
    r = length(pivots);
end

function Z = gf2_null(A)
    A = mod(full(A),2);
    [R, pivots] = gf2_rref(A);

    n = size(A,2);
    free_cols = setdiff(1:n, pivots);

    if isempty(free_cols)
        Z = zeros(n,0);
        return;
    end

    Z = zeros(n, length(free_cols));

    for j = 1:length(free_cols)
        f = free_cols(j);
        z = zeros(n,1);
        z(f) = 1;

        for p = 1:length(pivots)
            col = pivots(p);
            z(col) = R(p,f);
        end

        Z(:,j) = mod(z,2);
    end
end

function [R, pivots] = gf2_rref(A)
    A = mod(full(A),2);
    [m,n] = size(A);

    R = A;
    pivots = [];
    row = 1;

    for col = 1:n
        if row > m
            break;
        end

        pivot_rel = find(R(row:m,col), 1);

        if isempty(pivot_rel)
            continue;
        end

        pivot = pivot_rel + row - 1;

        if pivot ~= row
            tmp = R(row,:);
            R(row,:) = R(pivot,:);
            R(pivot,:) = tmp;
        end

        for rr = 1:m
            if rr ~= row && R(rr,col) == 1
                R(rr,:) = mod(R(rr,:) + R(row,:), 2);
            end
        end

        pivots(end+1) = col; %#ok<AGROW>
        row = row + 1;
    end
end
