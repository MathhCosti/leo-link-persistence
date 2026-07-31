%% compare_vrel_vrad.m
% Compare temporellement :
%
%   v_rel(t)  : norme moyenne de la vitesse relative entre les extrémités
%               des liens existants à l'instant t ;
%
%   v_rad(t)  : composante radiale sortante moyenne de cette vitesse
%               relative, c'est-à-dire la composante qui augmente la
%               distance entre les deux satellites :
%
%       v_rad,out = max( dot(v_j-v_i, (r_j-r_i)/||r_j-r_i||), 0 )
%
% Le script distingue aussi les liens au bord qui sont des ponts
% topologiques à l'instant t, afin de mesurer directement
% la corrélation entre vitesse radiale et caractère critique.
%
% Entrée :
%   analysis_temp_results.mat
%
% Sortie :
%   vrel_vrad_results.mat

clear; clc; close all;

%% ============================================================
%  1. Chargement
%% ============================================================
script_dir = fileparts(mfilename('fullpath'));
input_file = fullfile(script_dir, 'analysis_temp_results.mat');

if isempty(input_file)
    error(['Fichier introuvable. Place analysis_temp_results.mat ', ...
           'dans le dossier du script ou dans son dossier parent.']);
end

S = load(input_file);

required_fields = {'Positions', 'Adjacency', 'dt', 'dmax'};

for k_field = 1:numel(required_fields)
    if ~isfield(S, required_fields{k_field})
        error('Le fichier doit contenir la variable %s.', ...
              required_fields{k_field});
    end
end

Positions = S.Positions;
Adjacency = S.Adjacency;
dt = S.dt;
dmax = S.dmax;

if ~iscell(Positions) || ~iscell(Adjacency)
    error('Positions et Adjacency doivent être des cellules temporelles.');
end

Nt = min(numel(Positions), numel(Adjacency));

if Nt < 2
    error('Il faut au moins deux instants temporels.');
end

%% ============================================================
%  2. Vitesses instantanées des satellites
%
%  Différence centrée aux instants intérieurs ;
%  différence avant/arrière aux extrémités.
%% ============================================================

Velocities = cell(Nt,1);

for k = 1:Nt
    Pk = Positions{k};

    if size(Pk,2) ~= 3
        error('Positions{%d} doit être une matrice N x 3.', k);
    end

    if k == 1
        Velocities{k} = (Positions{2} - Positions{1}) / dt;

    elseif k == Nt
        Velocities{k} = (Positions{Nt} - Positions{Nt-1}) / dt;

    else
        Velocities{k} = (Positions{k+1} - Positions{k-1}) / (2*dt);
    end
end

%% ============================================================
%  2.b Vitesse orbitale empirique et modèle théorique de v_rel
%
%  Pour deux directions de vitesse indépendantes et uniformes,
%  séparées par un angle theta uniforme sur [0,pi] :
%
%    v_rel(theta) = 2*v_orb*sin(theta/2)
%
%  donc
%
%    E[v_rel] = (1/pi) int_0^pi 2*v_orb*sin(theta/2) dtheta
%             = (4/pi)*v_orb.
%% ============================================================

v_orb_mean_t = nan(Nt,1);

for k = 1:Nt
    speed_sat = vecnorm(Velocities{k},2,2);
    v_orb_mean_t(k) = mean(speed_sat,'omitnan');
end

v_orb_mean_global = mean(v_orb_mean_t,'omitnan');

% Modèle théorique évalué avec la vitesse orbitale empirique.
vrel_model_t = (4/pi) * v_orb_mean_t;
vrel_model_global = (4/pi) * v_orb_mean_global;

%% ============================================================
%  3. Calcul sur les liens existants
%% ============================================================

vrel_mean_t      = nan(Nt,1);
vrad_signed_t    = nan(Nt,1);
vrad_abs_t       = nan(Nt,1);
vrad_out_mean_t       = nan(Nt,1);
vrad_out_border_t     = nan(Nt,1);
vrel_border_t              = nan(Nt,1);
vrad_out_bridge_border_t   = nan(Nt,1);
vrel_bridge_border_t       = nan(Nt,1);
n_links_t                  = zeros(Nt,1);
n_border_links_t           = zeros(Nt,1);
n_bridge_links_t           = zeros(Nt,1);
n_bridge_border_links_t    = zeros(Nt,1);

% Stockage de toutes les observations, pour des moyennes globales
vrel_all      = [];
vrad_signed_all = [];
vrad_abs_all  = [];
vrad_out_all         = [];
vrad_out_border_all  = [];
vrel_border_all             = [];
vrad_out_bridge_border_all  = [];
vrel_bridge_border_all      = [];

for k = 1:Nt

    P = Positions{k};
    V = Velocities{k};

    A = logical(spones(Adjacency{k}));
    A = A | A.';
    A(1:size(A,1)+1:end) = false;

    [I,J] = find(triu(A,1));
    n_links_t(k) = numel(I);

    % Détection des ponts topologiques à l'instant t.
    % is_bridge_edge(q) correspond à l'arête (I(q),J(q)).
    is_bridge_edge = find_bridges_tarjan(A, I, J);
    n_bridge_links_t(k) = nnz(is_bridge_edge);

    if isempty(I)
        continue;
    end

    dr = P(J,:) - P(I,:);
    distance = vecnorm(dr,2,2);

    valid = distance > 0;
    I = I(valid);
    J = J(valid);
    dr = dr(valid,:);
    distance = distance(valid);
    is_bridge_edge = is_bridge_edge(valid);

    if isempty(I)
        continue;
    end

    e_rad = dr ./ distance;

    dv = V(J,:) - V(I,:);

    % Norme de la vitesse relative
    vrel = vecnorm(dv,2,2);

    % Composante radiale signée :
    % > 0 : éloignement ; < 0 : rapprochement
    vrad_signed = sum(dv .* e_rad, 2);

    % Valeur absolue de la composante radiale
    vrad_abs = abs(vrad_signed);

    % Composante radiale sortante, pertinente pour p_break
    vrad_out = max(vrad_signed, 0);

    % --------------------------------------------------------
    % Liens situés dans la couche de rupture au bord.
    %
    % Un lien est considéré "au bord" si, avec sa vitesse
    % radiale sortante actuelle, il peut atteindre dmax pendant dt :
    %
    %   distance + v_rad,out*dt >= dmax
    %
    % Cette définition sélectionne directement la population
    % cinématiquement susceptible de se rompre au prochain pas.
    % --------------------------------------------------------
    is_border = ...
        (vrad_out > 0) & ...
        (distance + vrad_out*dt >= dmax);

    n_border_links_t(k) = nnz(is_border);

    % Ponts appartenant à la couche de rupture
    is_bridge_border = is_border & is_bridge_edge;
    n_bridge_border_links_t(k) = nnz(is_bridge_border);

    if any(is_border)
        vrad_out_border_t(k) = mean(vrad_out(is_border));
        vrel_border_t(k) = mean(vrel(is_border));

        vrad_out_border_all = ...
            [vrad_out_border_all; vrad_out(is_border)]; %#ok<AGROW>
        vrel_border_all = ...
            [vrel_border_all; vrel(is_border)]; %#ok<AGROW>
    end

    if any(is_bridge_border)
        vrad_out_bridge_border_t(k) = ...
            mean(vrad_out(is_bridge_border));
        vrel_bridge_border_t(k) = ...
            mean(vrel(is_bridge_border));

        vrad_out_bridge_border_all = ...
            [vrad_out_bridge_border_all; ...
             vrad_out(is_bridge_border)]; %#ok<AGROW>

        vrel_bridge_border_all = ...
            [vrel_bridge_border_all; ...
             vrel(is_bridge_border)]; %#ok<AGROW>
    end

    vrel_mean_t(k)     = mean(vrel);
    vrad_signed_t(k)   = mean(vrad_signed);
    vrad_abs_t(k)      = mean(vrad_abs);
    vrad_out_mean_t(k) = mean(vrad_out);

    vrel_all = [vrel_all; vrel]; %#ok<AGROW>
    vrad_signed_all = [vrad_signed_all; vrad_signed]; %#ok<AGROW>
    vrad_abs_all = [vrad_abs_all; vrad_abs]; %#ok<AGROW>
    vrad_out_all = [vrad_out_all; vrad_out]; %#ok<AGROW>
end

%% ============================================================
%  4. Moyennes globales
%
%  Moyenne agrégée : chaque lien observé a le même poids.
%  Moyenne temporelle : chaque instant a le même poids.
%% ============================================================

vrel_mean_global = mean(vrel_all, 'omitnan');
vrad_signed_mean_global = mean(vrad_signed_all, 'omitnan');
vrad_abs_mean_global = mean(vrad_abs_all, 'omitnan');
vrad_out_mean_global = mean(vrad_out_all, 'omitnan');

vrad_out_border_global = mean(vrad_out_border_all, 'omitnan');
vrel_border_global = mean(vrel_border_all, 'omitnan');

vrad_out_bridge_border_global = ...
    mean(vrad_out_bridge_border_all, 'omitnan');
vrel_bridge_border_global = ...
    mean(vrel_bridge_border_all, 'omitnan');

vrel_mean_time = mean(vrel_mean_t, 'omitnan');
vrad_signed_mean_time = mean(vrad_signed_t, 'omitnan');
vrad_abs_mean_time = mean(vrad_abs_t, 'omitnan');
vrad_out_mean_time = mean(vrad_out_mean_t, 'omitnan');
vrad_out_border_time = mean(vrad_out_border_t, 'omitnan');
vrel_border_time = mean(vrel_border_t, 'omitnan');

vrad_out_bridge_border_time = ...
    mean(vrad_out_bridge_border_t, 'omitnan');
vrel_bridge_border_time = ...
    mean(vrel_bridge_border_t, 'omitnan');

%% ============================================================
%  5. Axe temporel
%% ============================================================

if isfield(S,'time_values') && numel(S.time_values) >= Nt
    time_values = S.time_values(1:Nt).';
else
    time_values = (0:Nt-1).' * dt;
end

%% ============================================================
%  6. Affichage
%% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' COMPARAISON v_rel ET v_rad SUR LES LIENS EXISTANTS\n');
fprintf('============================================================\n');
fprintf('Fichier chargé                         : %s\n', input_file);
fprintf('Nombre d''instants                     : %d\n', Nt);
fprintf('Nombre total de liens observés         : %d\n', numel(vrel_all));
fprintf('------------------------------------------------------------\n');
fprintf('Moyennes agrégées par lien\n');
fprintf('  <v_rel>                              : %.8f km/s\n', ...
        vrel_mean_global);
fprintf('  <v_rad signé>                        : %.8f km/s\n', ...
        vrad_signed_mean_global);
fprintf('  <|v_rad|>                            : %.8f km/s\n', ...
        vrad_abs_mean_global);
fprintf('  <v_rad sortant>                      : %.8f km/s\n', ...
        vrad_out_mean_global);
fprintf('  <v_rad sortant | lien au bord>       : %.8f km/s\n', ...
        vrad_out_border_global);
fprintf('  <v_rel | lien au bord>               : %.8f km/s\n', ...
        vrel_border_global);
fprintf('  Nombre de liens de bord observés     : %d\n', ...
        numel(vrad_out_border_all));
fprintf('  <v_rad sortant | pont au bord>       : %.8f km/s\n', ...
        vrad_out_bridge_border_global);
fprintf('  <v_rel | pont au bord>               : %.8f km/s\n', ...
        vrel_bridge_border_global);
fprintf('  Nombre de ponts au bord observés     : %d\n', ...
        numel(vrad_out_bridge_border_all));
fprintf('  Rapport vitesse pont/bord            : %.8f\n', ...
        vrad_out_bridge_border_global / vrad_out_border_global);
fprintf('------------------------------------------------------------\n');
fprintf('Moyennes temporelles des courbes\n');
fprintf('  mean_t(v_rel)                        : %.8f km/s\n', ...
        vrel_mean_time);
fprintf('  mean_t(v_rad signé)                  : %.8f km/s\n', ...
        vrad_signed_mean_time);
fprintf('  mean_t(|v_rad|)                      : %.8f km/s\n', ...
        vrad_abs_mean_time);
fprintf('  mean_t(v_rad sortant)                : %.8f km/s\n', ...
        vrad_out_mean_time);
fprintf('  mean_t(v_rad sortant | bord)         : %.8f km/s\n', ...
        vrad_out_border_time);
fprintf('  mean_t(v_rel | bord)                 : %.8f km/s\n', ...
        vrel_border_time);
fprintf('  mean_t(v_rad sortant | pont au bord) : %.8f km/s\n', ...
        vrad_out_bridge_border_time);
fprintf('  mean_t(v_rel | pont au bord)         : %.8f km/s\n', ...
        vrel_bridge_border_time);
fprintf('============================================================\n\n');

%% ============================================================
%  7. Courbe principale : v_rel et v_rad sortant
%% ============================================================

figure;
plot(time_values, vrel_mean_t, 'LineWidth', 1.5);
hold on;
plot(time_values, vrad_out_mean_t, 'LineWidth', 1.5);
plot(time_values, vrad_out_border_t, 'LineWidth', 1.5);
plot(time_values, vrad_out_bridge_border_t, 'LineWidth', 1.5);

yline(vrel_mean_global, '--', ...
    sprintf('<v_{rel}> = %.3f km/s', vrel_mean_global), ...
    'LabelHorizontalAlignment', 'left');

yline(vrad_out_mean_global, ':', ...
    sprintf('<v_{rad,out}> = %.3f km/s', vrad_out_mean_global), ...
    'LabelHorizontalAlignment', 'left');

yline(vrad_out_border_global, '-.', ...
    sprintf('<v_{rad,out}|bord> = %.3f km/s', ...
    vrad_out_border_global), ...
    'LabelHorizontalAlignment', 'left');

yline(vrad_out_bridge_border_global, '--', ...
    sprintf('<v_{rad,out}|pont,bord> = %.3f km/s', ...
    vrad_out_bridge_border_global), ...
    'LabelHorizontalAlignment', 'right');

grid on;
xlabel('Temps (s)');
ylabel('Vitesse moyenne (km/s)');
title('Comparaison de v_{rel} et de la vitesse radiale sortante');
legend('v_{rel}(t)', 'v_{rad,out}(t)', ...
       'v_{rad,out}(t) des liens au bord', ...
       'v_{rad,out}(t) des ponts au bord', ...
       'Moyenne v_{rel}', 'Moyenne v_{rad,out}', ...
       'Moyenne v_{rad,out} au bord', ...
       'Moyenne v_{rad,out} des ponts au bord', ...
       'Location', 'best');

%% ============================================================
%  8. Courbe complémentaire : toutes les définitions radiales
%% ============================================================

figure;
plot(time_values, vrel_mean_t, 'LineWidth', 1.4);
hold on;
plot(time_values, vrad_signed_t, 'LineWidth', 1.2);
plot(time_values, vrad_abs_t, 'LineWidth', 1.2);
plot(time_values, vrad_out_mean_t, 'LineWidth', 1.2);
grid on;
xlabel('Temps (s)');
ylabel('Vitesse moyenne (km/s)');
title('v_{rel} et différentes définitions de v_{rad}');
legend('v_{rel}', ...
       'v_{rad} signé', ...
       '|v_{rad}|', ...
       'v_{rad,out}', ...
       'Location', 'best');

%% ============================================================
%  9. Comparaison spécifique sur les liens au bord
%% ============================================================

figure;
plot(time_values, vrel_border_t, 'LineWidth', 1.4);
hold on;
plot(time_values, vrad_out_border_t, 'LineWidth', 1.4);
plot(time_values, vrad_out_bridge_border_t, 'LineWidth', 1.4);

yline(vrel_border_global, '--', ...
    sprintf('<v_{rel}|bord> = %.3f km/s', vrel_border_global), ...
    'LabelHorizontalAlignment', 'left');

yline(vrad_out_border_global, ':', ...
    sprintf('<v_{rad,out}|bord> = %.3f km/s', ...
    vrad_out_border_global), ...
    'LabelHorizontalAlignment', 'left');

yline(vrad_out_bridge_border_global, '-.', ...
    sprintf('<v_{rad,out}|pont,bord> = %.3f km/s', ...
    vrad_out_bridge_border_global), ...
    'LabelHorizontalAlignment', 'right');

grid on;
xlabel('Temps (s)');
ylabel('Vitesse moyenne (km/s)');
title('Vitesses des liens situés dans la couche de rupture');
legend('v_{rel} des liens au bord', ...
       'v_{rad,out} des liens au bord', ...
       'v_{rad,out} des ponts au bord', ...
       'Moyenne v_{rel} au bord', ...
       'Moyenne v_{rad,out} au bord', ...
       'Moyenne v_{rad,out} des ponts au bord', ...
       'Location', 'best');

%% ============================================================
%  10. Vérification du modèle v_rel = 4/pi * v_orb
%% ============================================================

ratio_vrel_model_t = vrel_mean_t ./ vrel_model_t;
ratio_vrel_model_t(vrel_model_t <= 0) = NaN;

ratio_vrel_vorb_t = vrel_mean_t ./ v_orb_mean_t;
ratio_vrel_vorb_t(v_orb_mean_t <= 0) = NaN;

ratio_vrel_model_global = ...
    vrel_mean_global / vrel_model_global;

ratio_vrel_vorb_global = ...
    vrel_mean_global / v_orb_mean_global;

figure;
plot(time_values, vrel_mean_t, 'LineWidth', 1.5);
hold on;
plot(time_values, vrel_model_t, '--', 'LineWidth', 1.5);
grid on;
xlabel('Temps (s)');
ylabel('Vitesse moyenne (km/s)');
title('Vérification de \langle v_{rel}angle = 4/\pi \langle v_{orb}angle');
legend('\langle v_{rel}angle empirique', ...
       '(4/\pi)\langle v_{orb}angle empirique', ...
       'Location', 'best');
hold off;

figure;
plot(time_values, ratio_vrel_vorb_t, 'LineWidth', 1.4);
hold on;
yline(4/pi, '--', ...
    sprintf('4/\pi = %.4f', 4/pi), ...
    'LineWidth', 1.5);
yline(ratio_vrel_vorb_global, ':', ...
    sprintf('Rapport global = %.4f', ratio_vrel_vorb_global), ...
    'LineWidth', 1.5);
grid on;
xlabel('Temps (s)');
ylabel('\langle v_{rel}angle / \langle v_{orb}angle');
title('Rapport entre vitesse relative et vitesse orbitale');
legend('Rapport empirique', '4/\pi', 'Rapport global', ...
       'Location', 'best');
hold off;

fprintf('\n');
fprintf('============================================================\n');
fprintf(' VERIFICATION DU MODELE v_rel = 4/pi * v_orb\n');
fprintf('============================================================\n');
fprintf('<v_orb> empirique                    : %.8f km/s\n', ...
    v_orb_mean_global);
fprintf('<v_rel> empirique                    : %.8f km/s\n', ...
    vrel_mean_global);
fprintf('(4/pi)<v_orb>                       : %.8f km/s\n', ...
    vrel_model_global);
fprintf('<v_rel>/<v_orb>                     : %.8f\n', ...
    ratio_vrel_vorb_global);
fprintf('4/pi                                : %.8f\n', 4/pi);
fprintf('<v_rel>/[(4/pi)<v_orb>]             : %.8f\n', ...
    ratio_vrel_model_global);
fprintf('============================================================\n');

%% ============================================================
%  11. Rapport radial / relatif
%% ============================================================

ratio_vrad_vrel_t = vrad_out_mean_t ./ vrel_mean_t;
ratio_vrad_vrel_t(vrel_mean_t <= 0) = NaN;

ratio_vrad_vrel_global = ...
    vrad_out_mean_global / vrel_mean_global;

figure;
plot(time_values, ratio_vrad_vrel_t, 'LineWidth', 1.4);
hold on;
yline(ratio_vrad_vrel_global, '--', ...
    sprintf('Rapport global = %.3f', ratio_vrad_vrel_global), ...
    'LabelHorizontalAlignment', 'left');
grid on;
xlabel('Temps (s)');
ylabel('v_{rad,out}/v_{rel}');
title('Part radiale sortante de la vitesse relative');

%% ============================================================
%  12. Rapport entre ponts au bord et ensemble des liens au bord
%% ============================================================

ratio_bridge_border_t = ...
    vrad_out_bridge_border_t ./ vrad_out_border_t;

ratio_bridge_border_t(vrad_out_border_t <= 0) = NaN;

ratio_bridge_border_global = ...
    vrad_out_bridge_border_global / vrad_out_border_global;

figure;
plot(time_values, ratio_bridge_border_t, 'LineWidth', 1.4);
hold on;
yline(ratio_bridge_border_global, '--', ...
    sprintf('Rapport global = %.3f', ratio_bridge_border_global), ...
    'LabelHorizontalAlignment', 'left');
grid on;
xlabel('Temps (s)');
ylabel('<v_{rad,out}|pont,bord> / <v_{rad,out}|bord>');
title('Survitesse radiale des ponts situés au bord');

%% ============================================================
%  13. Sauvegarde
%% ============================================================

output_file = fullfile(script_dir, 'vrel_vrad_emp_results.mat');

save(output_file, ...
    'time_values', 'n_links_t', ...
    'v_orb_mean_t', ...
    'v_orb_mean_global', ...
    'vrel_model_t', ...
    'vrel_model_global', ...
    'ratio_vrel_model_t', ...
    'ratio_vrel_model_global', ...
    'ratio_vrel_vorb_t', ...
    'ratio_vrel_vorb_global', ...
    'vrel_mean_t', ...
    'vrad_signed_t', ...
    'vrad_abs_t', ...
    'vrad_out_mean_t', ...
    'vrad_out_border_t', ...
    'vrel_border_t', ...
    'n_border_links_t', ...
    'n_bridge_links_t', ...
    'n_bridge_border_links_t', ...
    'vrad_out_bridge_border_t', ...
    'vrel_bridge_border_t', ...
    'vrel_mean_global', ...
    'vrad_signed_mean_global', ...
    'vrad_abs_mean_global', ...
    'vrad_out_mean_global', ...
    'vrad_out_border_global', ...
    'vrel_border_global', ...
    'vrad_out_bridge_border_global', ...
    'vrel_bridge_border_global', ...
    'vrel_mean_time', ...
    'vrad_signed_mean_time', ...
    'vrad_abs_mean_time', ...
    'vrad_out_mean_time', ...
    'vrad_out_border_time', ...
    'vrel_border_time', ...
    'vrad_out_bridge_border_time', ...
    'vrel_bridge_border_time', ...
    'ratio_vrad_vrel_t', ...
    'ratio_vrad_vrel_global', ...
    'ratio_bridge_border_t', ...
    'ratio_bridge_border_global');

fprintf('Résultats sauvegardés dans %s\n', output_file);


%% ============================================================
%  FONCTION LOCALE : détection des ponts par l'algorithme de Tarjan
%% ============================================================

function is_bridge_edge = find_bridges_tarjan(A, I, J)

    n = size(A,1);
    m = numel(I);

    is_bridge_edge = false(m,1);

    if m == 0
        return;
    end

    % Liste d'adjacence contenant aussi l'indice de chaque arête
    neighbors = cell(n,1);
    edge_ids  = cell(n,1);

    for e = 1:m
        u = I(e);
        v = J(e);

        neighbors{u}(end+1) = v; %#ok<AGROW>
        edge_ids{u}(end+1)  = e; %#ok<AGROW>

        neighbors{v}(end+1) = u; %#ok<AGROW>
        edge_ids{v}(end+1)  = e; %#ok<AGROW>
    end

    visited = false(n,1);
    discovery = zeros(n,1);
    low = zeros(n,1);
    timer = 0;

    for root = 1:n
        if ~visited(root)
            dfs(root, 0);
        end
    end

    function dfs(u, parent_edge)

        visited(u) = true;
        timer = timer + 1;
        discovery(u) = timer;
        low(u) = timer;

        for q = 1:numel(neighbors{u})

            v = neighbors{u}(q);
            e = edge_ids{u}(q);

            if e == parent_edge
                continue;
            end

            if ~visited(v)
                dfs(v, e);

                low(u) = min(low(u), low(v));

                if low(v) > discovery(u)
                    is_bridge_edge(e) = true;
                end

            else
                low(u) = min(low(u), discovery(v));
            end
        end
    end
end
