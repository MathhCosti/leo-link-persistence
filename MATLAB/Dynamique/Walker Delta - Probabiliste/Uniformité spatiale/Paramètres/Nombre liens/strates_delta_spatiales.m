function strates = strates_delta_spatiales(lambda, R, alpha_max, inc, varargin)
%STRATES_DELTA_SPATIALES
% Décomposition en strates autour des latitudes extrêmes +/-inc.
%
% beta mesure ici la distance angulaire à la frontière de la bande :
%
%   beta = inc - abs(latitude)
%
% beta = 0   : frontière orbitale |latitude| = inc
% beta = inc : équateur
%
% Les surfaces comprennent simultanément les hémisphères Nord et Sud.
%
% Entrées :
%   lambda    : densité surfacique initiale dans la bande [sat/km^2]
%   R         : rayon orbital [km]
%   alpha_max : angle central maximal de liaison [rad]
%   inc       : inclinaison commune [rad], avec 0 < inc <= pi/2
%
% Options :
%   'beta_step' : largeur angulaire des strates, défaut alpha_max/2
%   'beta_max'  : distance maximale à la frontière, défaut inc
%   'verbose'   : affiche la table récapitulative, défaut false
%
% Sortie :
%   strates.all_table
%   strates.keep_table
%   strates.active_table
%   strates.A_link
%   strates.beta_step
%   strates.beta_max
%   strates.beta_stop
%   strates.inc
%
% Remarque :
%   A_source repose sur une approximation de transport analogue au modèle
%   Walker Star. Pour une densité temporelle effective, il est préférable
%   de compter directement les satellites présents dans chaque strate.

    %% Vérification des entrées
    parser = inputParser;

    parser.addRequired('lambda', ...
        @(x) isnumeric(x) && isscalar(x) && x >= 0);

    parser.addRequired('R', ...
        @(x) isnumeric(x) && isscalar(x) && x > 0);

    parser.addRequired('alpha_max', ...
        @(x) isnumeric(x) && isscalar(x) && x > 0);

    parser.addRequired('inc', ...
        @(x) isnumeric(x) && isscalar(x) && x > 0 && x <= pi/2);

    parser.addParameter('beta_step', alpha_max/2, ...
        @(x) isnumeric(x) && isscalar(x) && x > 0);

    parser.addParameter('beta_max', inc, ...
        @(x) isnumeric(x) && isscalar(x) && x > 0);

    parser.addParameter('verbose', false, ...
        @(x) islogical(x) || isnumeric(x));

    parser.parse(lambda, R, alpha_max, inc, varargin{:});

    beta_step = parser.Results.beta_step;
    beta_max = min(parser.Results.beta_max, inc);
    verbose = logical(parser.Results.verbose);

    %% Frontières des strates
    beta_edges = 0:beta_step:beta_max;

    if isempty(beta_edges) || beta_edges(end) < beta_max - 1e-12
        beta_edges(end+1) = beta_max; %#ok<AGROW>
    else
        beta_edges(end) = beta_max;
    end

    beta_edges = min(beta_edges, inc);
    beta_edges = unique(beta_edges, 'stable');

    n_strates = length(beta_edges) - 1;

    beta_in = beta_edges(1:end-1)';
    beta_out = beta_edges(2:end)';
    beta_mid = 0.5*(beta_in + beta_out);

    %% Latitudes correspondant aux strates
    latitude_outer = inc - beta_in;
    latitude_inner = inc - beta_out;
    latitude_mid = inc - beta_mid;

    %% Aire de voisinage d'un lien
    A_link = 2*pi*R^2*(1-cos(alpha_max));

    %% Surface réelle des strates finales
    % Deux bandes symétriques autour de +inc et -inc.
    A_zone = 4*pi*R^2 .* ...
        (sin(latitude_outer)-sin(latitude_inner));

    %% Surface initiale alimentant chaque strate
    % Distribution initiale uniforme dans la bande.
    A_source = 4*pi*R^2 .* ...
        (sin(beta_out)-sin(beta_in));

    %% Nombre moyen de satellites transportés
    mu_strate = lambda .* A_source;

    %% Densité locale maximale estimée
    lambda_locale = mu_strate ./ A_zone;

    %% Degré moyen local
    k_mean_local = lambda_locale .* A_link;

    keep = true(n_strates,1);
    active_mask = true(n_strates,1);

    if n_strates > 0
        beta_stop = beta_out(end);
    else
        beta_stop = 0;
    end

    %% Table récapitulative
    index = (1:n_strates)';

    all_table = table( ...
        index, ...
        beta_in, beta_out, beta_mid, ...
        latitude_inner, latitude_outer, latitude_mid, ...
        rad2deg(beta_in), rad2deg(beta_out), ...
        rad2deg(latitude_inner), rad2deg(latitude_outer), ...
        A_source, A_zone, mu_strate, ...
        lambda_locale, k_mean_local, keep, active_mask, ...
        'VariableNames', { ...
        'index', ...
        'beta_in', 'beta_out', 'beta_mid', ...
        'latitude_inner', 'latitude_outer', 'latitude_mid', ...
        'beta_in_deg', 'beta_out_deg', ...
        'latitude_inner_deg', 'latitude_outer_deg', ...
        'A_source', 'A_zone', 'mu', ...
        'lambda_local', 'k_mean_local', 'keep', 'active'});

    %% Sortie
    strates = struct();

    strates.all_table = all_table;
    strates.keep_table = all_table;
    strates.active_table = all_table;

    strates.A_link = A_link;
    strates.beta_step = beta_step;
    strates.beta_max = beta_max;
    strates.beta_stop = beta_stop;
    strates.inc = inc;

    %% Affichage
    if verbose
        fprintf('\n=== Strates Walker Delta spatiales ===\n');
        fprintf('Inclinaison : %.2f deg\n', rad2deg(inc));
        fprintf('Aire totale de la bande : %.3e km^2\n', ...
            4*pi*R^2*sin(inc));
        fprintf('A_link : %.3e km^2\n', A_link);
        fprintf('beta_step : %.4f rad = %.2f deg\n\n', ...
            beta_step, rad2deg(beta_step));

        disp(all_table(:, { ...
            'index', ...
            'latitude_inner_deg', ...
            'latitude_outer_deg', ...
            'mu', ...
            'lambda_local', ...
            'k_mean_local'}));
    end
end
