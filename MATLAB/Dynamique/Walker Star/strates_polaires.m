function strates = strates_polaires(lambda, R, alpha_max, varargin)
%STRATES_POLAIRES_TOUTES Découpe la région polaire en strates jusqu'à pi/2.
%
%   strates = strates_polaires_toutes(lambda, R, alpha_max)
%
%   Cette version ne fait aucune sélection par percolation : toutes les
%   strates entre beta = 0 et beta_max <= pi/2 sont conservées.
%
%   Modèle utilisé :
%   - Une strate polaire est une couronne beta_in <= beta <= beta_out.
%   - La surface réelle de la couronne polaire vaut :
%       A_zone = 2*pi*R^2*(cos(beta_in) - cos(beta_out))
%   - La demi-bande équatoriale initiale qui alimente cette couronne vaut :
%       A_source = 2*pi*R^2*(sin(beta_out) - sin(beta_in))
%   - Le nombre moyen de satellites dans la strate est :
%       mu = lambda * A_source
%   - La densité locale dans la strate au maximum est :
%       lambda_local = mu / A_zone
%   - L'aire angulaire de voisinage d'un lien vaut :
%       A_link = 2*pi*R^2*(1 - cos(alpha_max))
%   - Le degré moyen local informatif vaut :
%       k_mean = lambda_local * A_link
%
%   Options :
%   - 'beta_step' : largeur angulaire des strates. Défaut alpha_max/2.
%   - 'beta_max'  : angle maximal étudié depuis le pôle. Défaut pi/2.
%                   La valeur est automatiquement bornée par pi/2.
%   - 'verbose'   : afficher la table. Défaut false.
%
%   Sortie strates : structure contenant
%   - all_table       : table de toutes les strates
%   - active_table    : identique à all_table, car aucune strate n'est rejetée
%   - keep_table      : identique à all_table
%   - A_link          : aire de voisinage de lien
%   - beta_step       : pas angulaire utilisé
%   - beta_max        : angle maximal réellement utilisé
%   - beta_stop       : frontière externe de la dernière strate

    %% Options
    parser = inputParser;
    parser.addRequired('lambda', @(x) isnumeric(x) && isscalar(x) && x >= 0);
    parser.addRequired('R', @(x) isnumeric(x) && isscalar(x) && x > 0);
    parser.addRequired('alpha_max', @(x) isnumeric(x) && isscalar(x) && x > 0);
    parser.addParameter('beta_step', alpha_max/2, @(x) isnumeric(x) && isscalar(x) && x > 0);
    parser.addParameter('beta_max', pi/2, @(x) isnumeric(x) && isscalar(x) && x > 0);
    parser.addParameter('verbose', false, @(x) islogical(x) || isnumeric(x));
    parser.parse(lambda, R, alpha_max, varargin{:});

    beta_step = parser.Results.beta_step;
    beta_max  = min(parser.Results.beta_max, pi/2);
    verbose   = logical(parser.Results.verbose);

    %% Frontières des strates
    % On ajoute pi/2 comme dernière frontière si le pas ne tombe pas
    % exactement dessus, sans jamais dépasser pi/2.
    beta_edges = 0:beta_step:beta_max;

    if isempty(beta_edges) || beta_edges(end) < beta_max - 1e-12
        beta_edges(end+1) = beta_max; %#ok<AGROW>
    else
        beta_edges(end) = beta_max;
    end

    % Sécurité numérique : aucune frontière ne doit dépasser pi/2.
    beta_edges = min(beta_edges, pi/2);
    beta_edges = unique(beta_edges, 'stable');

    n_strates = length(beta_edges) - 1;

    beta_in  = beta_edges(1:end-1)';
    beta_out = beta_edges(2:end)';
    beta_mid = 0.5*(beta_in + beta_out);

    %% Aires et grandeurs locales
    A_link = 2*pi*R^2*(1 - cos(alpha_max));

    A_zone = 2*pi*R^2*(cos(beta_in) - cos(beta_out));
    A_source = 2*pi*R^2*(sin(beta_out) - sin(beta_in));

    mu_strate = lambda * A_source;
    lambda_locale = mu_strate ./ A_zone;
    k_mean_local = lambda_locale * A_link;

    % Aucune comparaison à la percolation : tout est conservé.
    keep = true(n_strates, 1);
    active_mask = true(n_strates, 1);

    if n_strates > 0
        beta_stop = beta_out(end);
    else
        beta_stop = 0;
    end

    %% Table récapitulative
    index = (1:n_strates)';
    all_table = table(index, beta_in, beta_out, beta_mid, ...
        rad2deg(beta_in), rad2deg(beta_out), ...
        A_source, A_zone, mu_strate, lambda_locale, k_mean_local, keep, active_mask, ...
        'VariableNames', {'index', 'beta_in', 'beta_out', 'beta_mid', ...
        'beta_in_deg', 'beta_out_deg', ...
        'A_source', 'A_zone', 'mu', 'lambda_local', 'k_mean_local', 'keep', 'active'});

    %% Sortie
    strates = struct();
    strates.all_table = all_table;
    strates.keep_table = all_table;
    strates.active_table = all_table;
    strates.A_link = A_link;
    strates.beta_step = beta_step;
    strates.beta_max = beta_max;
    strates.beta_stop = beta_stop;

    %% Affichage optionnel
    if verbose
        fprintf('\n=== Strates polaires conservees jusqu''a pi/2 ===\n');
        fprintf('A_link = %.3e km^2\n', A_link);
        fprintf('beta_step = %.4f rad = %.2f deg\n', beta_step, rad2deg(beta_step));
        fprintf('beta_stop = %.4f rad = %.2f deg\n\n', beta_stop, rad2deg(beta_stop));

        disp(all_table(:, {'index', 'beta_in_deg', 'beta_out_deg', 'mu', 'k_mean_local', 'active'}));
    end
end
