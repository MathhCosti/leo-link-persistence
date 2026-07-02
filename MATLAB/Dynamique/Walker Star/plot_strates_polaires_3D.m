function plot_strates_polaires_3D(R, alpha_max, strates_polaires, varargin)
%PLOT_STRATES_POLAIRES_3D Représente en 3D les strates polaires conservées.
%
%   plot_strates_polaires_3D(R, alpha_max, strates_polaires)
%
%   Affiche la sphère orbitale et les strates polaires actives/conservées
%   issues de strates_polaires_degre_moyen.m. Les strates sont tracées au
%   pôle Nord et au pôle Sud sous forme de couronnes sphériques
%   transparentes.
%
%   Entrées :
%   - R                : rayon orbital en km
%   - alpha_max        : angle de liaison maximal en rad
%   - strates_polaires : structure retournée par strates_polaires_degre_moyen
%
%   Options :
%   - 'show_inactive' : affiche aussi les strates non conservées, en gris.
%                       Défaut false.
%   - 'sphere_alpha'  : transparence de la sphère orbitale. Défaut 0.04.
%   - 'strate_alpha'  : transparence des strates conservées. Défaut 0.28.
%   - 'n_theta'       : résolution angulaire en theta. Défaut 8.
%   - 'n_phi'         : résolution angulaire en phi. Défaut 160.
%
%   Exemple :
%       plot_strates_polaires_3D(R, alpha_max, strates_polaires);

    %% Options
    parser = inputParser;
    parser.addRequired('R', @(x) isnumeric(x) && isscalar(x) && x > 0);
    parser.addRequired('alpha_max', @(x) isnumeric(x) && isscalar(x) && x > 0);
    parser.addRequired('strates_polaires', @isstruct);
    parser.addParameter('show_inactive', false, @(x) islogical(x) || isnumeric(x));
    parser.addParameter('sphere_alpha', 0.04, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
    parser.addParameter('strate_alpha', 0.28, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
    parser.addParameter('n_theta', 8, @(x) isnumeric(x) && isscalar(x) && x >= 2);
    parser.addParameter('n_phi', 160, @(x) isnumeric(x) && isscalar(x) && x >= 20);
    parser.parse(R, alpha_max, strates_polaires, varargin{:});

    show_inactive = logical(parser.Results.show_inactive);
    sphere_alpha  = parser.Results.sphere_alpha;
    strate_alpha  = parser.Results.strate_alpha;
    n_theta       = round(parser.Results.n_theta);
    n_phi         = round(parser.Results.n_phi);

    if ~isfield(strates_polaires, 'all_table') || ~istable(strates_polaires.all_table)
        error('strates_polaires doit contenir le champ all_table retourné par strates_polaires_degre_moyen.m');
    end

    T = strates_polaires.all_table;

    %% Figure
    figure;
    hold on;

    %% Sphère orbitale de référence
    [Xs, Ys, Zs] = sphere(80);
    surf(R*Xs, R*Ys, R*Zs, ...
        'FaceAlpha', sphere_alpha, ...
        'EdgeAlpha', 0.07, ...
        'FaceColor', [0.75 0.75 0.75], ...
        'EdgeColor', [0.45 0.45 0.45]);

    %% Équateur
    phi_eq = linspace(0, 2*pi, 400);
    plot3(R*cos(phi_eq), R*sin(phi_eq), zeros(size(phi_eq)), 'k:', 'LineWidth', 0.8);

    %% Cercle correspondant à alpha_max pour repère
    x_alpha = R*sin(alpha_max)*cos(phi_eq);
    y_alpha = R*sin(alpha_max)*sin(phi_eq);
    z_alpha_N = R*cos(alpha_max)*ones(size(phi_eq));
    z_alpha_S = -R*cos(alpha_max)*ones(size(phi_eq));

    plot3(x_alpha, y_alpha, z_alpha_N, 'k--', 'LineWidth', 0.9);
    plot3(x_alpha, y_alpha, z_alpha_S, 'k--', 'LineWidth', 0.9);

    %% Choix des strates à afficher
    if show_inactive
        rows_to_plot = true(height(T), 1);
    else
        rows_to_plot = T.active;
    end

    if ~any(rows_to_plot)
        warning('Aucune strate conservée à afficher.');
    end

    n_active = max(1, sum(T.active));
    colors_active = lines(n_active);
    active_counter = 0;

    %% Tracé des strates
    phi_grid = linspace(0, 2*pi, n_phi);

    legend_handles = gobjects(0);
    legend_labels = {};

    for r = 1:height(T)
        if ~rows_to_plot(r)
            continue;
        end

        beta_in  = T.beta_in(r);
        beta_out = T.beta_out(r);

        if T.active(r)
            active_counter = active_counter + 1;
            col = colors_active(active_counter, :);
            face_alpha = strate_alpha;
            edge_alpha = 0.15;
            edge_col = col;
            line_style = '-';
        else
            col = [0.6 0.6 0.6];
            face_alpha = 0.08;
            edge_alpha = 0.05;
            edge_col = [0.55 0.55 0.55];
            line_style = ':';
        end

        beta_grid = linspace(beta_in, beta_out, n_theta);
        [PHI, BETA] = meshgrid(phi_grid, beta_grid);

        %% Pôle Nord : beta mesuré depuis le pôle Nord
        XN = R*sin(BETA).*cos(PHI);
        YN = R*sin(BETA).*sin(PHI);
        ZN = R*cos(BETA);

        hN = surf(XN, YN, ZN, ...
            'FaceColor', col, ...
            'FaceAlpha', face_alpha, ...
            'EdgeColor', edge_col, ...
            'EdgeAlpha', edge_alpha);

        %% Pôle Sud : beta mesuré depuis le pôle Sud
        XS = R*sin(BETA).*cos(PHI);
        YS = R*sin(BETA).*sin(PHI);
        ZS = -R*cos(BETA);

        surf(XS, YS, ZS, ...
            'FaceColor', col, ...
            'FaceAlpha', face_alpha, ...
            'EdgeColor', edge_col, ...
            'EdgeAlpha', edge_alpha);

        %% Frontières de la strate, Nord et Sud
        for b = [beta_in beta_out]
            xb = R*sin(b)*cos(phi_eq);
            yb = R*sin(b)*sin(phi_eq);
            zbN = R*cos(b)*ones(size(phi_eq));
            zbS = -R*cos(b)*ones(size(phi_eq));

            plot3(xb, yb, zbN, line_style, 'Color', col, 'LineWidth', 1.1);
            plot3(xb, yb, zbS, line_style, 'Color', col, 'LineWidth', 1.1);
        end

        if T.active(r)
            legend_handles(end+1) = hN; %#ok<AGROW>
            legend_labels{end+1} = sprintf('Strate %d : %.1f°–%.1f°', ...
                T.index(r), T.beta_in_deg(r), T.beta_out_deg(r)); %#ok<AGROW>
        end
    end

    %% Axes et titres
    axis equal;
    grid on;
    xlabel('x (km)');
    ylabel('y (km)');
    zlabel('z (km)');

    title(sprintf(['Strates polaires conservées — \\alpha_{max}=%.2f° ; ', ...
        '\\beta_{stop}=%.2f°'], ...
        rad2deg(alpha_max), rad2deg(strates_polaires.beta_stop)));

    view(35, 22);
    rotate3d on;

    if ~isempty(legend_handles)
        legend(legend_handles, legend_labels, 'Location', 'bestoutside');
    end

    hold off;
end
