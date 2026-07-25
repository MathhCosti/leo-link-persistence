function [G, positions, A, params] = graphe_3D_delta(lambda, inc_deg, dmax, t, h, rng_seed)
%GRAPHE_3D_DELTA Genere et affiche un graphe LEO Delta stochastique.
%
%   [G, positions, A, params] = graphe_3D_delta(lambda, inc_deg, dmax, t, h)
%
%   Le modele utilise :
%     - des orbites circulaires de rayon R = R_earth + h ;
%     - une inclinaison commune inc_deg ;
%     - un RAAN Omega uniforme sur [0,2*pi[ pour chaque satellite ;
%     - une phase orbitale initiale u0 uniforme sur [0,2*pi[ ;
%     - un mouvement u(t) = u0 + omega*t.
%
%   Entrees (toutes optionnelles) :
%     lambda   : densite satellitaire sur la sphere orbitale [sat/km^2]
%                defaut : 4e-7
%     inc_deg  : inclinaison commune [deg]
%                defaut : 58
%     dmax     : distance maximale d'un lien ISL [km]
%                defaut : 1500
%     t        : instant du graphe [s]
%                defaut : 0
%     h        : altitude orbitale [km]
%                defaut : 550
%     rng_seed : graine aleatoire ; [] conserve l'etat courant
%                defaut : []
%
%   Sorties :
%     G         : objet MATLAB graph
%     positions : matrice N x 3 des positions satellitaires [km]
%     A         : matrice d'adjacence sparse
%     params    : structure contenant les parametres et tirages orbitaux
%
%   Exemple :
%     [G,pos,A,p] = graphe_3D_delta(4e-7,58,1500,0,550,3);

    %% Valeurs par defaut
    if nargin < 1 || isempty(lambda),   lambda = 4e-7; end
    if nargin < 2 || isempty(inc_deg),  inc_deg = 58; end
    if nargin < 3 || isempty(dmax),     dmax = 2500; end
    if nargin < 4 || isempty(t),        t = 0; end
    if nargin < 5 || isempty(h),        h = 550; end
    if nargin < 6,                      rng_seed = []; end

    validateattributes(lambda,  {'numeric'}, {'scalar','real','finite','positive'});
    validateattributes(inc_deg, {'numeric'}, {'scalar','real','finite','>=',0,'<=',90});
    validateattributes(dmax,    {'numeric'}, {'scalar','real','finite','positive'});
    validateattributes(t,       {'numeric'}, {'scalar','real','finite','nonnegative'});
    validateattributes(h,       {'numeric'}, {'scalar','real','finite','positive'});

    if ~isempty(rng_seed)
        rng(rng_seed);
    end

    %% Parametres physiques
    R_earth = 6371;             % km
    R = R_earth + h;            % km
    mu = 398600;                % km^3/s^2
    omega = sqrt(mu/R^3);       % rad/s
    inc = deg2rad(inc_deg);

    %% Nombre de satellites du PPP sur la sphere orbitale
    surface_sphere = 4*pi*R^2;
    N_mean = lambda*surface_sphere;
    N = poissrnd(N_mean);

    if N < 2
        error(['Le tirage a produit N = %d satellite(s). ', ...
               'Augmenter lambda ou changer la graine aleatoire.'], N);
    end

    %% Modele Delta stochastique en uniformite orbitale
    Omega = 2*pi*rand(N,1);     % RAAN aleatoires
    u0 = 2*pi*rand(N,1);        % phases initiales aleatoires
    u_t = mod(u0 + omega*t, 2*pi);

    positions = walker_delta_positions(R, inc, Omega, u_t);
    x = positions(:,1);
    y = positions(:,2);
    z = positions(:,3);

    %% Graphe des liens intersatellites
    D = squareform(pdist(positions));
    A = sparse((D <= dmax) & (D > 0));
    G = graph(A);

    %% Preparation des segments des liens
    [row, col] = find(triu(A,1));
    E = numel(row);

    Xlinks = NaN(3*E,1);
    Ylinks = NaN(3*E,1);
    Zlinks = NaN(3*E,1);

    Xlinks(1:3:end) = x(row);
    Xlinks(2:3:end) = x(col);
    Ylinks(1:3:end) = y(row);
    Ylinks(2:3:end) = y(col);
    Zlinks(1:3:end) = z(row);
    Zlinks(2:3:end) = z(col);

    %% Affichage 3D statique
    figure;
    hold on;

    % Liens ISL
    plot3(Xlinks,Ylinks,Zlinks,'k-','LineWidth',0.5, ...
        'DisplayName','Liens ISL');

    % Satellites
    scatter3(x,y,z,35,'filled','DisplayName','Satellites');

    axis equal;
    grid on;
    xlabel('x (km)');
    ylabel('y (km)');
    zlabel('z (km)');
    title(sprintf(['Graphe LEO Delta stochastique | i = %.1f deg | ', ...
                   't = %.0f s | N = %d | E = %d'], ...
                   inc_deg,t,N,E));
    view(3);
    rotate3d on;
    hold off;

    %% Parametres retournes
    params = struct();
    params.R_earth = R_earth;
    params.h = h;
    params.R = R;
    params.mu = mu;
    params.omega = omega;
    params.lambda = lambda;
    params.N_mean = N_mean;
    params.N = N;
    params.inc_deg = inc_deg;
    params.inc = inc;
    params.dmax = dmax;
    params.t = t;
    params.Omega = Omega;
    params.u0 = u0;
    params.u_t = u_t;
    params.num_edges = E;

    fprintf('Modele Delta stochastique : N = %d, E = %d, i = %.1f deg, t = %.0f s.\n', ...
        N,E,inc_deg,t);
end

function positions = walker_delta_positions(R,inc,Omega,u)
%WALKER_DELTA_POSITIONS Positions cartesiennes d'orbites circulaires.
    x = R*(cos(Omega).*cos(u) - sin(Omega).*sin(u).*cos(inc));
    y = R*(sin(Omega).*cos(u) + cos(Omega).*sin(u).*cos(inc));
    z = R*(sin(u).*sin(inc));
    positions = [x y z];
end
