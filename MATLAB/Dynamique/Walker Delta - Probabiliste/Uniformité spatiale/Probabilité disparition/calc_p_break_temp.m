function [p_break_t,details] = calc_p_break_temp(vrel_t,dmax,dt)
%CALC_P_BREAK_TEMP Probabilite theorique temporelle de rupture Walker Delta.
%
% Version autonome : aucun fichier .mat n'est charge.
%
% Les facteurs topologiques du Walker Delta sont directement renseignes
% dans la section "Parametres topologiques".
%
% Modele :
%
%   q_link(t) = 2 (v_rel(t)/pi) dt / dmax
%
%   q_bridge(t)
%     = q_link(t) * p_bridge_bord / chi_bridge
%
%   p_break(t)
%     = 1 - (1-q_bridge(t))^Bbar
%
% avec :
%   chi_bridge    = P(pont | lien)
%   p_bridge_bord = P(pont | D=dmax,lien)
%   Bbar          = nombre moyen de ponts par composante
%
% Entrees :
%   vrel_t : vitesse relative par transition [km/s]
%   dmax   : portee maximale [km]
%   dt     : pas temporel [s]
%
% Sorties :
%   p_break_t : probabilite theorique temporelle de rupture
%   details   : structure contenant les valeurs intermediaires

    %% Verification des entrees
    validateattributes(vrel_t,{'numeric'}, ...
        {'real','finite','nonnegative'},mfilename,'vrel_t',1);
    validateattributes(dmax,{'numeric'}, ...
        {'scalar','real','finite','positive'},mfilename,'dmax',2);
    validateattributes(dt,{'numeric'}, ...
        {'scalar','real','finite','positive'},mfilename,'dt',3);

    input_size = size(vrel_t);
    vrel_t = double(vrel_t);

    %% =========================================================
    %  Parametres topologiques du Walker Delta
    %
    %  Valeurs obtenues pour :
    %    N = 204
    %    R = 6921 km
    %    inclinaison = 58 deg
    %    dmax = 1500 km
    %
    %  A modifier si les parametres du reseau changent.
    %% =========================================================

    chi_bridge = 0.1065318529;

    p_bridge_bord = 0.1689645497;

    mean_bridges_per_component = 1.1449766365;

    %% Verification de coherence
    if chi_bridge <= 0 || chi_bridge > 1
        error('chi_bridge doit appartenir a ]0,1].');
    end

    if p_bridge_bord < 0 || p_bridge_bord > 1
        error('p_bridge_bord doit appartenir a [0,1].');
    end

    if mean_bridges_per_component < 0
        error('Le nombre moyen de ponts doit etre positif ou nul.');
    end

    %% Probabilite de rupture d'un lien quelconque
    %
    % La vitesse radiale sortante moyenne est approximee par v_rel/pi.
    v_rad_t = vrel_t/pi;

    q_break_link_t = ...
        2 .* v_rad_t .* dt ./ dmax;

    q_break_link_t = ...
        min(max(q_break_link_t,0),1);

    %% Conditionnement par le fait que le lien soit un pont
    bridge_boundary_factor = ...
        p_bridge_bord / chi_bridge;

    q_break_bridge_t = ...
        q_break_link_t .* bridge_boundary_factor;

    q_break_bridge_t = ...
        min(max(q_break_bridge_t,0),1-eps);

    %% Rupture d'une composante
    p_break_t = ...
        1 - (1-q_break_bridge_t).^mean_bridges_per_component;

    p_break_t = ...
        min(max(p_break_t,0),1-eps);

    p_break_t = reshape(p_break_t,input_size);

    %% Diagnostics
    if nargout >= 2
        details = struct();

        details.chi_bridge = chi_bridge;
        details.p_bridge_bord = p_bridge_bord;
        details.bridge_boundary_factor = bridge_boundary_factor;
        details.mean_bridges_per_component = ...
            mean_bridges_per_component;

        details.v_rad_t = reshape(v_rad_t,input_size);
        details.q_break_link_t = ...
            reshape(q_break_link_t,input_size);
        details.q_break_bridge_t = ...
            reshape(q_break_bridge_t,input_size);

        details.p_break_linear_t = reshape( ...
            min(max(mean_bridges_per_component .* ...
                    q_break_bridge_t,0),1),input_size);
    end

    %% Affichage
    fprintf('\n--- calc_p_break_temp autonome ---\n');
    fprintf('chi_bridge                         : %.10f\n', ...
        chi_bridge);
    fprintf('p_bridge_bord                      : %.10f\n', ...
        p_bridge_bord);
    fprintf('Facteur conditionnel frontiere     : %.10f\n', ...
        bridge_boundary_factor);
    fprintf('Ponts moyens par composante        : %.10f\n', ...
        mean_bridges_per_component);
    fprintf('p_break min / moyenne / max        : %.6f / %.6f / %.6f\n', ...
        min(p_break_t,[],'all'), ...
        mean(p_break_t,'all','omitnan'), ...
        max(p_break_t,[],'all'));
end
