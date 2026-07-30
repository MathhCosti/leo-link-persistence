function [p_break_t,details] = calc_p_break_th( ...
    vrel_t,lambda_eff_t,E_edges_t,beta0_t,dmax,dt)
%CALC_P_BREAK_TH Formule lineaire du modele spatial ecrit dans le LaTeX.
%
%   p_break(t) =
%       E[|E(t)|] / E[beta0(t)]
%       * 2 v_rel(t) dt / (pi dmax)
%       * exp[-lambda_eff(t) A_inter(dmax)]
%
% avec :
%
%   A_inter(dmax)
%       = dmax^2 * (2*pi/3 - sqrt(3)/2).
%
% Entrees :
%   vrel_t      : vitesse relative theorique [km/s]
%   lambda_eff_t: intensite effective [satellites/km^2]
%   E_edges_t   : nombre moyen theorique de liens
%   beta0_t     : nombre moyen theorique de composantes
%   dmax        : portee maximale [km]
%   dt          : pas temporel [s]
%
% Sorties :
%   p_break_t   : probabilite theorique de rupture par composante
%   details     : variables intermediaires

    validateattributes(vrel_t,{'numeric'}, ...
        {'real','finite','nonnegative'},mfilename,'vrel_t',1);
    validateattributes(lambda_eff_t,{'numeric'}, ...
        {'real','finite','nonnegative'},mfilename,'lambda_eff_t',2);
    validateattributes(E_edges_t,{'numeric'}, ...
        {'real','finite','nonnegative'},mfilename,'E_edges_t',3);
    validateattributes(beta0_t,{'numeric'}, ...
        {'real','finite','positive'},mfilename,'beta0_t',4);
    validateattributes(dmax,{'numeric'}, ...
        {'scalar','real','finite','positive'},mfilename,'dmax',5);
    validateattributes(dt,{'numeric'}, ...
        {'scalar','real','finite','positive'},mfilename,'dt',6);

    input_size = size(vrel_t);

    vrel_t = double(vrel_t(:));
    lambda_eff_t = double(lambda_eff_t(:));
    E_edges_t = double(E_edges_t(:));
    beta0_t = double(beta0_t(:));

    n = numel(vrel_t);

    if numel(lambda_eff_t) ~= n || ...
            numel(E_edges_t) ~= n || ...
            numel(beta0_t) ~= n
        error(['vrel_t, lambda_eff_t, E_edges_t et beta0_t ', ...
               'doivent avoir la meme longueur.']);
    end

    A_inter_at_dmax = ...
        dmax^2 * (2*pi/3 - sqrt(3)/2);

    q_break_link_t = ...
        2 .* vrel_t .* dt ./ (pi*dmax);

    p_bridge_bord_t = ...
        exp(-lambda_eff_t .* A_inter_at_dmax);

    mean_links_per_component_t = ...
        E_edges_t ./ beta0_t;

    p_break_raw_t = ...
        mean_links_per_component_t .* ...
        q_break_link_t .* ...
        p_bridge_bord_t;

    p_break_t = min(max(p_break_raw_t,0),1);
    p_break_t = reshape(p_break_t,input_size);

    if nargout >= 2
        details = struct();
        details.A_inter_at_dmax = A_inter_at_dmax;
        details.q_break_link_t = ...
            reshape(q_break_link_t,input_size);
        details.p_bridge_bord_t = ...
            reshape(p_bridge_bord_t,input_size);
        details.mean_links_per_component_t = ...
            reshape(mean_links_per_component_t,input_size);
        details.p_break_raw_t = ...
            reshape(p_break_raw_t,input_size);
    end
end
