function p_break_t = calc_p_break_temp(vrel_t, dmax, dt)
% calc_p_break_temp
% Calcule p_break(t) a partir de la vitesse relative moyenne des liens.
%
% Formule :
%   p_break(t) = 2 v_rel(t) dt / (pi dmax)
%
% Unites attendues :
%   vrel_t : km / s
%   dmax   : km
%   dt     : s
%
% Sortie :
%   p_break_t : probabilite par pas de temps, bornee dans [0,1]

    vrel_t = vrel_t(:);

    if dmax <= 0 || dt <= 0
        error('dmax et dt doivent etre strictement positifs.');
    end

    p_break_t = 2 .* vrel_t .* dt ./ (pi .* dmax);
    p_break_t = min(max(p_break_t, 0), 1);
end
