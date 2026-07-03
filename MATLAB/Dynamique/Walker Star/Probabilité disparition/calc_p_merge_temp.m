function p_merge_t = calc_p_merge_temp(lambda_eff_t, vrel_t, dmax, dt)
% calc_p_merge_temp
% Calcule p_merge(t) a partir de la densite effective lambda_eff(t)
% et de la vitesse relative moyenne v_rel(t).
%
% Formule :
%   p_merge(t) = 1 - exp(-2 lambda_eff(t) dmax v_rel(t) dt)
%
% Unites attendues :
%   lambda_eff_t : satellites / km^2
%   vrel_t       : km / s
%   dmax         : km
%   dt           : s
%
% Sortie :
%   p_merge_t    : probabilite par pas de temps, bornee dans [0,1]

    lambda_eff_t = lambda_eff_t(:);
    vrel_t = vrel_t(:);

    if numel(lambda_eff_t) ~= numel(vrel_t)
        error('lambda_eff_t et vrel_t doivent avoir la meme longueur.');
    end
    if dmax <= 0 || dt <= 0
        error('dmax et dt doivent etre strictement positifs.');
    end

    expo_arg = 2 .* 0.005 * lambda_eff_t .* dmax .* vrel_t .* dt;
    expo_arg = max(expo_arg, 0);

    p_merge_t = 1 - exp(-expo_arg);
    p_merge_t = min(max(p_merge_t, 0), 1);
end
