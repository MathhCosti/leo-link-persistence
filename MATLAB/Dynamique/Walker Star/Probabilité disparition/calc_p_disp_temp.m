function p_disp_t = calc_p_disp_temp(p_merge_t, p_break_t)
% calc_p_disp_temp
% Combine p_merge(t) et p_break(t) pour obtenir p_disp(t).
%
% Formule :
%   p_disp(t) = 1 - (1 - p_merge(t)) (1 - p_break(t))
%
% Sortie :
%   p_disp_t : probabilite par pas de temps, bornee dans [0,1]

    p_merge_t = p_merge_t(:);
    p_break_t = p_break_t(:);

    if numel(p_merge_t) ~= numel(p_break_t)
        error('p_merge_t et p_break_t doivent avoir la meme longueur.');
    end

    p_merge_t = min(max(p_merge_t, 0), 1);
    p_break_t = min(max(p_break_t, 0), 1);

    p_disp_t = 1 - (1 - p_merge_t) .* (1 - p_break_t);
    p_disp_t = min(max(p_disp_t, 0), 1);
end
