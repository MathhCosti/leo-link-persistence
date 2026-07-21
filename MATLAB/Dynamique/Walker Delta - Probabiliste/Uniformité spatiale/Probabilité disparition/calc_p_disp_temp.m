function p_disp_t = calc_p_disp_temp(p_merge_t, p_break_t)
%CALC_P_DISP_TEMP Combine p_merge(t) et p_break(t).
%
% Formule de modification topologique totale :
%
%   p_disp(t) = 1 - (1-p_merge(t))(1-p_break(t)).
%
% Cette formule suppose l'independance des mecanismes de fusion et de
% rupture. Pour la seule mort des barres H0, utiliser plutot p_merge_t.

    input_size = size(p_merge_t);

    p_merge_t = double(p_merge_t(:));
    p_break_t = double(p_break_t(:));

    if numel(p_merge_t) ~= numel(p_break_t)
        error('p_merge_t et p_break_t doivent avoir la meme longueur.');
    end

    if any(~isfinite(p_merge_t)) || any(~isfinite(p_break_t))
        error('Les probabilites doivent etre finies.');
    end

    p_merge_t = min(max(p_merge_t,0),1);
    p_break_t = min(max(p_break_t,0),1);

    p_disp_t = 1 - (1-p_merge_t).*(1-p_break_t);
    p_disp_t = min(max(p_disp_t,0),1);

    p_disp_t = reshape(p_disp_t,input_size);
end
