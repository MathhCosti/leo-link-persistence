function [p_break, n_break, n_components, is_break] = ...
    calc_p_break_emp(A0, A1)
%CALC_P_BREAK_EMP Probabilite empirique de rupture entre deux graphes.
%
% Une composante de A0 est dite rompue si ses sommets appartiennent
% a au moins deux composantes distinctes dans A1.

    labels0 = conncomp(graph(A0)).';
    labels1 = conncomp(graph(A1)).';

    n0 = max(labels0);
    n1 = max(labels1);

    overlap = accumarray([labels0, labels1], 1, [n0, n1]);

    is_break = false(n0,1);

    for c = 1:n0
        is_break(c) = nnz(overlap(c,:) > 0) >= 2;
    end

    n_components = n0;
    n_break = sum(is_break);

    if n_components > 0
        p_break = n_break / n_components;
    else
        p_break = NaN;
    end
end
