function [p_merge, n_merge, n_components, is_merge] = ...
    calc_p_merge_emp(A0, A1)
%CALC_P_MERGE_EMP Probabilite empirique de fusion entre deux graphes.
%
% Une composante de A0 participe a une fusion si une composante de A1
% qu'elle intersecte contient aussi des sommets issus d'au moins une
% autre composante de A0.

    labels0 = conncomp(graph(A0)).';
    labels1 = conncomp(graph(A1)).';

    n0 = max(labels0);
    n1 = max(labels1);

    overlap = accumarray([labels0, labels1], 1, [n0, n1]);

    is_merge = false(n0,1);

    for c = 1:n0
        target_ids = find(overlap(c,:) > 0);

        for d = target_ids
            if nnz(overlap(:,d) > 0) >= 2
                is_merge(c) = true;
                break;
            end
        end
    end

    n_components = n0;
    n_merge = sum(is_merge);

    if n_components > 0
        p_merge = n_merge / n_components;
    else
        p_merge = NaN;
    end
end
