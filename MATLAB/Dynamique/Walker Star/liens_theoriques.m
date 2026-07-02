function [L_max_strates, details_table] = liens_theoriques(strates_polaires)
%LIENS_MAX_DEPUIS_STRATES_DEGRE_MOYEN Estime Lmax avec la formule du degre moyen uniquement.
%
%   [Lpoiss, Ldet, T] = liens_max_depuis_strates_degre_moyen(strates_polaires)
%
%   Approximation utilisee :
%       pour chaque strate active i, contenant en moyenne mu_i satellites,
%       et de degre moyen local k_i,
%
%           L_i^{pole} ~= (1/2) * mu_i * k_i.
%
%   Le facteur 1/2 evite de compter deux fois les liens par somme des degres.
%   On multiplie ensuite par 2 pour tenir compte des deux poles.
%
%   Cette version ne traite pas la calotte alpha_max/2 comme une clique :
%   toutes les strates sont traitees de la meme maniere par la formule
%   du degre moyen local.
%
%   Remarque :
%   Cette approximation est une approximation champ moyen par strate. Elle ne
%   distingue pas explicitement les liens internes a une strate des liens entre
%   strates voisines.

    if ~isfield(strates_polaires, 'active_table') || isempty(strates_polaires.active_table)
        L_max_strates = 0;
        details_table = table();
        return;
    end

    T = strates_polaires.active_table;

    mu = T.mu;
    k_mean = T.k_mean_local;

    %% Approximation par degre moyen uniquement
    % Par pole :
    %   L_i = 1/2 * mu_i * k_i
    %
    % Les deux versions sont identiques ici, car on n'applique plus de
    % correction poissonienne/deterministe particuliere.
    L_pole     = 0.5 * mu .* k_mean;

    % Deux poles
    L_two_poles     = 2 * L_pole;
    L_max_strates     = sum(L_two_poles);

    %% Table de details
    details_table = T(:, {'index', 'beta_in', 'beta_out', 'beta_in_deg', 'beta_out_deg', ...
                          'mu', 'k_mean_local'});

    details_table.L_pole = L_pole;
    details_table.L_two_poles = L_two_poles;
end
