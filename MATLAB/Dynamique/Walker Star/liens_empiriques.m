function [L_incident, n_sat] = liens_empiriques(A, z_t, R, strates_polaires)
%LIENS_EMPIRIQUES_PAR_STRATE Compte les liens empiriques dans les strates polaires.
%
%   [L_incident, n_sat] = liens_empiriques_par_strate(A, z_t, R, strates_polaires)
%
%   Entrées :
%   - A                : matrice d'adjacence du graphe à un instant donné.
%   - z_t              : coordonnées z des satellites à cet instant.
%   - R                : rayon orbital.
%   - strates_polaires : structure retournée par strates_polaires_degre_moyen.
%
%   Sorties, pour chaque strate active et pour les deux pôles réunis :
%   - L_incident(i) : contribution par degré de la strate i, c'est-à-dire
%                     1/2 * somme des degrés des satellites de la strate i.
%                     Cette quantité attribue la moitié d'un lien inter-strates
%                     à chacune des deux strates concernées.
%   - n_sat(i)      : nombre de satellites dans la strate i.
%
%   La distance angulaire au pôle le plus proche est calculée par
%       beta = acos(abs(z)/R).
%   Ainsi, les strates Nord et Sud sont comptées ensemble.

    if ~isfield(strates_polaires, 'active_table') || isempty(strates_polaires.active_table)
        L_incident = zeros(0,1);
        n_sat = zeros(0,1);
        return;
    end

    T = strates_polaires.active_table;
    n_strates = height(T);

    L_incident = zeros(n_strates,1);
    n_sat = zeros(n_strates,1);

    % Distance angulaire au pôle le plus proche.
    z_norm = abs(z_t(:)) ./ R;
    z_norm = max(min(z_norm, 1), 0); % sécurité numérique
    beta_pole = acos(z_norm);

    % Degré empirique de chaque satellite dans le graphe complet.
    deg = sum(A, 2);

    for i = 1:n_strates
        beta_in = T.beta_in(i);
        beta_out = T.beta_out(i);

        if i < n_strates
            mask = (beta_pole >= beta_in) & (beta_pole < beta_out);
        else
            % On inclut la frontière externe pour la dernière strate active.
            mask = (beta_pole >= beta_in) & (beta_pole <= beta_out);
        end

        n_sat(i) = sum(mask);

        % Contribution par degré : 1/2 somme des degrés des satellites de la strate.
        % Cette grandeur est celle qui se compare le mieux au modèle
        % L_i ≈ 1/2 * mu_i * k_i.
        L_incident(i) = 0.5 * sum(deg(mask));
    end
end
