function [N_bord, N_balayage, K_bord] = nb_bordure_composante(n, lambda, ell, facteur_balayage)
%NB_BORDURE_COMPOSANTE Estime le nombre de points en bordure d'une composante.
%
%   [N_bord, N_balayage, K_bord] = nb_bordure_composante(n, lambda, ell)
%
%   Entrees :
%       n      : taille de la composante connexe, c.-a-d. nombre de points.
%                Peut etre un scalaire ou un vecteur.
%
%       lambda : densite surfacique moyenne de points [points/km^2].
%                Peut etre un scalaire ou un vecteur compatible avec n.
%
%       ell    : epaisseur effective de la couche de bordure [km].
%                Peut etre un scalaire ou un vecteur compatible avec n.
%                Choix raisonne possible : ell = 1/sqrt(lambda), qui
%                correspond a l'espacement inter-satellites local typique.
%
%       facteur_balayage : fraction geometrique de la bordure contribuant
%                          a l'aire nouvellement balayee.
%                          Par defaut : 1/pi, issu de la moyenne de
%                          (cos(theta))_+ sur les orientations.
%
%   Sorties :
%       N_bord      : nombre moyen de points en bordure
%       N_balayage  : nombre moyen de points de bordure contribuant au balayage
%       K_bord      : coefficient K tel que N_bord = K*sqrt(n)
%
%   Modele :
%       N_bord ≈ K_bord * sqrt(n)
%       K_bord = 2*sqrt(pi)*ell*sqrt(lambda)
%
%       N_balayage ≈ facteur_balayage * N_bord
%
%   Remarque :
%       Si ell = 1/sqrt(lambda) et facteur_balayage = 1/pi, alors
%       N_balayage(n) = 2/sqrt(pi) * sqrt(n), et donc
%       chi_merge = N_balayage(n)/n = 2/sqrt(pi*n).

    if nargin < 4 || isempty(facteur_balayage)
        facteur_balayage = 1/pi;
    end

    % Verifications simples, compatibles avec des vecteurs.
    if any(n(:) < 0)
        error('La taille de composante n doit etre positive ou nulle.');
    end

    if any(lambda(:) <= 0)
        error('La densite lambda doit etre strictement positive.');
    end

    if any(ell(:) <= 0)
        error('L''epaisseur ell doit etre strictement positive.');
    end

    if any(facteur_balayage(:) < 0) || any(facteur_balayage(:) > 1)
        error('Le facteur de balayage doit etre compris entre 0 et 1.');
    end

    % Coefficient de bordure. Les operations element par element permettent
    % d'utiliser des scalaires ou des vecteurs compatibles.
    K_bord = 2 .* sqrt(pi) .* ell .* sqrt(lambda);

    % Nombre de points en bordure.
    N_bord = K_bord .* sqrt(n);

    % Nombre de points contribuant a l'aire balayee.
    N_balayage = facteur_balayage .* N_bord;
end
