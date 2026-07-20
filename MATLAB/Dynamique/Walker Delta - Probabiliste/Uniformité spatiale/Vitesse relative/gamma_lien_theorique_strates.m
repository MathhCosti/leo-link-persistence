function [E_sin_half_th, gamma_eff_th, gamma_eff_th_deg, vrel_th, ...
          alt_centers_deg, s_alt_th, Phi_alt, link_counts_alt_time] = ...
          gamma_lien_theorique_strates(link_counts_alt_time, alt_edges_deg, R, dmax, v_orb)
%GAMMA_LIEN_THEORIQUE_STRATES Approximation theorique de gamma_eff(t) par strates.
%
% Cette fonction ne recalcule pas le nombre de liens. Elle prend en entree
% le nombre de liens observes ou modelises dans chaque strate d'altitude :
%
%   link_counts_alt_time(k,m) = nombre de liens au temps k dans la strate m.
%
% Elle calcule ensuite un facteur angulaire theorique par strate
%
%   s_m = E[ sin(gamma/2) | lien, altitude dans strate m ]
%
% a partir d'une approximation geometrique locale :
%
%   Phi_m = min(pi, alpha_max / cos(ell_m)),
%
%   s_m = sin(ell_m) * (2/Phi_m) * (1 - cos(Phi_m/2)).
%
% La variation temporelle theorique vient alors uniquement du poids des
% strates :
%
%   E_sin_half_th(t) = sum_m L_m(t) s_m / sum_m L_m(t).
%
% Puis :
%
%   gamma_eff_th(t) = 2 asin(E_sin_half_th(t)),
%   vrel_th(t)      = 2 v_orb E_sin_half_th(t).
%
% Entrees :
%   link_counts_alt_time : matrice Nt x M
%   alt_edges_deg        : bords des strates en degres, taille 1 x (M+1)
%   R                    : rayon orbital en km
%   dmax                 : distance maximale de lien en km
%   v_orb                : vitesse orbitale en km/s
%
% Sorties :
%   E_sin_half_th        : serie temporelle theorique de E[sin(gamma/2)]
%   gamma_eff_th         : angle effectif theorique en radians
%   gamma_eff_th_deg     : angle effectif theorique en degres
%   vrel_th              : vitesse relative theorique des liens en km/s
%   alt_centers_deg      : centres des strates en degres
%   s_alt_th             : facteur angulaire theorique par strate
%   Phi_alt              : ouverture longitudinale effective par strate en rad
%   link_counts_alt_time : renvoye tel quel pour sauvegarde

    if nargin < 5
        error('Usage: gamma_lien_theorique_strates(link_counts_alt_time, alt_edges_deg, R, dmax, v_orb)');
    end

    alt_edges_deg = alt_edges_deg(:).';
    alt_centers_deg = 0.5 * (alt_edges_deg(1:end-1) + alt_edges_deg(2:end));
    ell = deg2rad(alt_centers_deg);

    % Angle central maximal de communication.
    alpha_max = 2 * asin(min(1, dmax/(2*R)));

    % Ouverture effective en longitude dans une strate.
    % Pres de l'equateur : Phi ~ alpha_max.
    % Pres des poles : cos(ell) -> 0 donc Phi est tronque a pi.
    cosell = max(cos(ell), 1e-12);
    Phi_alt = min(pi, alpha_max ./ cosell);

    % Facteur angulaire theorique par strate :
    % s_m = E[sin(gamma/2) | lien, strate m]
    s_alt_th = sin(ell) .* (2 ./ Phi_alt) .* (1 - cos(Phi_alt/2));

    % Evite d'eventuelles valeurs hors [0,1] dues aux approximations.
    s_alt_th = max(0, min(1, s_alt_th));

    Ltot = sum(link_counts_alt_time, 2);
    Nt = size(link_counts_alt_time, 1);

    E_sin_half_th = NaN(Nt,1);

    valid = Ltot > 0;
    if any(valid)
        E_sin_half_th(valid) = (link_counts_alt_time(valid,:) * s_alt_th(:)) ./ Ltot(valid);
    end

    E_sin_half_th = max(0, min(1, E_sin_half_th));

    gamma_eff_th = 2 * asin(E_sin_half_th);
    gamma_eff_th_deg = rad2deg(gamma_eff_th);
    vrel_th = 2 * v_orb * E_sin_half_th;
end
