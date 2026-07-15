function [p_merge_t, chi_merge_t, E_t, beta0_geom] = calc_p_merge_temp(lambda_eff_t, vrel_t, dmax, dt, liens_file)
% calc_p_merge_temp
% Calcule p_merge(t) a partir de la densite effective lambda_eff(t)
% et de la vitesse relative moyenne v_rel(t), avec correction topologique
% par fraction de satellites en bordure des composantes.
%
% Formule brute :
%   A_sweep(t)  = 2 dmax v_rel(t) dt / pi
%   p_merge(t) = 1 - exp(-lambda_eff(t) A_sweep(t))
%
% Correction de merge :
%   beta0_geom ~= 1 + N1 + N2 + N3
%   n_comp     ~= N / beta0_geom
%
% Le nombre de satellites en bordure est calcule via la fonction externe :
%   nb_bordure_composante(n, lambda, ell, facteur_balayage)
%
% puis :
%   chi_merge(t) = N_balayage(n_comp, lambda_eff(t)) / n_comp
%
% Cette correction signifie que seuls les satellites en bordure d'une
% composante contribuent a l'aire nouvellement balayee, et qu'en moyenne
% seule une fraction facteur_balayage de la bordure est orientee favorablement.
%
% Approximation beta0 geom :
%   N1 = N(1-p)^(N-1)
%   N2 = C(N,2) p (1-c2 p)^(N-2)
%   N3 = C(N,3) c3_conn p^2 (1-c3_union p)^(N-3)
%
% Unites attendues :
%   lambda_eff_t : satellites / km^2
%   vrel_t       : km / s
%   dmax         : km
%   dt           : s
%   liens_file   : fichier .mat contenant N_mean_theory et p_link_uniform,
%                  ou les variables permettant de les reconstruire.
%
% Sorties :
%   p_merge_t   : probabilite par pas de temps, bornee dans [0,1]
%   chi_merge_t : facteur correctif de merge par bordure de composante
%   E_t         : nombre de liens temporel lu/interpole si disponible
%                 retourne NaN sinon. Il n'est plus utilise dans p_merge.
%   beta0_geom  : beta0 theorique utilise, avec isoles + dimeres + trimeres

    if nargin < 5 || isempty(liens_file)
        script_dir = fileparts(mfilename('fullpath'));
        liens_file = fullfile(script_dir, '..', 'Nombre liens', 'liens_inter.mat');
    end

    lambda_eff_t = lambda_eff_t(:);
    vrel_t = vrel_t(:);

    if numel(lambda_eff_t) ~= numel(vrel_t)
        error('lambda_eff_t et vrel_t doivent avoir la meme longueur.');
    end
    if any(lambda_eff_t < 0)
        error('lambda_eff_t doit etre positive ou nulle.');
    end
    if dmax <= 0 || dt <= 0
        error('dmax et dt doivent etre strictement positifs.');
    end

    % S'assure que la fonction nb_bordure_composante.m est visible si elle
    % est placee dans le meme dossier que ce fichier.
    script_dir = fileparts(mfilename('fullpath'));
    if ~isempty(script_dir)
        addpath(script_dir);
    end

    if exist('nb_bordure_composante', 'file') ~= 2
        error(['La fonction nb_bordure_composante.m est introuvable. ', ...
               'Place-la dans le meme dossier que calc_p_merge_temp.m ', ...
               'ou ajoute son dossier au path MATLAB.']);
    end

    n_target = numel(vrel_t);

    % Aire balayee brute par un satellite pendant dt.
    % Le facteur 1/pi correspond a l'approximation de la composante radiale
    % effective de la vitesse relative utilisee dans le modele.
    A_sweep = 2 .* dmax .* vrel_t .* dt ./ pi;

    % Chargement des parametres globaux et calcul de beta0 avec c2/c3.
    [N, p_link, E_t, beta0_geom, N1_theory, N2_theory, N3_theory] = ...
        local_params_beta0_geom(liens_file, n_target);

    % Taille moyenne d'une composante.
    beta0_geom_safe = max(beta0_geom, eps);
    n_comp_mean = N ./ beta0_geom_safe;
    n_comp_mean = max(n_comp_mean, 1);

    % Facteur de bordure pour les fusions.
    % Choix de modelisation non arbitraire : ell_bord(t) est l'espacement
    % inter-satellites local typique, soit ell = 1/sqrt(lambda_eff(t)).
    % Avec facteur_balayage = 1/pi, on obtient analytiquement :
    %   N_balayage(n) = 2/sqrt(pi) * sqrt(n)
    %   chi_merge = 2/sqrt(pi*n).
    facteur_balayage = 1/pi;

    N_bord_mean_t = zeros(n_target, 1);
    N_balayage_mean_t = zeros(n_target, 1);
    K_bord_t = zeros(n_target, 1);
    ell_bord_t = NaN(n_target, 1);

    for k = 1:n_target
        if lambda_eff_t(k) > 0
            ell_bord_t(k) = 1 ./ (2*sqrt(lambda_eff_t(k)));
            [N_bord_mean_t(k), N_balayage_mean_t(k), K_bord_t(k)] = ...
                nb_bordure_composante(n_comp_mean, lambda_eff_t(k), ell_bord_t(k), facteur_balayage);
        else
            N_bord_mean_t(k) = 0;
            N_balayage_mean_t(k) = 0;
            K_bord_t(k) = 0;
        end
    end

    % Rapport : satellites utiles au balayage / satellites dans la composante.
    % Avec ell = 1/sqrt(lambda_eff), ce rapport est independant de lambda_eff
    % sauf pour les instants ou lambda_eff vaut zero.
    chi_merge_t = N_balayage_mean_t ./ n_comp_mean;

    % Un facteur correctif est une fraction effective de satellites utiles.
    chi_merge_t = min(max(chi_merge_t, 0), 1);

    % Correction appliquee a l'aire balayee, et non directement a p_merge.
    A_sweep_corr = chi_merge_t .* A_sweep;

    expo_arg = 0.5 * lambda_eff_t .* A_sweep_corr;
    expo_arg = max(expo_arg, 0);

    p_merge_t = 1 - exp(-expo_arg);
    p_merge_t = min(max(p_merge_t, 0), 1);

    fprintf('\n--- calc_p_merge_temp ---\n');
    fprintf('N utilise : %.3f\n', N);
    fprintf('p_link utilise : %.6f\n', p_link);
    fprintf('N1 theorie : %.3f\n', N1_theory);
    fprintf('N2 theorie geom : %.3f\n', N2_theory);
    fprintf('N3 theorie geom : %.3f\n', N3_theory);
    fprintf('beta0 geom utilise : %.3f\n', beta0_geom);
    fprintf('taille moyenne composante N/beta0 : %.3f\n', n_comp_mean);
    fprintf('chi_merge theorique simplifie 2/sqrt(pi*n) : %.6f\n', ...
        2 / sqrt(pi * n_comp_mean));
    fprintf('ell_bord(t) min / mean / max : %.3f / %.3f / %.3f km\n', ...
        min(ell_bord_t,[],'omitnan'), mean(ell_bord_t,'omitnan'), max(ell_bord_t,[],'omitnan'));
    fprintf('facteur_balayage : %.3f\n', facteur_balayage);
    fprintf('N_balayage moyen min / mean / max : %.6f / %.6f / %.6f\n', ...
        min(N_balayage_mean_t), mean(N_balayage_mean_t, 'omitnan'), max(N_balayage_mean_t));
    fprintf('chi_merge min / mean / max : %.6f / %.6f / %.6f\n', ...
        min(chi_merge_t), mean(chi_merge_t, 'omitnan'), max(chi_merge_t));
end

function [N, p_link, E_t, beta0_geom, N1_theory, N2_theory, N3_theory] = ...
    local_params_beta0_geom(liens_file, n_target)
% Charge N et p_link, puis construit beta0_geom = 1 + N1 + N2 + N3.
% E_t est conserve uniquement pour diagnostic / compatibilite.

    if ~isfile(liens_file)
        error('Fichier %s introuvable. Place liens_inter.mat dans le dossier courant ou passe son chemin en argument.', liens_file);
    end

    data = load(liens_file);

    % Nombre de liens temporel : conserve pour diagnostic, pas utilise pour p_merge.
    if isfield(data, 'mean_edges')
        E_raw = double(data.mean_edges(:));
    elseif isfield(data, 'num_edges_all')
        E_raw = mean(double(data.num_edges_all), 1, 'omitnan').';
    else
        E_raw = NaN;
    end

    if isfield(data, 'N_mean_theory')
        N = double(data.N_mean_theory);
    elseif isfield(data, 'N_all')
        N = mean(double(data.N_all(:)), 'omitnan');
    elseif isfield(data, 'lambda') && isfield(data, 'surface_sphere')
        N = double(data.lambda) * double(data.surface_sphere);
    else
        error('Impossible de determiner N depuis %s.', liens_file);
    end

    if isfield(data, 'p_link_uniform')
        p_link = double(data.p_link_uniform);
    elseif isfield(data, 'dmax') && isfield(data, 'R')
        alpha_max = 2 * asin(min(double(data.dmax) / (2 * double(data.R)), 1));
        p_link = (1 - cos(alpha_max)) / 2;
    else
        error('Impossible de determiner p_link depuis %s.', liens_file);
    end

    p_link = min(max(p_link, 0), 1);

    % Interpolation des liens si disponibles.
    if all(isnan(E_raw))
        E_t = NaN(n_target, 1);
    elseif numel(E_raw) == n_target
        E_t = E_raw;
    else
        x_raw = linspace(0, 1, numel(E_raw));
        x_tgt = linspace(0, 1, n_target);
        E_t = interp1(x_raw, E_raw, x_tgt, 'linear', 'extrap').';
    end

    % Constantes geometriques.
    c2_union = 1 + 3*sqrt(3)/(4*pi);  % ~= 1.4135, aire union dimere
    c3_conn  = 1 + 3*sqrt(3)/(2*pi);  % ~= 1.827, connexite interne triplet
    c3_union = 1.80;                  % coefficient effectif d'aire d'union

    q1_ext = max(1 - p_link, 0);
    q2_ext = max(1 - c2_union*p_link, 0);
    q3_ext = max(1 - c3_union*p_link, 0);

    % Composantes de taille 1.
    N1_theory = N * q1_ext^(N - 1);

    % Composantes de taille 2.
    if N >= 2
        N2_theory = local_nchoosek_real(N, 2) * p_link * q2_ext^(N - 2);
    else
        N2_theory = 0;
    end

    % Composantes de taille 3.
    if N >= 3
        p_conn_3_geom = c3_conn * p_link^2;
        p_conn_3_geom = min(max(p_conn_3_geom, 0), 1);
        N3_theory = local_nchoosek_real(N, 3) * p_conn_3_geom * q3_ext^(N - 3);
    else
        N3_theory = 0;
    end

    beta0_geom = 1 + N1_theory + N2_theory + N3_theory;
    beta0_geom = min(max(beta0_geom, 1), N);
end

function c = local_nchoosek_real(n, k)
% Extension numeriquement stable de nchoosek pour n non entier eventuel.
% Pour n entier, donne la meme valeur que nchoosek(n,k).
    if n < k
        c = 0;
    else
        c = exp(gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1));
    end
end
