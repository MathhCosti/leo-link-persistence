function [p_break_t, chi_t, E_t, beta0_iso, p_break_raw] = calc_p_break_temp(vrel_t, dmax, dt, liens_file)
% calc_p_break_temp
% Calcule p_break(t) a partir de la vitesse relative moyenne des liens,
% avec correction topologique temporelle issue du nombre de liens contenu
% dans liens_inter.mat.
%
% Formule brute :
%   p_break_raw(t) = 2 v_rel(t) dt / (pi dmax)
%
% Correction :
%   chi(t) = (N - beta0_iso) / |E(t)|
%   p_break_corr(t) = p_break_raw(t) chi(t)
%
% beta0_iso est estime par la theorie des satellites isoles :
%   beta0_iso = 1 + N (1 - p_link)^(N-1)
%
% |E(t)| est lu dans liens_inter.mat, variable mean_edges.
%
% Unites attendues :
%   vrel_t     : km / s
%   dmax       : km
%   dt         : s
%   liens_file : fichier .mat contenant mean_edges, N_mean_theory, p_link_uniform
%
% Sorties :
%   p_break_t   : probabilite corrigee par pas de temps, bornee dans [0,1]
%   chi_t       : facteur correctif temporel
%   E_t         : nombre de liens temporel utilise
%   beta0_iso   : beta0 theorique isoles utilise
%   p_break_raw : probabilite brute non corrigee

    script_dir = fileparts(mfilename('fullpath'));
    liens_file = fullfile(script_dir, '..', 'Nombre liens', 'liens_inter.mat');
    vrel_t = vrel_t(:);

    if dmax <= 0 || dt <= 0
        error('dmax et dt doivent etre strictement positifs.');
    end

    p_break_raw = 2 .* vrel_t.* dt ./ (pi .* dmax * pi);
    p_break_raw = min(max(p_break_raw, 0), 1);

    [chi_t, E_t, beta0_iso] = local_chi_from_liens(liens_file, numel(vrel_t));

    p_break_t = p_break_raw .* chi_t;
    p_break_t = min(max(p_break_t, 0), 1);
end

function [chi_t, E_t, beta0_iso] = local_chi_from_liens(liens_file, n_target)
% Charge le nombre de liens temporel et construit
% chi(t) = (N - beta0_iso) / |E(t)|.

    if ~isfile(liens_file)
        error('Fichier %s introuvable. Place liens_inter.mat dans le dossier courant ou passe son chemin en argument.', liens_file);
    end

    data = load(liens_file);

    if isfield(data, 'mean_edges')
        E_raw = double(data.mean_edges(:));
    elseif isfield(data, 'num_edges_all')
        E_raw = mean(double(data.num_edges_all), 1, 'omitnan').';
    else
        error('Le fichier %s doit contenir mean_edges ou num_edges_all.', liens_file);
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
        alpha_max = 2 * asin(double(data.dmax) / (2 * double(data.R)));
        p_link = (1 - cos(alpha_max)) / 2;
    else
        error('Impossible de determiner p_link depuis %s.', liens_file);
    end

    beta0_iso = 1 + N * (1 - p_link)^(N - 1);

    % Interpolation si les longueurs temporelles ne coincident pas.
    if numel(E_raw) == n_target
        E_t = E_raw;
    else
        x_raw = linspace(0, 1, numel(E_raw));
        x_tgt = linspace(0, 1, n_target);
        E_t = interp1(x_raw, E_raw, x_tgt, 'linear', 'extrap').';
    end

    % Evite division par zero et bornage physique.
    E_t = max(E_t, eps);
    chi_t = (N - beta0_iso) ./ E_t;
    chi_t = min(max(chi_t, 0), 1);
end
