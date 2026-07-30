function [p_merge_t, chi_merge_t, E_t, beta0_geom] = ...
    calc_p_merge_th(lambda_eff_t, vrel_t, dmax, dt, liens_file)
%CALC_P_MERGE_TEMP Probabilite theorique temporelle de fusion.
%
% Cette version utilise le facteur topologique
%
%   chi_merge(t) = (N - beta0_geom) / E(t),
%
% ou E(t) est le nombre moyen de liens au temps t. Si E(t) n'est pas
% disponible dans le fichier, on utilise
%
%   E_th = N(N-1)p_link/2.
%
% La probabilite de fusion est ensuite
%
%   p_merge(t) = 1 - exp[-2 dmax lambda_eff(t)
%                        chi_merge(t) v_rel(t) dt].
%
% Entrees :
%   lambda_eff_t : densite effective par transition [satellites/km^2]
%   vrel_t       : vitesse relative par transition [km/s]
%   dmax         : portee maximale [km]
%   dt           : pas temporel [s]
%   liens_file   : optionnel, chemin vers le fichier contenant N, p_link
%                  et eventuellement le nombre moyen de liens temporel.
%
% Sorties :
%   p_merge_t    : probabilite temporelle de fusion
%   chi_merge_t  : facteur topologique temporel
%   E_t          : nombre de liens utilise
%   beta0_geom   : approximation theorique de beta0

    %% Verification des entrees
    validateattributes(lambda_eff_t,{'numeric'}, ...
        {'real','finite','nonnegative'},mfilename,'lambda_eff_t',1);
    validateattributes(vrel_t,{'numeric'}, ...
        {'real','finite','nonnegative'},mfilename,'vrel_t',2);
    validateattributes(dmax,{'numeric'}, ...
        {'scalar','real','finite','positive'},mfilename,'dmax',3);
    validateattributes(dt,{'numeric'}, ...
        {'scalar','real','finite','positive'},mfilename,'dt',4);

    input_size = size(lambda_eff_t);

    lambda_eff_t = double(lambda_eff_t(:));
    vrel_t = double(vrel_t(:));

    if numel(lambda_eff_t) ~= numel(vrel_t)
        error('lambda_eff_t et vrel_t doivent avoir la meme longueur.');
    end

    n_target = numel(lambda_eff_t);

    %% Localisation du fichier de liens
    if nargin < 5 || isempty(liens_file)
        script_dir = fileparts(mfilename('fullpath'));
        liens_file = fullfile(script_dir, '..', 'Paramètres', 'Nombre liens', 'liens_quadrature_results.mat');
    end

    %% Recuperation des parametres topologiques
    [N, p_link, E_t, beta0_geom, ...
        N1_theory, N2_theory, N3_theory, edge_source] = ...
        local_params_beta0_geom(liens_file,n_target);

    %% Facteur topologique de fusion
    numerator = max(N-beta0_geom,0);

    chi_merge_t = zeros(n_target,1);

    valid_edges = isfinite(E_t) & (E_t > 0);
    chi_merge_t(valid_edges) = numerator ./ E_t(valid_edges);

    chi_merge_t = min(max(chi_merge_t,0),1);

    %% Probabilite de fusion
    exponent = ...
        2 .* dmax .* lambda_eff_t .* ...
        chi_merge_t .* vrel_t .* dt;

    exponent = max(exponent,0);

    p_merge_t = 1-exp(-exponent);
    p_merge_t = min(max(p_merge_t,0),1);

    %% Conservation de la forme des entrees
    p_merge_t = reshape(p_merge_t,input_size);
    chi_merge_t = reshape(chi_merge_t,input_size);

    %% Affichage
    fprintf('\n--- calc_p_merge_temp : facteur topologique ---\n');
    fprintf('N utilise                         : %.6f\n',N);
    fprintf('p_link utilise                    : %.10f\n',p_link);
    fprintf('E[N1]                             : %.6f\n',N1_theory);
    fprintf('E[N2]                             : %.6f\n',N2_theory);
    fprintf('E[N3]                             : %.6f\n',N3_theory);
    fprintf('beta0_geom                        : %.6f\n',beta0_geom);
    fprintf('N-beta0_geom                      : %.6f\n',numerator);
    fprintf('Source de E(t)                    : %s\n',edge_source);
    fprintf('E(t) min / moyenne / max          : %.6f / %.6f / %.6f\n', ...
        min(E_t,[],'omitnan'), ...
        mean(E_t,'omitnan'), ...
        max(E_t,[],'omitnan'));
    fprintf('chi_merge min / moyenne / max     : %.6f / %.6f / %.6f\n', ...
        min(chi_merge_t,[],'all'), ...
        mean(chi_merge_t,'all','omitnan'), ...
        max(chi_merge_t,[],'all'));
    fprintf('p_merge min / moyenne / max       : %.6f / %.6f / %.6f\n', ...
        min(p_merge_t,[],'all'), ...
        mean(p_merge_t,'all','omitnan'), ...
        max(p_merge_t,[],'all'));
end

function [N,p_link,E_t,beta0_geom,N1,N2,N3,edge_source] = ...
    local_params_beta0_geom(liens_file,n_target)

    if ~isfile(liens_file)
        error('Fichier introuvable : %s',liens_file);
    end

    data = load(liens_file);

    %% N
    if isfield(data,'N_mean_theory')
        N = double(data.N_mean_theory);
    elseif isfield(data,'N_all')
        N = mean(double(data.N_all(:)),'omitnan');
    elseif isfield(data,'N')
        N = double(data.N);
    elseif isfield(data,'lambda') && isfield(data,'surface_sphere')
        N = double(data.lambda)*double(data.surface_sphere);
    else
        error('Impossible de determiner N depuis %s.',liens_file);
    end

    %% p_link
    if isfield(data,'p_link_uniform')
        p_link = double(data.p_link_uniform);
    elseif isfield(data,'p_link')
        p_link = double(data.p_link);
    elseif isfield(data,'dmax') && isfield(data,'R')
        alpha_max = 2*asin(min(double(data.dmax)/(2*double(data.R)),1));
        p_link = (1-cos(alpha_max))/2;
    else
        error('Impossible de determiner p_link depuis %s.',liens_file);
    end

    p_link = min(max(p_link,0),1);

    %% E(t)
    if isfield(data,'mean_edges')
        E_raw = double(data.mean_edges(:));
        edge_source = 'mean_edges du fichier';

    elseif isfield(data,'num_edges_all')
        E_raw = mean(double(data.num_edges_all),1,'omitnan').';
        edge_source = 'moyenne de num_edges_all';

    elseif isfield(data,'E_delta')
        E_raw = double(data.E_delta);
        edge_source = 'E_delta du fichier';

    else
        E_raw = N*(N-1)/2*p_link;
        edge_source = 'N(N-1)p_link/2';
    end

    if isscalar(E_raw)
        E_t = repmat(E_raw,n_target,1);
    elseif numel(E_raw) == n_target
        E_t = E_raw;
    else
        x_raw = linspace(0,1,numel(E_raw));
        x_target = linspace(0,1,n_target);
        E_t = interp1(x_raw,E_raw,x_target,'linear','extrap').';
    end

    %% Approximation de beta0 par N1, N2 et N3
    c2_union = 1 + 3*sqrt(3)/(4*pi);
    c3_conn  = 1 + 3*sqrt(3)/(2*pi);
    c3_union = 1.80;

    N1 = N*max(1-p_link,0)^(N-1);

    if N >= 2
        N2 = local_nchoosek_real(N,2)*p_link ...
            * max(1-c2_union*p_link,0)^(N-2);
    else
        N2 = 0;
    end

    if N >= 3
        p_conn_3 = min(max(c3_conn*p_link^2,0),1);
        N3 = local_nchoosek_real(N,3)*p_conn_3 ...
            * max(1-c3_union*p_link,0)^(N-3);
    else
        N3 = 0;
    end

    beta0_geom = 1 + N1 + N2 + N3;
    beta0_geom = min(max(beta0_geom,1),N);
end

function c = local_nchoosek_real(n,k)
    if n < k
        c = 0;
    else
        c = exp(gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1));
    end
end
