function [time_values, chi_t, E_t, beta0_iso] = plot_chi_temp_liens(liens_file)
% plot_chi_temp_liens
% Trace l'evolution temporelle du facteur correctif topologique :
%
%   chi(t) = (N - beta0_iso) / |E(t)|
%
% ou |E(t)| est le nombre moyen de liens temporel lu dans liens_inter.mat.
% beta0_iso est estime avec la theorie des satellites isoles :
%
%   beta0_iso = 1 + N (1 - p_link)^(N-1)
%
% Utilisation :
%   plot_chi_temp_liens
%   plot_chi_temp_liens('liens_inter.mat')
%   [t, chi_t, E_t, beta0_iso] = plot_chi_temp_liens('liens_inter.mat');

script_dir = fileparts(mfilename('fullpath'));
liens_file = fullfile(script_dir, '..', 'Nombre liens', 'liens_inter.mat');

    data = load(liens_file);

    %% Nombre de liens temporel |E(t)|
    if isfield(data, 'mean_edges')
        E_t = double(data.mean_edges(:));
    elseif isfield(data, 'num_edges_all')
        E_t = mean(double(data.num_edges_all), 1, 'omitnan').';
    else
        error('Le fichier %s doit contenir mean_edges ou num_edges_all.', liens_file);
    end

    %% Axe temporel
    if isfield(data, 'time_values')
        time_values = double(data.time_values(:));
    elseif isfield(data, 'dt')
        dt = double(data.dt);
        time_values = (0:numel(E_t)-1).' * dt;
    else
        time_values = (0:numel(E_t)-1).';
    end

    % Si besoin, on ajuste la longueur de l'axe temporel.
    if numel(time_values) ~= numel(E_t)
        time_values = linspace(time_values(1), time_values(end), numel(E_t)).';
    end

    %% Nombre moyen de satellites N
    if isfield(data, 'N_mean_theory')
        N = double(data.N_mean_theory);
    elseif isfield(data, 'N_all')
        N = mean(double(data.N_all(:)), 'omitnan');
    elseif isfield(data, 'lambda') && isfield(data, 'surface_sphere')
        N = double(data.lambda) * double(data.surface_sphere);
    else
        error('Impossible de determiner N depuis %s.', liens_file);
    end

    %% Probabilite de lien uniforme p_link
    if isfield(data, 'p_link_uniform')
        p_link = double(data.p_link_uniform);
    elseif isfield(data, 'dmax') && isfield(data, 'R')
        alpha_max = 2 * asin(double(data.dmax) / (2 * double(data.R)));
        p_link = (1 - cos(alpha_max)) / 2;
    else
        error('Impossible de determiner p_link depuis %s.', liens_file);
    end

    %% Approximation beta0 isoles
    beta0_iso = 1 + N * (1 - p_link)^(N - 1);

    %% Facteur correctif temporel
    E_t = max(E_t, eps);
    chi_t = (N - beta0_iso) ./ E_t;
    chi_t = min(max(chi_t, 0), 1);

    %% Affichages utiles
    fprintf('N moyen utilise                 : %.3f\n', N);
    fprintf('p_link uniforme utilise          : %.6f\n', p_link);
    fprintf('beta0_iso theorique              : %.3f\n', beta0_iso);
    fprintf('Liens moyens |E(t)| : min/mean/max = %.3f / %.3f / %.3f\n', ...
        min(E_t), mean(E_t, 'omitnan'), max(E_t));
    fprintf('chi(t) : min/mean/max = %.4f / %.4f / %.4f\n', ...
        min(chi_t), mean(chi_t, 'omitnan'), max(chi_t));

    %% Graphe chi(t)
    figure;
    plot(time_values, chi_t, 'LineWidth', 1.6);
    grid on;
    xlabel('Temps (s)');
    ylabel('\chi(t) = (N - \beta_0^{iso}) / |E(t)|');
    title('Facteur correctif topologique \chi(t)');

    %% Graphe compare avec |E(t)|
    figure;
    yyaxis left;
    plot(time_values, chi_t, 'LineWidth', 1.6);
    ylabel('\chi(t)');

    yyaxis right;
    plot(time_values, E_t, '--', 'LineWidth', 1.3);
    ylabel('|E(t)|');

    grid on;
    xlabel('Temps (s)');
    title('Evolution de \chi(t) et du nombre de liens');
    legend('\chi(t)', '|E(t)|', 'Location', 'best');
end
