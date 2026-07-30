function [t, vrel_mean, vrad_mean, ratio_mean] = plot_vrel_vs_vrad(pos_all, vel_all, t, dmax, varargin)
%PLOT_VREL_VS_VRAD Compare total relative speed and useful radial speed.
%
%   [t, vrel_mean, vrad_mean, ratio_mean] = plot_vrel_vs_vrad(pos_all, vel_all, t, dmax)
%
%   pos_all : positions satellites
%       Accepted formats:
%           - N x 3 x T
%           - T x N x 3
%           - cell array of length T, each cell N x 3
%
%   vel_all : velocities satellites, same format as pos_all.
%
%   t       : time vector, length T. If empty, t = 1:T.
%   dmax    : distance threshold. If empty or omitted, all pairs are used.
%
%   Optional name-value arguments:
%       'Pairs'       : 'all', 'linked', or 'near_boundary'
%                       all           : all pairs i<j
%                       linked        : only pairs with d_ij <= dmax
%                       near_boundary : pairs with dmax-width <= d_ij <= dmax+width
%                       default: 'linked' if dmax is provided, otherwise 'all'
%       'Width'       : width for near_boundary. Default: 0.05*dmax
%       'UseAbsRadial': true to plot mean(abs(v_rad)), false to plot mean(max(0,-v_rad)).
%                       default: true
%       'MakeFigure'  : true/false. Default: true
%
%   Outputs:
%       vrel_mean  : mean ||v_j-v_i|| over selected pairs at each time
%       vrad_mean  : mean useful radial component over selected pairs at each time
%       ratio_mean : vrad_mean ./ vrel_mean
%
%   Definitions:
%       r_ij = x_j - x_i
%       e_ij = r_ij / ||r_ij||
%       vrel_ij = v_j - v_i
%       v_rad = vrel_ij . e_ij
%
%   For merge, the useful component is max(0,-v_rad).
%   For a symmetric diagnostic, use abs(v_rad).

p = inputParser;
addParameter(p, 'Pairs', '', @(s)ischar(s) || isstring(s));
addParameter(p, 'Width', [], @(x)isempty(x) || (isscalar(x) && x >= 0));
addParameter(p, 'UseAbsRadial', true, @(x)islogical(x) || isnumeric(x));
addParameter(p, 'MakeFigure', true, @(x)islogical(x) || isnumeric(x));
parse(p, varargin{:});

pairs_mode = char(p.Results.Pairs);
width = p.Results.Width;
use_abs_radial = logical(p.Results.UseAbsRadial);
make_figure = logical(p.Results.MakeFigure);

if nargin < 4
    dmax = [];
end

if isempty(pairs_mode)
    if isempty(dmax)
        pairs_mode = 'all';
    else
        pairs_mode = 'linked';
    end
end

if strcmpi(pairs_mode, 'near_boundary') && isempty(dmax)
    error('dmax is required when Pairs = near_boundary.');
end

if isempty(width)
    if isempty(dmax)
        width = 0;
    else
        width = 0.05 * dmax;
    end
end

T = local_num_times(pos_all);
if nargin < 3 || isempty(t)
    t = 1:T;
end
t = t(:);

if numel(t) ~= T
    error('Length of t must match the number of time steps in pos_all/vel_all.');
end

vrel_mean = nan(T,1);
vrad_mean = nan(T,1);
ratio_mean = nan(T,1);
nb_pairs_used = zeros(T,1);

for k = 1:T
    X = local_get_slice(pos_all, k); % N x 3
    V = local_get_slice(vel_all, k); % N x 3

    N = size(X,1);
    if size(X,2) ~= 3 || size(V,2) ~= 3 || size(V,1) ~= N
        error('Each position/velocity slice must be N x 3.');
    end

    vals_vrel = [];
    vals_vrad = [];

    for i = 1:N-1
        rij = X(i+1:N,:) - X(i,:);
        dij = sqrt(sum(rij.^2,2));
        valid = dij > 0;

        switch lower(pairs_mode)
            case 'all'
                mask = valid;
            case 'linked'
                mask = valid & (dij <= dmax);
            case 'near_boundary'
                mask = valid & (dij >= dmax-width) & (dij <= dmax+width);
            otherwise
                error('Unknown Pairs mode. Use all, linked, or near_boundary.');
        end

        if ~any(mask)
            continue;
        end

        rij = rij(mask,:);
        dij = dij(mask);
        vij = V(i+1:N,:) - V(i,:);
        vij = vij(mask,:);

        eij = rij ./ dij;
        vrel = sqrt(sum(vij.^2,2));
        vrad_signed = sum(vij .* eij, 2);

        if use_abs_radial
            vrad = abs(vrad_signed);
        else
            % Useful for merge: only approaching pairs.
            vrad = max(0, -vrad_signed);
        end

        vals_vrel = [vals_vrel; vrel]; %#ok<AGROW>
        vals_vrad = [vals_vrad; vrad]; %#ok<AGROW>
    end

    nb_pairs_used(k) = numel(vals_vrel);
    if nb_pairs_used(k) > 0
        vrel_mean(k) = mean(vals_vrel, 'omitnan');
        vrad_mean(k) = mean(vals_vrad, 'omitnan');
        ratio_mean(k) = vrad_mean(k) / vrel_mean(k);
    end
end

if make_figure
    figure;
    plot(t, vrel_mean, 'LineWidth', 1.5); hold on;
    plot(t, vrad_mean, 'LineWidth', 1.5);
    grid on;
    xlabel('Temps (s)');
    ylabel('Vitesse moyenne');
    if use_abs_radial
        legend('||v_{rel}|| moyen', '|v_{rel} \cdot \hat r| moyen', 'Location', 'best');
        title('Comparaison vitesse relative totale et composante radiale');
    else
        legend('||v_{rel}|| moyen', '[-v_{rel} \cdot \hat r]_+ moyen', 'Location', 'best');
        title('Comparaison vitesse relative totale et composante radiale utile pour fusion');
    end

    figure;
    plot(t, ratio_mean, 'LineWidth', 1.5);
    grid on;
    xlabel('Temps (s)');
    ylabel('v_{rad} / ||v_{rel}||');
    title('Fraction radiale effective de la vitesse relative');
end

fprintf('\n=== Diagnostic vitesse radiale ===\n');
fprintf('Mode de paires : %s\n', pairs_mode);
fprintf('Nombre moyen de paires utilisees : %.2f\n', mean(nb_pairs_used));
fprintf('Moyenne ||v_rel||           : %.6g\n', mean(vrel_mean, 'omitnan'));
fprintf('Moyenne composante radiale  : %.6g\n', mean(vrad_mean, 'omitnan'));
fprintf('Ratio radial moyen          : %.6g\n', mean(ratio_mean, 'omitnan'));

end

function T = local_num_times(A)
    if iscell(A)
        T = numel(A);
    elseif ndims(A) == 3
        sz = size(A);
        if sz(2) == 3
            % N x 3 x T
            T = sz(3);
        elseif sz(3) == 3
            % T x N x 3
            T = sz(1);
        else
            error('3D arrays must be N x 3 x T or T x N x 3.');
        end
    else
        error('pos_all and vel_all must be cell arrays or 3D arrays.');
    end
end

function X = local_get_slice(A, k)
    if iscell(A)
        X = A{k};
    else
        sz = size(A);
        if sz(2) == 3
            X = A(:,:,k);       % N x 3 x T
        elseif sz(3) == 3
            X = squeeze(A(k,:,:)); % T x N x 3 -> N x 3
        else
            error('3D arrays must be N x 3 x T or T x N x 3.');
        end
    end
end
