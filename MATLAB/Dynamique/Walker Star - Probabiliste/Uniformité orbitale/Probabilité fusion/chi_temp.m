clear; clc; close all;

%% ============================================================
%  LIENS-PONTS — WALKER-DELTA A UNIFORMITE ORBITALE
%
%  Entrées :
%    - plink_results.mat
%    - betti_results.mat
%
%  Le facteur correctif est défini par :
%
%      chi_bridge = (N - beta0) / E
%
%  où :
%    - N - beta0 est le nombre d'arêtes d'une forêt couvrante,
%      donc le nombre minimal de liens indispensables pour relier
%      les composantes ;
%    - E est le nombre moyen total de liens.
%
%  Pour le modèle Delta orbital :
%
%      E_delta = N(N-1)/2 * p_link_delta
%
%  Le fichier betti_results.mat fournit :
%    - beta0_emp    : réalisations empiriques de beta0 ;
%    - beta0_th_123 : approximation théorique jusqu'aux trimères.
%
%  Attention :
%  betti_results.mat ne contient pas le nombre empirique d'arêtes
%  de chaque réalisation. La courbe empirique ci-dessous utilise
%  donc E_delta théorique comme dénominateur commun.
%% ============================================================

%% 1. Fichiers d'entrée
script_dir = fileparts(mfilename('fullpath'));

plink_file = fullfile(script_dir, '..', 'plink_results.mat');
betti_file = fullfile(script_dir, '..', 'betti_results.mat');

if ~isfile(plink_file)
    error('Fichier introuvable : %s', plink_file);
end

if ~isfile(betti_file)
    error('Fichier introuvable : %s', betti_file);
end

Splink = load(plink_file);
Sbetti = load(betti_file);

%% 2. Vérification des variables nécessaires
required_plink = {'N', 'p_link_delta'};
for k = 1:numel(required_plink)
    if ~isfield(Splink, required_plink{k})
        error('Variable ''%s'' absente de %s.', ...
            required_plink{k}, plink_file);
    end
end

required_betti = {'N', 'beta0_emp', 'beta0_th_123'};
for k = 1:numel(required_betti)
    if ~isfield(Sbetti, required_betti{k})
        error('Variable ''%s'' absente de %s.', ...
            required_betti{k}, betti_file);
    end
end

%% 3. Paramètres communs
N_plink = double(Splink.N);
N_betti = double(Sbetti.N);

if N_plink ~= N_betti
    error(['Incohérence entre les fichiers :\n' ...
        'N = %g dans plink_results.mat\n' ...
        'N = %g dans betti_results.mat'], ...
        N_plink, N_betti);
end

N = N_plink;
p_link_delta = double(Splink.p_link_delta);

% Paramètres facultatifs
inc_deg = get_optional_scalar(Splink, Sbetti, 'inc_deg', NaN);
R = get_optional_scalar(Splink, Sbetti, 'R', NaN);
dmax = get_optional_scalar(Splink, Sbetti, 'dmax', NaN);
alpha_max = get_optional_scalar(Splink, Sbetti, 'alpha_max', NaN);

%% 4. Nombre moyen théorique de liens Delta
E_delta_from_plink = N*(N-1)/2 * p_link_delta;

if isfield(Splink, 'E_delta')
    E_delta_saved = double(Splink.E_delta);

    relative_E_error = abs(E_delta_saved - E_delta_from_plink) ...
        / max(abs(E_delta_from_plink), eps);

    if relative_E_error > 1e-10
        warning(['E_delta sauvegardé et E_delta recalculé diffèrent.\n' ...
            'Valeur sauvegardée : %.12g\n' ...
            'Valeur recalculée  : %.12g'], ...
            E_delta_saved, E_delta_from_plink);
    end

    % On utilise la valeur directement issue du p_link chargé.
    E_delta = E_delta_from_plink;
else
    E_delta_saved = NaN;
    relative_E_error = NaN;
    E_delta = E_delta_from_plink;
end

if E_delta <= 0
    error('Le nombre moyen de liens E_delta doit être strictement positif.');
end

%% 5. Bettis empiriques et théoriques
beta0_emp = double(Sbetti.beta0_emp(:));
beta0_th_123 = double(Sbetti.beta0_th_123);

beta0_emp_mean = mean(beta0_emp);
beta0_emp_std = std(beta0_emp);
beta0_emp_median = median(beta0_emp);

% Valeurs facultatives des contributions théoriques
N1_th_local = get_field_or_nan(Sbetti, 'N1_th_local');
N2_th_geom = get_field_or_nan(Sbetti, 'N2_th_geom');
N3_th_geom = get_field_or_nan(Sbetti, 'N3_th_geom');
C_macro = get_field_or_nan(Sbetti, 'C_macro');

%% 6. Nombre de liens-ponts équivalents
%
% Pour un graphe possédant beta0 composantes, une forêt couvrante
% contient exactement N-beta0 arêtes.
%
% Cette quantité ne donne pas nécessairement le nombre exact de ponts
% du graphe lorsque les composantes contiennent des cycles. Elle constitue
% le nombre minimal de liens topologiquement indispensables.

n_bridge_equiv_emp = N - beta0_emp;
n_bridge_equiv_emp_mean = mean(n_bridge_equiv_emp);
n_bridge_equiv_emp_std = std(n_bridge_equiv_emp);

n_bridge_equiv_th = N - beta0_th_123;

%% 7. Facteur correctif chi_bridge
chi_bridge_emp = n_bridge_equiv_emp / E_delta;
chi_bridge_emp_mean = mean(chi_bridge_emp);
chi_bridge_emp_std = std(chi_bridge_emp);
chi_bridge_emp_median = median(chi_bridge_emp);

chi_bridge_th_123 = n_bridge_equiv_th / E_delta;

% Bornage numérique naturel
chi_bridge_emp = min(max(chi_bridge_emp, 0), 1);
chi_bridge_emp_mean = min(max(chi_bridge_emp_mean, 0), 1);
chi_bridge_emp_median = min(max(chi_bridge_emp_median, 0), 1);
chi_bridge_th_123 = min(max(chi_bridge_th_123, 0), 1);

%% 8. Comparaison des résultats
fprintf('\n=== Liens-ponts — Walker-Delta orbital ===\n');
fprintf('N                                      : %d\n', N);
fprintf('p_link_delta                           : %.12g\n', p_link_delta);
fprintf('E_delta = N(N-1)p_link/2               : %.12g\n', E_delta);

if isfinite(E_delta_saved)
    fprintf('E_delta sauvegardé                     : %.12g\n', ...
        E_delta_saved);
end

fprintf('\n--- Betti 0 ---\n');
fprintf('beta0 empirique moyen                  : %.8f\n', ...
    beta0_emp_mean);
fprintf('beta0 empirique écart-type             : %.8f\n', ...
    beta0_emp_std);
fprintf('beta0 théorique jusqu''aux trimères     : %.8f\n', ...
    beta0_th_123);

fprintf('\n--- Liens-ponts équivalents N-beta0 ---\n');
fprintf('Valeur empirique moyenne               : %.8f\n', ...
    n_bridge_equiv_emp_mean);
fprintf('Valeur théorique jusqu''aux trimères    : %.8f\n', ...
    n_bridge_equiv_th);

fprintf('\n--- Facteur correctif chi_bridge ---\n');
fprintf('chi_bridge empirique moyen             : %.10f\n', ...
    chi_bridge_emp_mean);
fprintf('chi_bridge empirique écart-type        : %.10f\n', ...
    chi_bridge_emp_std);
fprintf('chi_bridge théorique jusqu''aux trimères: %.10f\n', ...
    chi_bridge_th_123);

%% 9. Figure : distribution empirique de beta0
figure;
histogram(beta0_emp, 'Normalization', 'probability');
hold on;
grid on;

xline(beta0_emp_mean, '--', ...
    sprintf('Moyenne empirique = %.3f', beta0_emp_mean), ...
    'LineWidth', 1.6, ...
    'LabelHorizontalAlignment', 'left');

xline(beta0_th_123, '-.', ...
    sprintf('Théorie 1+2+3 = %.3f', beta0_th_123), ...
    'LineWidth', 1.8, ...
    'LabelHorizontalAlignment', 'right');

xlabel('\beta_0');
ylabel('Fréquence');
title('Distribution empirique de \beta_0 — Walker-Delta orbital');
legend('Réalisations empiriques', ...
    'Moyenne empirique', ...
    'Théorie jusqu''aux trimères', ...
    'Location', 'best');
hold off;

%% 10. Figure : distribution empirique de chi_bridge
figure;
histogram(chi_bridge_emp, 'Normalization', 'probability');
hold on;
grid on;

xline(chi_bridge_emp_mean, '--', ...
    sprintf('Moyenne empirique = %.4f', chi_bridge_emp_mean), ...
    'LineWidth', 1.6, ...
    'LabelHorizontalAlignment', 'left');

xline(chi_bridge_th_123, '-.', ...
    sprintf('Théorie 1+2+3 = %.4f', chi_bridge_th_123), ...
    'LineWidth', 1.8, ...
    'LabelHorizontalAlignment', 'right');

xlabel('\chi_{bridge} = (N-\beta_0)/E_\Delta');
ylabel('Fréquence');
title('Facteur correctif des liens-ponts — Walker-Delta orbital');
legend('Réalisations empiriques', ...
    'Moyenne empirique', ...
    'Théorie jusqu''aux trimères', ...
    'Location', 'best');
hold off;

%% 11. Figure : comparaison synthétique
figure;

values_beta0 = [beta0_emp_mean, beta0_th_123];
bar(values_beta0);
grid on;
set(gca, 'XTick', 1:2, ...
    'XTickLabel', {'Empirique', 'Théorie 1+2+3'});
ylabel('\beta_0');
title('Comparaison de \beta_0 — Walker-Delta orbital');

figure;

values_chi = [chi_bridge_emp_mean, chi_bridge_th_123];
bar(values_chi);
grid on;
set(gca, 'XTick', 1:2, ...
    'XTickLabel', {'Empirique', 'Théorie 1+2+3'});
ylabel('\chi_{bridge}');
ylim([0, max(0.05, 1.15*max(values_chi))]);
title('Comparaison du facteur correctif — Walker-Delta orbital');

%% 12. Profil local de p_link(phi), s'il est disponible
if isfield(Sbetti, 'phi_vals') && isfield(Sbetti, 'p_link_phi')
    phi_vals = double(Sbetti.phi_vals(:));
    p_link_phi = double(Sbetti.p_link_phi(:));

    if numel(phi_vals) == numel(p_link_phi)
        figure;
        plot(phi_vals*180/pi, p_link_phi, 'LineWidth', 1.8);
        hold on;
        grid on;

        yline(p_link_delta, '--', ...
            sprintf('p_{link}^{Delta} = %.4g', p_link_delta), ...
            'LineWidth', 1.5);

        xlabel('Latitude orbitale \phi [deg]');
        ylabel('p_{link}(\phi)');
        title('Probabilité locale de lien — Walker-Delta orbital');
        legend('p_{link}(\phi)', ...
            'Moyenne orbitale chargée', ...
            'Location', 'best');
        hold off;
    end
end

%% 13. Tableau récapitulatif
T = table( ...
    N, p_link_delta, E_delta, ...
    beta0_emp_mean, beta0_emp_std, beta0_th_123, ...
    n_bridge_equiv_emp_mean, n_bridge_equiv_th, ...
    chi_bridge_emp_mean, chi_bridge_emp_std, chi_bridge_th_123, ...
    'VariableNames', { ...
    'N', 'p_link_delta', 'E_delta', ...
    'beta0_emp_mean', 'beta0_emp_std', 'beta0_th_123', ...
    'n_bridge_emp_mean', 'n_bridge_th_123', ...
    'chi_bridge_emp_mean', 'chi_bridge_emp_std', ...
    'chi_bridge_th_123'});

disp(T);

%% ============================================================
%  FONCTIONS LOCALES
%% ============================================================

function value = get_field_or_nan(S, field_name)
    if isfield(S, field_name)
        value = double(S.(field_name));
    else
        value = NaN;
    end
end

function value = get_optional_scalar(S1, S2, field_name, default_value)
    if isfield(S1, field_name)
        value = double(S1.(field_name));
    elseif isfield(S2, field_name)
        value = double(S2.(field_name));
    else
        value = default_value;
    end
end
