clear; clc; close all;

%% ============================================================
% DUREE MOYENNE D'UNE ROUTE SOL-SATELLITES-SOL
%
% Le modele utilise :
%
%   beta_route = beta_assign,A + beta_assign,B
%                + H_mean * beta_break_ISL
%
% avec
%
%   beta_assign,A = 1/T_assign,A
%   beta_assign,B = 1/T_assign,B
%
% et, si q_break_link est une probabilite de rupture sur Delta_t :
%
%   beta_break_ISL = -log(1-q_break_link)/Delta_t.
%
% Finalement :
%
%   T_route_mean = 1/beta_route.
%
% Hypotheses :
%   - H est remplace par sa valeur moyenne ;
%   - les ruptures des liens sont independantes ;
%   - la route est fixe : la rupture d'un seul ISL interrompt la route ;
%   - les durees d'assignation sont modelisees par des lois exponentielles.
%% ============================================================

%% Choix utilisateur
user_lat_A_deg = 0;
user_lat_B_deg = 45;

% 'corrected', 'empirical' ou 'old'
assignment_source = 'corrected';

% 'empirical' ou 'theoretical'
hop_source = 'theoretical';

%% Fichiers de resultats
script_dir = fileparts(mfilename('fullpath'));

assignment_file = fullfile(script_dir,'verification_duree_assignation_corrigee.mat');

shortest_path_file = fullfile(script_dir,'distribution_shortest_path_delta.mat');

pbreak_file = fullfile(script_dir,'..', 'Probabilité disparition', 'pbreak_theorique_walker_delta_results.mat');

assert(isfile(assignment_file), ...
    'Fichier d''assignation introuvable : %s',assignment_file);
assert(isfile(shortest_path_file), ...
    'Fichier shortest path introuvable : %s',shortest_path_file);
assert(isfile(pbreak_file), ...
    'Fichier p_break introuvable : %s',pbreak_file);

assignment_data = load(assignment_file);
path_data = load(shortest_path_file);
break_data = load(pbreak_file);

%% ============================================================
% 1. Nombre moyen de sauts
%% ============================================================

switch lower(hop_source)
    case 'empirical'
        assert(isfield(path_data,'mean_H'), ...
            'Le fichier shortest path ne contient pas mean_H.');
        H_mean = double(path_data.mean_H);
        H_label = 'empirique';

    case 'theoretical'
        assert(isfield(path_data,'mean_H_theory'), ...
            ['Le fichier shortest path ne contient pas ' ...
             'mean_H_theory.']);
        H_mean = double(path_data.mean_H_theory);
        H_label = 'theorique Delta';

    otherwise
        error('hop_source doit valoir empirical ou theoretical.');
end

validateattributes(H_mean,{'numeric'}, ...
    {'scalar','real','finite','positive'});

%% ============================================================
% 2. Durees moyennes d'assignation
%% ============================================================

required_assignment_fields = {'lambda_values','user_lat_deg'};
for k = 1:numel(required_assignment_fields)
    assert(isfield(assignment_data,required_assignment_fields{k}), ...
        'Champ manquant dans le fichier d''assignation : %s', ...
        required_assignment_fields{k});
end

lambda_values = double(assignment_data.lambda_values(:));
user_lat_deg = double(assignment_data.user_lat_deg(:));

switch lower(assignment_source)
    case 'corrected'
        field_assign = 'MeanAssign_corrected';
        assignment_label = 'heuristique corrigee';

    case 'empirical'
        field_assign = 'MeanAssign_emp';
        assignment_label = 'empirique';

    case 'old'
        field_assign = 'MeanAssign_old';
        assignment_label = 'ancienne heuristique';

    otherwise
        error(['assignment_source doit valoir corrected, ' ...
               'empirical ou old.']);
end

assert(isfield(assignment_data,field_assign), ...
    'Le fichier d''assignation ne contient pas %s.',field_assign);

Tassign_matrix = double(assignment_data.(field_assign));

% La densite utilisee dans le code shortest path est privilegiee.
if isfield(path_data,'lambda')
    lambda_target = double(path_data.lambda);
else
    warning(['Le fichier shortest path ne contient pas lambda. ' ...
             'La premiere densite du fichier d''assignation est utilisee.']);
    lambda_target = lambda_values(1);
end

[lambda_error,il] = min(abs(lambda_values-lambda_target));

if lambda_error > max(1e-12,1e-6*abs(lambda_target))
    warning(['Aucune densite d''assignation exactement egale a ' ...
             'lambda=%.4e. La valeur la plus proche %.4e est utilisee.'], ...
             lambda_target,lambda_values(il));
end

[lat_error_A,qA] = min(abs(user_lat_deg-user_lat_A_deg));
[lat_error_B,qB] = min(abs(user_lat_deg-user_lat_B_deg));

if lat_error_A > 1e-9
    warning('Latitude A %.2f deg absente : %.2f deg utilisee.', ...
        user_lat_A_deg,user_lat_deg(qA));
end

if lat_error_B > 1e-9
    warning('Latitude B %.2f deg absente : %.2f deg utilisee.', ...
        user_lat_B_deg,user_lat_deg(qB));
end

Tassign_A = Tassign_matrix(il,qA);
Tassign_B = Tassign_matrix(il,qB);

validateattributes(Tassign_A,{'numeric'}, ...
    {'scalar','real','finite','positive'});
validateattributes(Tassign_B,{'numeric'}, ...
    {'scalar','real','finite','positive'});

beta_assign_A = 1/Tassign_A;
beta_assign_B = 1/Tassign_B;

%% ============================================================
% 3. Taux de rupture d'un lien ISL
%% ============================================================

assert(isfield(break_data,'Delta_t'), ...
    'Le fichier p_break ne contient pas Delta_t.');

Delta_t_break = double(break_data.Delta_t);

% Pour une route fixe, toute rupture d'un lien utilise par la route
% provoque sa disparition. On utilise donc q_break_link, et non
% p_break_delta qui correspond a une rupture au niveau d'une composante.
if isfield(break_data,'q_break_link')
    q_break_ISL = double(break_data.q_break_link);
    qbreak_label = 'q_break_link';
elseif isfield(break_data,'q_break_link_raw')
    q_break_ISL = double(break_data.q_break_link_raw);
    qbreak_label = 'q_break_link_raw';
else
    error(['Le fichier p_break ne contient ni q_break_link ' ...
           'ni q_break_link_raw.']);
end

validateattributes(q_break_ISL,{'numeric'}, ...
    {'scalar','real','finite','>=',0,'<',1});
validateattributes(Delta_t_break,{'numeric'}, ...
    {'scalar','real','finite','positive'});

% Conversion exacte probabilite discrete -> taux exponentiel.
beta_break_ISL = -log1p(-q_break_ISL)/Delta_t_break;

% Approximation lineaire, utile lorsque q_break_ISL << 1.
beta_break_ISL_linear = q_break_ISL/Delta_t_break;

%% ============================================================
% 4. Duree moyenne de la route
%% ============================================================

beta_route_GSL = beta_assign_A+beta_assign_B;
beta_route_ISL = H_mean*beta_break_ISL;
beta_route_total = beta_route_GSL+beta_route_ISL;

Troute_mean = 1/beta_route_total;

% Version lineaire pour verifier que q_break est suffisamment petit.
beta_route_total_linear = ...
    beta_route_GSL+H_mean*beta_break_ISL_linear;
Troute_mean_linear = 1/beta_route_total_linear;

% Probabilite de survie de la route pendant un pas Delta_t_break.
p_survival_GSL_one_step = ...
    exp(-beta_route_GSL*Delta_t_break);

p_survival_ISL_one_step = ...
    (1-q_break_ISL)^H_mean;

p_survival_route_one_step = ...
    p_survival_GSL_one_step*p_survival_ISL_one_step;

p_break_route_one_step = 1-p_survival_route_one_step;

%% ============================================================
% 5. Resultats
%% ============================================================

fprintf('\n============================================================\n');
fprintf('DUREE MOYENNE DE LA ROUTE\n');
fprintf('============================================================\n');
fprintf('Densite utilisee                         : %.4e km^-2\n', ...
    lambda_values(il));
fprintf('Latitude utilisateur A                   : %.1f deg\n', ...
    user_lat_deg(qA));
fprintf('Latitude utilisateur B                   : %.1f deg\n', ...
    user_lat_deg(qB));
fprintf('Source des durees d''assignation          : %s\n', ...
    assignment_label);
fprintf('T_assign,A                               : %.6f s\n', ...
    Tassign_A);
fprintf('T_assign,B                               : %.6f s\n', ...
    Tassign_B);
fprintf('beta_assign,A                            : %.10e s^-1\n', ...
    beta_assign_A);
fprintf('beta_assign,B                            : %.10e s^-1\n', ...
    beta_assign_B);
fprintf('Source de H moyen                        : %s\n',H_label);
fprintf('H moyen                                  : %.6f\n',H_mean);
fprintf('Probabilite de rupture ISL (%s) : %.10f\n', ...
    qbreak_label,q_break_ISL);
fprintf('Pas temporel de p_break                  : %.6f s\n', ...
    Delta_t_break);
fprintf('beta_break_ISL exact                     : %.10e s^-1\n', ...
    beta_break_ISL);
fprintf('Contribution GSL au taux de route        : %.10e s^-1\n', ...
    beta_route_GSL);
fprintf('Contribution ISL au taux de route        : %.10e s^-1\n', ...
    beta_route_ISL);
fprintf('Taux total de disparition de la route    : %.10e s^-1\n', ...
    beta_route_total);
fprintf('Probabilite de rupture sur un pas        : %.10f\n', ...
    p_break_route_one_step);
fprintf('T_route moyen                            : %.6f s\n', ...
    Troute_mean);
fprintf('T_route moyen, approximation lineaire    : %.6f s\n', ...
    Troute_mean_linear);
fprintf('============================================================\n');

Results = table( ...
    lambda_values(il),user_lat_deg(qA),user_lat_deg(qB), ...
    Tassign_A,Tassign_B,H_mean,q_break_ISL,Delta_t_break, ...
    beta_assign_A,beta_assign_B,beta_break_ISL, ...
    beta_route_GSL,beta_route_ISL,beta_route_total, ...
    p_break_route_one_step,Troute_mean,Troute_mean_linear, ...
    'VariableNames',{ ...
    'Lambda','LatitudeA_deg','LatitudeB_deg', ...
    'TassignA_s','TassignB_s','MeanHopCount', ...
    'ISLBreakProbabilityPerStep','BreakTimeStep_s', ...
    'BetaAssignA_per_s','BetaAssignB_per_s', ...
    'BetaBreakISL_per_s','BetaRouteGSL_per_s', ...
    'BetaRouteISL_per_s','BetaRouteTotal_per_s', ...
    'RouteBreakProbabilityPerStep', ...
    'MeanRouteLifetime_s','MeanRouteLifetimeLinear_s'});

disp(Results);

%% ============================================================
% 6. Contributions au taux de disparition
%% ============================================================

figure;
bar(categorical({'Assignation A','Assignation B','ISL'}), ...
    [beta_assign_A,beta_assign_B,beta_route_ISL]);
grid on;
ylabel('Contribution au taux de disparition (s^{-1})');
title(sprintf('Decomposition du taux de route, T_{route}=%.2f s', ...
    Troute_mean));

%% ============================================================
% 7. Sauvegarde
%% ============================================================

writetable(Results,'duree_moyenne_route_results.csv');

save('duree_moyenne_route_results.mat', ...
    'Results','Troute_mean','Troute_mean_linear', ...
    'H_mean','Tassign_A','Tassign_B', ...
    'q_break_ISL','Delta_t_break', ...
    'beta_assign_A','beta_assign_B','beta_break_ISL', ...
    'beta_route_GSL','beta_route_ISL','beta_route_total', ...
    'p_break_route_one_step', ...
    'assignment_source','hop_source');

fprintf('\nResultats sauvegardes dans :\n');
fprintf('  duree_moyenne_route_results.csv\n');
fprintf('  duree_moyenne_route_results.mat\n');
