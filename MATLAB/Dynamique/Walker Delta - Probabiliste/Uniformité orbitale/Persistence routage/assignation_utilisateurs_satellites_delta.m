clear; clc; close all;

%% ============================================================
%  ASSIGNATION UTILISATEURS - SATELLITES DANS UN MODELE DELTA
%
%  Principe :
%  - satellites sur orbites circulaires d'inclinaison commune ;
%  - utilisateurs fixes sur la Terre, qui tourne dans le repere inertiel ;
%  - un satellite est admissible si son elevation est superieure ou egale
%    a l'elevation minimale ;
%  - parmi les satellites admissibles, l'utilisateur choisit celui dont
%    la distance sol-satellite est minimale ;
%  - l'identifiant 0 represente une absence de satellite admissible.
%
%  Sorties :
%  - AssignedSatellite(k,u) : satellite assigne a l'utilisateur u au temps k ;
%  - AssignedDistance(k,u)  : distance correspondante [km] ;
%  - AssignedElevation(k,u) : elevation correspondante [deg] ;
%  - AssignmentEpisodes     : intervalles continus d'assignation ;
%  - durees moyenne, minimale et maximale des assignations.
%% ============================================================

rng(1); % Reproductibilite

%% 1. PARAMETRES PHYSIQUES
R_earth = 6371;                 % [km]
h = 550;                        % [km]
R_orbit = R_earth + h;          % [km]
mu = 398600;                    % [km^3/s^2]
omega_sat = sqrt(mu/R_orbit^3); % [rad/s]

T_sidereal = 86164;             % [s]
omega_earth = 2*pi/T_sidereal;  % [rad/s]

%% 2. GENERATION DES SATELLITES
lambda_sat = 4e-7;                    % [satellites/km^2]
surface_orbital_sphere = 4*pi*R_orbit^2;
N_sat = poissrnd(lambda_sat*surface_orbital_sphere);

inc_deg = 58;
inc = deg2rad(inc_deg);

Omega = 2*pi*rand(N_sat,1); % RAAN
u0 = 2*pi*rand(N_sat,1);    % phase initiale

%% 3. GENERATION DES UTILISATEURS
N_users = 2;

lat_fixed_deg = 58;
lat_fixed_rad = lat_fixed_deg*pi/180;
% Uniformite surfacique dans la bande directement survolee [-inc,+inc].
user_lat = [lat_fixed_rad ; lat_fixed_rad];
% asin((2*rand(N_users,1)-1)*sin(inc));
user_lon0 = -pi + 2*pi*rand(N_users,1);

%% 4. PARAMETRES D'ASSIGNATION ET TEMPORELS
elevation_min_deg = 20;             % elevation minimale [deg]
elevation_min = deg2rad(elevation_min_deg);

dt = 10;                            % pas temporel [s]
Tmax = 12000;                        % duree totale [s]
time_values = (0:dt:Tmax).';
Nt = numel(time_values);

%% 5. STOCKAGE DES RESULTATS
% 0 signifie qu'aucun satellite ne respecte le seuil d'elevation.
AssignedSatellite = zeros(Nt,N_users);
AssignedDistance = NaN(Nt,N_users);
AssignedElevation = NaN(Nt,N_users);
VisibleSatelliteCount = zeros(Nt,N_users);

%% 6. CALCUL DES ASSIGNATIONS
for k = 1:Nt
    t = time_values(k);

    % Positions des satellites dans le repere inertiel.
    u_t = mod(u0 + omega_sat*t,2*pi);
    sat_positions = walker_delta_positions(R_orbit,inc,Omega,u_t);

    % Positions des utilisateurs dans le repere inertiel.
    user_lon_t = mod(user_lon0 + omega_earth*t + pi,2*pi)-pi;
    user_positions = ground_user_positions(R_earth,user_lat,user_lon_t);

    for user_id = 1:N_users
        r_user = user_positions(user_id,:);

        % Vecteurs ligne de visee utilisateur -> satellites.
        rho = sat_positions - r_user;
        distances = sqrt(sum(rho.^2,2));

        % Zenith local de l'utilisateur.
        local_vertical = r_user/R_earth;

        % sin(elevation) = rho_hat . verticale_locale.
        sin_elevation = (rho*local_vertical.') ./ distances;
        sin_elevation = max(-1,min(1,sin_elevation));
        elevations = asin(sin_elevation);

        admissible = elevations >= elevation_min;
        VisibleSatelliteCount(k,user_id) = nnz(admissible);

        if any(admissible)
            candidate_ids = find(admissible);
            [best_distance,local_index] = min(distances(candidate_ids));
            best_satellite = candidate_ids(local_index);

            AssignedSatellite(k,user_id) = best_satellite;
            AssignedDistance(k,user_id) = best_distance;
            AssignedElevation(k,user_id) = rad2deg(elevations(best_satellite));
        end
    end
end

%% 7. CONSTRUCTION DES INTERVALLES CONTINUS D'ASSIGNATION
% Chaque ligne de AssignmentEpisodes correspond a un intervalle pendant
% lequel un utilisateur conserve le meme satellite.
episode_user = [];
episode_satellite = [];
episode_start = [];
episode_end = [];
episode_duration = [];
episode_mean_distance = [];
episode_min_elevation = [];

N_intervals = Nt-1; % intervalles [t_k,t_{k+1}[ couvrant exactement Tmax

for user_id = 1:N_users
    sequence = AssignedSatellite(1:N_intervals,user_id);

    start_index = 1;
    for k = 2:(N_intervals+1)
        sequence_changed = (k > N_intervals) || sequence(k) ~= sequence(start_index);

        if sequence_changed
            end_index = k-1;
            satellite_id = sequence(start_index);

            % On ne compte comme assignation que les intervalles sat > 0.
            if satellite_id > 0
                % L'intervalle represente [t_debut, t_fin + dt[.
                duration = (end_index-start_index+1)*dt;

                episode_user(end+1,1) = user_id;
                episode_satellite(end+1,1) = satellite_id;
                episode_start(end+1,1) = time_values(start_index);
                episode_end(end+1,1) = time_values(end_index+1);
                episode_duration(end+1,1) = duration;
                episode_mean_distance(end+1,1) = ...
                    mean(AssignedDistance(start_index:end_index,user_id), ...
                    'omitnan');
                episode_min_elevation(end+1,1) = ...
                    min(AssignedElevation(start_index:end_index,user_id), ...
                    [],'omitnan');
            end

            start_index = k;
        end
    end
end

AssignmentEpisodes = table( ...
    episode_user,episode_satellite,episode_start,episode_end, ...
    episode_duration,episode_mean_distance,episode_min_elevation, ...
    'VariableNames',{ ...
    'UserID','SatelliteID','StartTime_s','EndTime_s', ...
    'Duration_s','MeanDistance_km','MinElevation_deg'});

%% 8. STATISTIQUES GLOBALES DES DUREES D'ASSIGNATION
if isempty(episode_duration)
    mean_assignment_duration = NaN;
    min_assignment_duration = NaN;
    max_assignment_duration = NaN;
else
    mean_assignment_duration = mean(episode_duration);
    min_assignment_duration = min(episode_duration);
    max_assignment_duration = max(episode_duration);
end

%% 9. STATISTIQUES PAR UTILISATEUR
UserID = (1:N_users).';
MeanAssignmentDuration_s = NaN(N_users,1);
MinAssignmentDuration_s = NaN(N_users,1);
MaxAssignmentDuration_s = NaN(N_users,1);
NumberOfAssignments = zeros(N_users,1);
OutageDuration_s = zeros(N_users,1);
OutageFraction = zeros(N_users,1);
MeanVisibleSatellites = mean(VisibleSatelliteCount,1).';

for user_id = 1:N_users
    rows = AssignmentEpisodes.UserID == user_id;
    durations = AssignmentEpisodes.Duration_s(rows);

    NumberOfAssignments(user_id) = numel(durations);

    if ~isempty(durations)
        MeanAssignmentDuration_s(user_id) = mean(durations);
        MinAssignmentDuration_s(user_id) = min(durations);
        MaxAssignmentDuration_s(user_id) = max(durations);
    end

    outage_steps = nnz(AssignedSatellite(1:N_intervals,user_id)==0);
    OutageDuration_s(user_id) = outage_steps*dt;
    OutageFraction(user_id) = outage_steps/N_intervals;
end

UserStatistics = table( ...
    UserID,rad2deg(user_lat),rad2deg(user_lon0), ...
    NumberOfAssignments,MeanAssignmentDuration_s, ...
    MinAssignmentDuration_s,MaxAssignmentDuration_s, ...
    OutageDuration_s,OutageFraction,MeanVisibleSatellites, ...
    'VariableNames',{ ...
    'UserID','Latitude_deg','InitialLongitude_deg', ...
    'NumberOfAssignments','MeanDuration_s','MinDuration_s', ...
    'MaxDuration_s','OutageDuration_s','OutageFraction', ...
    'MeanVisibleSatellites'});

%% 10. AFFICHAGE TEXTUEL
fprintf('\n============================================================\n');
fprintf('ASSIGNATION UTILISATEURS - SATELLITES\n');
fprintf('============================================================\n');
fprintf('Nombre de satellites        : %d\n',N_sat);
fprintf('Nombre d''utilisateurs       : %d\n',N_users);
fprintf('Inclinaison                 : %.1f deg\n',inc_deg);
fprintf('Elevation minimale          : %.1f deg\n',elevation_min_deg);
fprintf('Pas temporel                : %.1f s\n',dt);
fprintf('Duree simulee               : %.1f s\n',Tmax);

fprintf('\nDurees des intervalles continus d''assignation :\n');
fprintf('  Duree moyenne             : %.2f s\n',mean_assignment_duration);
fprintf('  Duree minimale            : %.2f s\n',min_assignment_duration);
fprintf('  Duree maximale            : %.2f s\n',max_assignment_duration);

fprintf('\nSatellites successivement assignes a chaque utilisateur :\n');
for user_id = 1:N_users
    rows = AssignmentEpisodes.UserID == user_id;
    user_episodes = AssignmentEpisodes(rows,:);

    fprintf('\nUtilisateur %d | latitude = %.2f deg | longitude initiale = %.2f deg\n', ...
        user_id,rad2deg(user_lat(user_id)),rad2deg(user_lon0(user_id)));

    if isempty(user_episodes)
        fprintf('  Aucun satellite assigne pendant la simulation.\n');
    else
        for j = 1:height(user_episodes)
            fprintf(['  Satellite %3d : [%7.1f s, %7.1f s] | ' ...
                'duree = %6.1f s | distance moyenne = %7.1f km | ' ...
                'elevation min = %5.1f deg\n'], ...
                user_episodes.SatelliteID(j), ...
                user_episodes.StartTime_s(j), ...
                user_episodes.EndTime_s(j), ...
                user_episodes.Duration_s(j), ...
                user_episodes.MeanDistance_km(j), ...
                user_episodes.MinElevation_deg(j));
        end
    end
end

fprintf('\nStatistiques par utilisateur :\n');
disp(UserStatistics);

%% 11. SAUVEGARDE
save('assignation_utilisateurs_satellites_delta.mat', ...
    'R_earth','h','R_orbit','mu','omega_sat','omega_earth', ...
    'lambda_sat','N_sat','inc_deg','inc','Omega','u0', ...
    'N_users','user_lat','user_lon0', ...
    'elevation_min_deg','elevation_min','dt','Tmax','time_values', ...
    'AssignedSatellite','AssignedDistance','AssignedElevation', ...
    'VisibleSatelliteCount','AssignmentEpisodes','UserStatistics', ...
    'mean_assignment_duration','min_assignment_duration', ...
    'max_assignment_duration');

writetable(AssignmentEpisodes,'assignment_episodes.csv');
writetable(UserStatistics,'user_assignment_statistics.csv');

fprintf('\nResultats sauvegardes dans :\n');
fprintf('  assignation_utilisateurs_satellites_delta.mat\n');
fprintf('  assignment_episodes.csv\n');
fprintf('  user_assignment_statistics.csv\n');

%% ============================================================
%  FONCTIONS LOCALES
%% ============================================================
function positions = walker_delta_positions(R,inc,Omega,u)
    x = R.*(cos(Omega).*cos(u) ...
        - sin(Omega).*sin(u).*cos(inc));
    y = R.*(sin(Omega).*cos(u) ...
        + cos(Omega).*sin(u).*cos(inc));
    z = R.*sin(u).*sin(inc);
    positions = [x,y,z];
end

function positions = ground_user_positions(R_earth,latitude,longitude)
    x = R_earth.*cos(latitude).*cos(longitude);
    y = R_earth.*cos(latitude).*sin(longitude);
    z = R_earth.*sin(latitude);
    positions = [x,y,z];
end
