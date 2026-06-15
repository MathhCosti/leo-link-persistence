% Tirer deux satellites distincts au hasard
pair = randperm(N, 2);
s = pair(1);
t = pair(2);

% Vérifier s'il existe un chemin de routage entre s et t
path = shortestpath(G, s, t);
nb_hops = length(path) - 1;

if isempty(path)
    disp("Pas de chemin de routage entre les deux satellites");
else
    disp("Chemin de routage trouvé :");
    disp(nb_hops);
    disp(path);
end