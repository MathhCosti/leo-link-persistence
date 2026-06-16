nTrials = 10000;   % nombre de paires testées
success = 0;

for k = 1:nTrials
    pair = randperm(N, 2);
    s = pair(1);
    t = pair(2);

    path = shortestpath(G, s, t);

    if ~isempty(path)
        success = success + 1;
    end
end

P_routing = success / nTrials;

fprintf("Probabilité de routage estimée = %.4f\n", P_routing);