dmax = 1500; % km, exemple

positions = [x y z];
D = squareform(pdist(positions));

A = (D <= dmax) & (D > 0); % matrice d'adjacence
G = graph(A);

figure;
plot(G);
title('Graphe des liens intersatellites');