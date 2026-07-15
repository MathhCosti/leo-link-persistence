dmax = 1500; % km
positions = [x y z];
D = squareform(pdist(positions));

A = (D <= dmax) & (D > 0); % matrice d'adjacence
G = graph(A);

figure;
hold on;

% Affichage des liens intersatellites en un seul plot3
[row, col] = find(triu(A, 1));   % arêtes i<j
E = length(row);

Xlinks = NaN(3*E, 1);
Ylinks = NaN(3*E, 1);
Zlinks = NaN(3*E, 1);

Xlinks(1:3:end) = x(row);
Xlinks(2:3:end) = x(col);

Ylinks(1:3:end) = y(row);
Ylinks(2:3:end) = y(col);

Zlinks(1:3:end) = z(row);
Zlinks(2:3:end) = z(col);

plot3(Xlinks, Ylinks, Zlinks, 'k-', 'LineWidth', 0.5);

% Affichage des satellites
scatter3(x, y, z, 35, 'filled');

axis equal;
grid on;
xlabel('x (km)');
ylabel('y (km)');
zlabel('z (km)');
title('Graphe sphérique des liens intersatellites');

view(3);
rotate3d on;
hold off;