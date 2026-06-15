dmax = 1500; % km, exemple

positions = [x y z];
D = squareform(pdist(positions));

A = (D <= dmax) & (D > 0); % matrice d'adjacence
G = graph(A);

figure;
hold on;

% Optionnel : afficher la sphère orbitale si R est déjà défini
[Xs, Ys, Zs] = sphere(80);
surf(R*Xs, R*Ys, R*Zs, ...
    'FaceAlpha', 0.05, ...
    'EdgeColor', 'none');

% Affichage des liens intersatellites
for i = 1:size(A,1)
    for j = i+1:size(A,2)
        if A(i,j)
            plot3([x(i), x(j)], ...
                  [y(i), y(j)], ...
                  [z(i), z(j)], ...
                  'k-', 'LineWidth', 0.7);
        end
    end
end

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