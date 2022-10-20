function plotPolygons(R, color)
% plotPolygons(R, color)
%
% Plots a set of polygons created by polyhedronToPolygon in 3D 
%
%Parameters: 
% R: Data structure returned from polyhedronToPolygon
%
% color (optional): Color to plot, given in any format supported by 
% matlab's plot function (e.g. 'b', '#00ff00', [1 0 1])

if nargin < 2
    color = 'k';
end
Q = zeros(0, 3);
for i = 1:length(R)
    for j = 1:length(R{i})
        Q = [Q;R{i}{j};R{i}{j}(1, :);[nan, nan, nan]];
    end
end
plot3(Q(:,1), Q(:,2), Q(:,3), 'Color', color);

axis equal