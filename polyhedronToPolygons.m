function R = polyhedronToPolygons(vertices, faces, zCoords, sliceDir)
% R = polyhedronToPolygons(vertices, faces, zCoords)
%
% Converts a polyhedron to a set of polygons with
% different z-coordinates. This representation (the set of polygons) is
% used for defining the outline of organs/regions of interest in the DICOM
% file format (within the RTstruct in the DICOM file).
%
% Input parameters:
% vertices: n by 3 array of the vertices of the polyhedron. There are n
% vertices, and each vertex has an x, y and z-coordinate.
%
% faces: m by 3 array of faces of the polyhedron. Each face has three
% coordinates, and each row of this array points to the index of the three
% vertices in the "vertices" array. See isosurface for more details on the
% representation of a polyhedron as vertices/faces.
%
% zCoords: Array of k z-coordinates for which to compute polygons.
%
% sliceDir (optional): Selects which direction is the "z" direction along
% which the slices are found (1, 2 or 3). Default: 3.
%
%
% Output parameters:
% R: k times 1 cell array containing the polygons. Polygons are defined by
% q by 3 arrays of vertices (x, y and z coordinates; the z-coordinate is
% constant for each polygon). There is an edge between each adjacent vertex
% in the array, and an edge between the last and first vertex. There may be
% multiple polygons in each plane (even if the polyhedron is connected),
% therefore each cell of R is itself a cell array that contains zero or
% more polygon definitions. 
%


if nargin < 4
    sliceDir = 3;
end

vertices = vertices(:, [mod(sliceDir, 3) + 1, mod(sliceDir+1, 3) + 1, mod(sliceDir + 2, 3) + 1]);
    


R = cell(length(zCoords), 1);
Fmat = faceMatrix(vertices, faces);
for i = 1:length(zCoords)
    z = zCoords(i);
    indices = (min(Fmat(:,:,3), [], 2) < z) & (max(Fmat(:,:,3),[], 2) > z);
    f = faces(indices, :);
    R{i} = cell(0);
    count1 = 1;
    while ~isempty(f)        
        face = f(1, :);
        f = f(2:end, :);
        F = vertices(face', :);
        below = find(F(:, 3) < z, 1);
        above = find(F(:, 3) >= z, 1);
        next = find(any(f == face(below), 2) & any(f == face(above), 2));
        if size(next) > 1
            error('Invalid polyhedron');
        end
        count2 = 1;        
        R{i}{count1}(count2, :) = crossing(vertices(face(below), :), vertices(face(above), :), z);
        while ~isempty(next)
            below = find(f(next,:) == face(below));
            above = find(f(next,:) == face(above));            
            face = f(next, :);            
            f = f([1:(next-1) (next+1):end], :);            
            count2 = count2 + 1;
            R{i}{count1}(count2, :) = crossing(vertices(face(below), :), vertices(face(above), :), z);
            other = 6 - below - above;
            if vertices(face(other), 3) < z
                below = other;
            else
                above = other;
            end
            next = find(any(f == face(below), 2) & any(f == face(above), 2));
            if size(next) > 1
                error('Invalid polyhedron');
            end
        end
        count1 = count1 + 1;
    end
end

if sliceDir ~= 3
    i1 = mod(2*sliceDir, 3) + 1; i2 = mod(2*sliceDir+1, 3) + 1; i3 = mod(2*sliceDir + 2, 3) + 1;
    for i = 1:length(R)
        for j = 1:length(R{i})
            R{i}{j} = R{i}{j}(:, [i1, i2, i3]);
        end
    end
end

end

function f = faceMatrix(vertices, faces)
f = zeros(size(faces, 1), 3, 3);
f(:, 1, :) = vertices(faces(:, 1), :);
f(:, 2, :) = vertices(faces(:, 2), :);
f(:, 3, :) = vertices(faces(:, 3), :);
end

function c = crossing(a, b, z)
d = b - a;
u = (z - a(3))/d(3);
c = a + u * d;
end
