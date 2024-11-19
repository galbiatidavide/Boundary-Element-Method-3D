TR = data.TR;
triangles = TR.triangles;

j = 2;
% get goordinates of triangle j
triangle = squeeze(triangles(j,:,:));
% vertexes
x1 = triangle(1,1);
x2 = triangle(2,1);
x3 = triangle(3,1);
y1 = triangle(1,2);
y2 = triangle(2,2);
y3 = triangle(3,2);
z1 = triangle(1,3);
z2 = triangle(2,3);
z3 = triangle(3,3);

% compute the area of the triangle (det(Jacobian))
area = triangle_area_3D([x1 y1 z1], [x2 y2 z2], [x3 y3 z3]);

% get the iglo
iglo = TR.ConnectivityList(j,:)
single_layer_potential_singular_quadratic_vertices(area, triangle, 1)
single_layer_potential_singular_quadratic_vertices(area, triangle, 2)
single_layer_potential_singular_quadratic_vertices(area, triangle, 3)
single_layer_potential_singular_quadratic_midpoints(area, triangle, 4)
single_layer_potential_singular_quadratic_midpoints(area, triangle, 5)
single_layer_potential_singular_quadratic_midpoints(area, triangle, 6)