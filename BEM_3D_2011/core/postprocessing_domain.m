function postprocessing_domain(u_boundary_nodal, dn_u_boundary_nodal, ...
    triangles, face_normals, n_elements, TR, z_coord, which_bem)

% POSTPROCESSING_DOMAIN Post-processing phase (domain on a x-y plane given 
%                       by the user)
% INPUTS
% u_boundary_nodal: the solution at the boundary at the nodal values
% dn_u_boundary_nodal: the normal derivative of solution at the boundary
%                      at the nodal values
% triangles: the triangles in the triangulation
% face_normals: the triangulation normals
% n_elements: the total number of elements
% TR: the triangulation
% z_coord: the z coordinate which defines the plane
% which_bem: 'constant', 'linear' or 'quadratic'

% generating points internal to the domain with z fixed
points_per_dimension = 20;
points = TR.Points;
epsilon = 1e-1*max(abs(points)); % to guarantee that you don't keep points 
                                 % on the boundary
coordinate_extremes = [min(points,[],1) + epsilon; max(points,[],1) - epsilon];
x_coord = linspace(coordinate_extremes(1,1), coordinate_extremes(2,1), ...
    points_per_dimension);
y_coord = linspace(coordinate_extremes(1,2), coordinate_extremes(2,2), ...
    points_per_dimension);
[XX, YY] = meshgrid(x_coord, y_coord);
XX = reshape(XX, points_per_dimension * points_per_dimension, 1, 1);
YY = reshape(YY, points_per_dimension * points_per_dimension, 1, 1);
ZZ = z_coord * ones(length(XX), 1);
internal_points = [XX YY ZZ];

% evaluating the solution at the internal points: choice between constant 
% or linear
if (strcmp(which_bem, 'constant'))
    solution = get_solution_domain_constant(u_boundary_nodal, ...
        dn_u_boundary_nodal, triangles, face_normals, internal_points, ...
        n_elements);
elseif (strcmp(which_bem, 'linear'))
    solution = get_solution_domain_linear(u_boundary_nodal, ...
        dn_u_boundary_nodal, triangles, face_normals, internal_points, ...
        n_elements, TR);
elseif (strcmp(which_bem, 'quadratic'))
    solution = get_solution_domain_quadratic(u_boundary_nodal, ...
        dn_u_boundary_nodal, triangles, face_normals, internal_points, ...
        n_elements, TR);
else
    error('Inexistent kind of BEM chosen');
end

% plotting solution on a x-y plane
figure
XX = reshape(XX, points_per_dimension, points_per_dimension);
YY = reshape(YY, points_per_dimension, points_per_dimension);
solution = reshape(solution, points_per_dimension, points_per_dimension);
surf(XX, YY, solution)
colorbar
xlabel('x')
ylabel('y')
zlabel(['u(x,y,', num2str(z_coord),')'])
title(['Solution on the plane: z = ', num2str(z_coord)])

end

