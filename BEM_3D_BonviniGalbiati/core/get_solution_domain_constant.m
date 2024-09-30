function [solution] = get_solution_domain_constant(u_boundary, dn_u_boundary, ...
    triangles, face_normals, points, n_elements)

% GET_SOLUTION_DOMAIN computes the solution inside the domain, given
%                     the solution at the boundary
% INPUTS
% u_boundary: solution at the boundary
% dn_u_boundary: normal derivative of the solution at the boundary
%                (convention: outward normal)
% triangles: 3D coordinates of the triangles in the triangulation
% face_normals: outward normals of each triangle
% points: the 3D coordinates of the collocation points (usually the
%         baricenters)
% n_elements: total number of elements in the triangulation
%
% OUTPUTS
% solution: the value of u evaluated at the points provided as argument

[G, H_hat] = assembly_constant(triangles, face_normals, points, n_elements);
single_layer_potential = G * dn_u_boundary;
double_layer_potential = H_hat * u_boundary;
solution = - double_layer_potential + single_layer_potential;

end

