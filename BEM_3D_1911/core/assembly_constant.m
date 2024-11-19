function [G, H_hat, H] = assembly_constant(triangles, face_normals, ...
    points, n_elements)

% ASSEMBLY_CONSTANT Assembles the global matrices using constant elements
%
% INPUTS
% triangles: 3D coordinates of the triangles in the triangulation
% face_normals: outward normals of each triangle
% points: the 3D coordinates of the collocation points (usually the
%         baricenters)
% n_elements: total number of elements in the triangulation
%
% OUTPUTS
% G: the matrix which corresponds to the single layer potential
% H_hat: the matrix which corresponds to the double layer potential
% H: H_hat with modified diagonal terms
%    (to obtain H * u_boundary = G * dn_u_boundary)

%% INITIALIZATION OF STRUCTURES
n_points = size(points, 1);
G = zeros(n_points, n_elements);
H_hat = zeros(n_points, n_elements);

%% LOOP OVER THE ELEMENTS
for j = 1:n_elements 
    % selecting current triangle
    triangle = squeeze(triangles(j,:,:));
    % computing integrals of the Green function and its gradient
    [int_green_fun, int_green_grad] = int_green3d_tri(points, ...
        triangle);
    int_green_fun = int_green_fun / (4*pi);
    int_green_grad = int_green_grad / (4*pi);
    % extracting the normal versor to the current triangle
    normal = face_normals(j,:)';
    % extracting the component in the normal direction
    int_green_directional_der = int_green_grad * normal;
    % filling the matrices
    G(:,j) = int_green_fun;
    H_hat(:,j) = int_green_directional_der;
end

%% FINAL COMPUTATIONS
if (nargout == 3)
    H = H_hat - 0.5 * eye(n_elements);
    H = - H;
end
H_hat = - H_hat;

end

