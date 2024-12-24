function [alpha] = solid_angle(points, triangles, face_normals, n_elements)

% SOLID_ANGLE computes the solid angle in a point of the boundary
% 
% INPUTS
% points: the 3D coordinates of the collocation points 
% triangles: 3D coordinates of the triangles in the triangulation
% face_normals: outward normals of each triangle
% n_elements: total number of elements in the triangulation
%
% OUTPUTS
% alpha: the computed solid angle

alpha = 0.0;
for j = 1:n_elements 
    % selecting current triangle
    triangle = squeeze(triangles(j,:,:));
    % computing integral of the Green function's gradient
    [~, int_green_grad] = int_green3d_tri(points, triangle);
    int_green_grad = int_green_grad / (4*pi);
    % extracting the normal versor to the current triangle
    normal = face_normals(j,:)';
    % extracting the component in the normal direction
    int_green_directional_der = int_green_grad * normal;
    alpha = alpha + int_green_directional_der;
end
alpha = alpha * (4*pi);

end

