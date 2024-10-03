function [G, H_hat, H, alpha] = assembly_linear(triangles, face_normals, ...
    points, n_elements, TR)

% ASSEMBLY_LINEAR Assembles the global matrices using linear elements
%                 (vectorization enabled, ONLY TRIANGULAR MESHES WITH FLAT 
%                  TRIANGLES!)
%
% INPUTS
% triangles: 3D coordinates of the triangles in the triangulation
% face_normals: outward normals of each triangle
% points: 3D coordinates of the collocation points
% n_elements: total number of elements in the triangulation
% TR: the triangulation
%
% OUTPUTS
% G: the matrix which corresponds to the single layer potential
% H_hat: the matrix which corresponds to the double layer potential
% H: H_hat with modified diagonal terms 
%    (to obtain H * u_boundary = G * dn_u_boundary)
% alpha: the interior solid angle at each node

%% INITIALIZATION OF STRUCTURES
n_points = size(points,1);
nodes = TR.Points;

n_nodes = size(nodes,1);
is_boundary_system = (n_nodes == n_points); % checks if the nodes are
                                            % the same as the points
if is_boundary_system
    is_boundary_system = (nodes == points);
end
H_hat = zeros(n_points, n_nodes);
G = zeros(n_points, n_nodes);
basis_functions = @(csi) [1 - csi(1) - csi(2); csi(1); csi(2)];
quadrature_nodes = [0.5 0;
                     0.5 0.5;
                     0 0.5];
basis_functions_eval = zeros(size(quadrature_nodes,1), 3);
for q = 1:size(quadrature_nodes,1)
   basis_functions_eval(q,:) = basis_functions(quadrature_nodes(q,:));
end

%% CONSTRUCT DATA STRUCTURES TO ENABLE VECTORIZATION

% creating a block diagonal matrix (each block is equal to
% basis_functions_eval)
basis_functions_eval_helper = repmat({basis_functions_eval}, 1, n_points);
basis_functions_eval_rep = blkdiag(basis_functions_eval_helper{:});
% each *_points_rep (where * can be x, y or z) is a column vector containing
% {[*_i; *_i; *_i]} for i=1:n_points. For instance:
% x_nodes_rep = [x_1; x_1; x_1; x_2; x_2; x_2; ...; x_{n_points};
%                x_{n_points}; x_{n_points}]
x_points_rep = repmat(points(:,1), 1, 3);
y_points_rep = repmat(points(:,2), 1, 3);
z_points_rep = repmat(points(:,3), 1, 3);
x_points_rep = x_points_rep';
y_points_rep = y_points_rep';
z_points_rep = z_points_rep';
x_points_rep = reshape(x_points_rep, 3*n_points, 1);
y_points_rep = reshape(y_points_rep, 3*n_points, 1);
z_points_rep = reshape(z_points_rep, 3*n_points, 1);

%% LOOP OVER THE ELEMENTS
for j = 1:n_elements
    % extracting information about the considered element
    iglo = TR.ConnectivityList(j,:);
    triangle = squeeze(triangles(j,:,:));
    normal = face_normals(j,:)';
    x = triangle(:,1);
    y = triangle(:,2);
    z = triangle(:,3);
    
    % constructing a data structure to enable vectorization: *_tria_rep
    % is a column vector which contains [*_{1}, *_{2}, *_{3}] 
    % concatenated n_points times; the index (which goes from 1 to 3) 
    % represent the local index of the vertexes of the considered element
    % For instance:
    % x_tria_rep = [x_1; x_2; x_3; x_1; x_2; x_3; ...; x_1; x_2; x_3];
    x_tria_rep = repmat(x, n_points, 1);
    y_tria_rep = repmat(y, n_points, 1);
    z_tria_rep = repmat(z, n_points, 1);
    
    % values for mid-point quadrature formula in a triangle
    ons = [1 1 1];
    weight = 1/3;
    area = 0.5 * sqrt(det([x';y';ons])^2 + ...
        det([y';z';ons])^2 + ...
        det([z';x';ons])^2);

 
    % parametrized coordinates for midpoint quadrature formula in a
    % triangle
    x_param = basis_functions_eval_rep * (x_tria_rep - x_points_rep);
    y_param = basis_functions_eval_rep * (y_tria_rep - y_points_rep); 
    z_param = basis_functions_eval_rep * (z_tria_rep - z_points_rep);
    r_param = sqrt(x_param.^2 + y_param.^2 + z_param.^2); 
    
    % computing the green function and its normal derivative 
    green_fun = 1/(4*pi) ./ r_param;
    green_der = - 1/(4*pi) * ([x_param y_param z_param] * normal) ...
                ./(r_param.^3);
            
    % reshaping to obtain a matrix of shape (n_points x 3)
    green_fun = reshape(green_fun, 3, n_points)';
    green_der = reshape(green_der, 3, n_points)';
    
    % setting to 0 the values at the singular nodes, if solving the boundary
    % system
    if is_boundary_system
        green_fun(iglo,:) = 0.0;
        green_der(iglo,:) = 0.0;
    end
        
    % updating the global matrixes
    G(:,iglo) = G(:,iglo) + weight * green_fun * basis_functions_eval * area;
    H_hat(:,iglo) = H_hat(:,iglo) + weight * green_der * basis_functions_eval * area;
    
    % correcting G with the singular values, if solving the boundary
    % system (singular values of H are equal to 0)
    if is_boundary_system
        for local_index=1:3
            G(iglo(local_index),iglo) = G(iglo(local_index),iglo) + ...
                single_layer_potential_singular(area, triangle, local_index);
        end
    end
end

%% COMPUTATION OF H
% computing the solid angle and H matrix if necessary
if nargout >= 3
    alpha = solid_angle(nodes, triangles, face_normals, n_elements);
    H = H_hat + diag(alpha/(4*pi));
end

end