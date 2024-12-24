function [G, H_hat, H, alpha] = assembly_quadratic_high_efficiency(triangles, face_normals, ...
    points, n_elements, TR)
%% INITIALIZATION OF STRUCTURES
n_points = size(points,1);
G = zeros(n_points,n_points);
H_hat = zeros(n_points,n_points);

%% define the quadrature
quadrature_nodes = [0 0; 1 0; 0 1; 0.5 0; 0.5 0.5; 0 0.5];
weights = [0, 0, 0, 1/3, 1/3, 1/3];

% quadrature_nodes = [0.5 0; 0.5 0.5; 0 0.5];
% weights = [1/3, 1/3, 1/3];

%% define the basis functions and evaluate them on the quadrature nodes
basis_functions = @(csi) ...
    [(1-csi(1)-csi(2)).*(1-2*csi(1)-2*csi(2)); ...
    csi(1).*(2*csi(1)-1); ...
    csi(2).*(2*csi(2)-1); ...
    4*csi(1).*(1-csi(1)-csi(2)); ...
    4*csi(1).*csi(2); ...
    4*csi(2).*(1-csi(1)-csi(2))];

% basis_functions_eval(i,j) has the evaluation of basis_function(j) in
% quadrature_node(i)
basis_functions_eval = zeros(size(quadrature_nodes,1), 6);
for q = 1:size(quadrature_nodes,1)
   basis_functions_eval(q,:) = basis_functions(quadrature_nodes(q,:));
end
%% build the singularity matrix 
% (a boolean matrix that gives me 1 where I will have singular integrals and 0 otherwise)
is_singular = build_singularity_matrix(TR.ConnectivityList, n_points);

%% assemble the matrices
for i = 1:n_points

    % get the coordinates of point i
    x_point = points(i,1);
    y_point = points(i,2);
    z_point = points(i,3);

    for j = 1:n_elements
        
        if point_belongs_to_elem(i,j,points,triangles) == false % skip sigular integral calculations

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
            % midpoints
            x4 = triangle(4,1);
            x5 = triangle(5,1);
            x6 = triangle(6,1);
            y4 = triangle(4,2);
            y5 = triangle(5,2);
            y6 = triangle(6,2);
            z4 = triangle(4,3);
            z5 = triangle(5,3);
            z6 = triangle(6,3);
        
            % compute the area of the triangle (det(Jacobian))
            area = triangle_area_3D([x1 y1 z1], [x2 y2 z2], [x3 y3 z3]);
        
            % get the iglo
            iglo = TR.ConnectivityList(j,:);
        
            % get the normal direction
            normal = face_normals(j,:)';
    
            for l = 1:6
                temp_H_hat = 0;
                temp_G = 0;
    
                for h = 1:size(quadrature_nodes,1)
    
                    % parametrize the distances
                    x_param_eval = ([x1 x2 x3 x4 x5 x6]-x_point) * basis_functions_eval(h,:)';
                    y_param_eval = ([y1 y2 y3 y4 y5 y6]-y_point) * basis_functions_eval(h,:)';
                    z_param_eval = ([z1 z2 z3 z4 z5 z6]-z_point) * basis_functions_eval(h,:)';
                    r = sqrt(x_param_eval^2 + y_param_eval^2 + z_param_eval^2);
    
                    % evaluate the fundamental solution v and its normal
                    % derivative v_n
                    v = 1/r; % green fun
                    v_n = -([x_param_eval y_param_eval z_param_eval]*normal)*(1/(r^3)); % green der
                    
                    temp_H_hat = temp_H_hat + basis_functions_eval(h,l)*v_n*weights(h)*area;
                    temp_G = temp_G + basis_functions_eval(h,l)*v*weights(h)*area;
    
                end
                
                H_hat(i,iglo(l)) = H_hat(i,iglo(l)) + (1/(8*pi))*temp_H_hat*2;
                G(i,iglo(l)) = G(i,iglo(l)) + (1/(8*pi))*temp_G*2;
    
            end
        end
    end
end

%% singular integrals
for j=1:n_elements

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
    iglo = TR.ConnectivityList(j,:);

    for local_index = 1:6
        if local_index <= 3
         G(iglo(local_index),iglo) = G(iglo(local_index),iglo) + ...
             single_layer_potential_singular_quadratic_vertices(area, triangle, local_index);
        else
         G(iglo(local_index),iglo) = G(iglo(local_index),iglo) + ...
            single_layer_potential_singular_quadratic_midpoints(area, triangle, local_index);
        end
    end
end

%% Computation of H with solid angle
if nargout >= 3
    alpha = solid_angle_quadratic(points, triangles, face_normals, n_elements);
    H = H_hat + diag(alpha/(4*pi));
    assignin("base", "H_hat",H_hat);
    assignin("base","H",H);
    assignin("base","G",G);
    assignin("base", "is_singular", is_singular);
end

end