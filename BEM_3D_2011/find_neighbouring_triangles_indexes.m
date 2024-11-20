function [triangle_indices] = find_neighbouring_triangles_indexes(triangles, midpoint)
    % Initialize an empty vector to store the indices of triangles
    triangle_indices = [];
    
    % Number of triangles in the mesh
    num_triangles = size(triangles, 1);
    
    % Loop over each triangle to check for the midpoint
    for j = 1:num_triangles
        % Extract the 6x3 matrix for the current triangle
        triangle = squeeze(triangles(j, :, :));
        
        % Check the last three points in the triangle (midpoints)
        for k = 4:6
            if isequal(triangle(k, :), midpoint)
                % If the midpoint is found, add the index of this triangle
                triangle_indices = [triangle_indices; j]; %#ok<AGROW> 
                break;  % Stop checking further points in this triangle
            end
        end
        
        % Stop if we have already found two triangles
        if length(triangle_indices) == 2
            break;
        end
    end
end
