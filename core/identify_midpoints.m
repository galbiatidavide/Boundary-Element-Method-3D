function [is_mid_point] = identify_midpoints(triangles, points)

% Initialize the midpoint indicator vector
num_points = size(points, 1);
is_mid_point = zeros(num_points, 1);

% Loop through each triangle to mark midpoints
num_triangles = size(triangles, 1);
for j = 1:num_triangles
    % Extract the indices of the last 3 points (midpoints) in this triangle
    triangle = squeeze(triangles(j, :, :));
    
    % Check each of the last 3 points in the current triangle
    for k = 4:6
        % Find the row in points corresponding to this midpoint
        point_idx = find(ismember(points, triangle(k, :), 'rows'));
        
        % Mark this point as a midpoint
        if ~isempty(point_idx)
            is_mid_point(point_idx) = 1;
        end
    end
end
end