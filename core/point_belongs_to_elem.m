function [bool] = point_belongs_to_elem(point_index, element_index, points, triangles)
    % point_belongs_to_elem checks if a given point belongs to a specified element.
    %
    % Inputs:
    % - point_index: index of the point to check (scalar)
    % - element_index: index of the triangle element to check (scalar)
    % - points: matrix of size (n_points x 3) containing (x, y, z) coordinates of each point
    % - triangles: 3D array of size (n_elem x 6 x 3), where each triangle element is defined
    %   by the coordinates of its 6 points in the 3D space
    %
    % Output:
    % - bool: logical value (true if point belongs to the element, false otherwise)
    
    % Get the coordinates of the specified point
    point_coords = points(point_index, :);
    
    % Get the coordinates of all 6 points of the specified element (triangle)
    elem_points = squeeze(triangles(element_index, :, :));
    
    % Initialize boolean as false
    bool = false;
    
    % Check if the point's coordinates match any of the element's points
    for i = 1:6
        if all(point_coords == elem_points(i, :))
            bool = true;
            return; % Exit function as soon as a match is found
        end
    end
end