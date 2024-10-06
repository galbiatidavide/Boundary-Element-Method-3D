function is_inside = point_inside_mesh(p_test, tri_connectivity, points, ray_direction)
    
    num_intersections = 0;  % Count the number of intersections

    for i = 1:size(tri_connectivity, 1)
        % Get the vertices of the current triangle (discarding midpoints)
        tri_nodes = tri_connectivity(i, :);
        v1 = points(tri_nodes(1), :);
        v2 = points(tri_nodes(2), :);
        v3 = points(tri_nodes(3), :);

        % Perform ray-triangle intersection check
        if ray_intersects_triangle(p_test, ray_direction, v1, v2, v3)
            num_intersections = num_intersections + 1;
        end
    end

    % If the number of intersections is odd, the point is inside; otherwise, outside
    is_inside = mod(num_intersections, 2) == 1;
end
