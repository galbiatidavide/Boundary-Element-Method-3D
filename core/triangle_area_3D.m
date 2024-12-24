function [area] = triangle_area_3D(A, B, C)
    % Calculate vectors AB and AC
    AB = B - A;
    AC = C - A;
    
    % Compute the cross product of AB and AC
    cross_product = cross(AB, AC);
    
    % Calculate the area of the triangle
    area = 0.5 * norm(cross_product);
end