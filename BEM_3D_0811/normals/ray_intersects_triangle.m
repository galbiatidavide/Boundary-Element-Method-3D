function intersects = ray_intersects_triangle(origin, direction, v1, v2, v3)
    epsilon = 1e-8;
    
    % Compute the edges of the triangle
    edge1 = v2 - v1;
    edge2 = v3 - v1;

    % Compute the determinant
    h = cross(direction, edge2);
    a = dot(edge1, h);

    if abs(a) < epsilon
        intersects = false;  % The ray is parallel to the triangle
        return;
    end

    f = 1.0 / a;
    s = origin - v1;
    u = f * dot(s, h);

    if u < 0.0 || u > 1.0
        intersects = false;  % Intersection point is outside the triangle
        return;
    end

    q = cross(s, edge1);
    v = f * dot(direction, q);

    if v < 0.0 || u + v > 1.0
        intersects = false;  % Intersection point is outside the triangle
        return;
    end

    % Compute the distance along the ray to the intersection point
    t = f * dot(edge2, q);

    if t > epsilon
        intersects = true;
    else
        intersects = false;
    end
end
