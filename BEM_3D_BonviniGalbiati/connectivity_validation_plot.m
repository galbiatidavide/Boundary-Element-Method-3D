% Extract the attributes from the TR struct
connectivity_list = TR.ConnectivityList; % The connectivity of the mesh
points = TR.Points;                         % The physical coordinates of the points

% Display sizes
disp('Connectivity List Size:');
disp(size(connectivity_list)); % Should be [24, 6]
disp('Points Size:');
disp(size(points));            % Should be [50, 3]

% Check unique indices in the connectivity list
unique_indices = unique(connectivity_list);
disp('Unique Indices in Connectivity List:');
disp(unique_indices);

% Validate index range
if any(unique_indices < 1) || any(unique_indices > size(points, 1))
    error('Invalid indices detected in the connectivity list.');
else
    % Create the surface plot
    figure;
    trisurf(connectivity_list, points(:,1), points(:,2), points(:,3));

    % Set plot properties
    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel('Z-axis');
    title('3D Mesh Plot');
    axis equal; % Maintain aspect ratio
    view(3);    % Set the view to 3D
    grid on;    % Turn on the grid
end
