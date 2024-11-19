function [singularity_matrix] = build_singularity_matrix(ConnectivityList, n_points)
    % Initialize an empty adjacency matrix of size n_points x n_points
    singularity_matrix = zeros(n_points, n_points);

    % Iterate over each row (element) in ConnectivityList
    for k = 1:size(ConnectivityList, 1)
        element_points = ConnectivityList(k, :);  % Points in the current element
        
        % For each pair of points in the current element
        for i = 1:length(element_points)
            for j = i+1:length(element_points)
                p1 = element_points(i);
                p2 = element_points(j);
                
                % Set both directions to 1 (p1 -> p2 and p2 -> p1)
                singularity_matrix(p1, p2) = 1;
                singularity_matrix(p2, p1) = 1;
                singularity_matrix(p1, p1) = 1;
                singularity_matrix(p2, p2) = 1;
            end
        end
    end
end
