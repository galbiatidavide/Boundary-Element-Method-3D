function outward_normals = compute_normals(meshData)
    % Funzione per correggere le normali in base alla loro direzione rispetto alla mesh
    %
    % INPUT:
    % meshData.ConnectivityList - Connettività degli elementi (triangoli)
    % meshData.Points - Coordinate dei punti
    %
    % OUTPUT:
    % outward_normals - Normali corrette che puntano verso l'esterno

    % Compute normal directions
    normals = calculate_normal_directions(meshData);

    % Numero di triangoli
    numTriangles = size(meshData.ConnectivityList, 1);

    % Inizializziamo la matrice per le normali corrette
    outward_normals = zeros(numTriangles, 3);

    % Estrazione dei triangoli e delle coordinate dei punti
    tri_connectivity = meshData.ConnectivityList(:, 1:3);  % Consideriamo solo i primi 3 vertici
    points = meshData.Points;  % Coordinate dei punti

    % Loop su ogni triangolo
    for i = 1:numTriangles
        % Otteniamo gli indici dei vertici del triangolo
        vertexIndices = tri_connectivity(i, :);

        % Estraiamo le coordinate dei tre vertici del triangolo
        v1 = points(vertexIndices(1), :);
        v2 = points(vertexIndices(2), :);
        v3 = points(vertexIndices(3), :);

        % Calcoliamo il centroide del triangolo
        element_centroid = (v1 + v2 + v3) / 3;
        
        % Use as epsilon 0.1*the shortest border of the triangle
        epsilon = 0.1*min([norm(v1 - v2), norm(v2 - v3), norm(v3 - v1)]);

        % Generiamo un punto vicino al triangolo lungo la normale
        test_point = element_centroid + normals(i, :) * epsilon;

        % Verifichiamo se il punto si trova dentro o fuori dal dominio
        if point_inside_mesh(test_point, tri_connectivity, points, normals(i, :))
            % Se il punto è dentro il dominio, invertiamo la normale
            outward_normals(i, :) = -normals(i, :);
        else
            % Altrimenti, la normale rimane invariata
            outward_normals(i, :) = normals(i, :);
        end
    end
end