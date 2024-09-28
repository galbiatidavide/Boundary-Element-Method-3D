function visualizeNormals(meshData)
    % Calcola i vettori normali uscenti per ogni triangolo della mesh
    normals = calculateNormals(meshData);
    
    % Numero di triangoli
    numTriangles = size(meshData.ConnectivityList, 1);
    
    % Inizializziamo la matrice per i centri dei triangoli (numTriangles x 3)
    centers = zeros(numTriangles, 3);
    
    % Calcoliamo i centri di ogni triangolo
    for i = 1:numTriangles
        % Otteniamo gli indici dei vertici del triangolo
        vertexIndices = meshData.ConnectivityList(i, :);
        
        % Estraiamo le coordinate dei tre vertici del triangolo
        v1 = meshData.Points(vertexIndices(1), :);
        v2 = meshData.Points(vertexIndices(2), :);
        v3 = meshData.Points(vertexIndices(3), :);
        
        % Calcoliamo il centro del triangolo come la media dei tre vertici
        centers(i, :) = (v1 + v2 + v3) / 3;
    end
    
    % Visualizzazione della mesh
    trisurf(meshData.ConnectivityList, ...
            meshData.Points(:, 1), meshData.Points(:, 2), meshData.Points(:, 3), ...
            'FaceColor', 'cyan', 'EdgeColor', 'black', 'FaceAlpha', 0.5);
    hold on;
    
    % Visualizzazione dei vettori normali
    quiver3(centers(:, 1), centers(:, 2), centers(:, 3), ...
            normals(:, 1), normals(:, 2), normals(:, 3), ...
            0.5, 'r', 'LineWidth', 1.5); % 0.5 è la scala dei vettori
    
    % Aggiustamenti della visualizzazione
    axis equal;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title('Mesh con Vettori Normali');
    hold off;
end
