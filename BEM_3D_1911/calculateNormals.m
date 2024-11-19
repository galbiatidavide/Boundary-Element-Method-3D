function normals = calculateNormals(meshData)
    % Calcola i vettori normali uscenti per ogni triangolo della mesh
    %
    % INPUT:
    % meshData.ConnectivityList - Connettività degli elementi (triangoli)
    % meshData.Points - Coordinate dei punti
    %
    % OUTPUT:
    % normals - Matrice (numTriangles x 3) contenente i vettori normali uscenti
    
    % Numero di triangoli
    numTriangles = size(meshData.ConnectivityList, 1);
    
    % Calcoliamo il centro della mesh (media di tutte le coordinate dei punti)
    meshCenter = mean(meshData.Points, 1);
    
    % Inizializziamo la matrice per i vettori normali (numTriangles x 3)
    normals = zeros(numTriangles, 3);
    
    % Loop su ogni triangolo per calcolare il vettore normale uscente
    for i = 1:numTriangles
        % Otteniamo gli indici dei vertici del triangolo
        vertexIndices = meshData.ConnectivityList(i, :);
        
        % Estraiamo le coordinate dei tre vertici del triangolo
        v1 = meshData.Points(vertexIndices(1), :);
        v2 = meshData.Points(vertexIndices(2), :);
        v3 = meshData.Points(vertexIndices(3), :);
        
        % Calcoliamo i vettori dei lati del triangolo
        vec1 = v2 - v1;
        vec2 = v3 - v1;
        
        % Calcoliamo il prodotto vettoriale (cross product) per ottenere il normale
        normal = cross(vec1, vec2);
        
        % Normalizziamo il vettore normale
        normal = normal / norm(normal);
         
%         % Calcoliamo il centroide del triangolo
%         centroid = (v1 + v2 + v3) / 3;
%         
%         % Creiamo il vettore che va dal centro della mesh al centroide del triangolo
%         vectorToCentroid = centroid - meshCenter;
%         
%         % Verifichiamo se la normale è uscente confrontandola con il vettore centro->centroide
%         % Se la normale punta verso l'interno (verso il centro della mesh), la invertiamo
%         if dot(normal, vectorToCentroid) < 0
%             normal = -normal;
%         end
        
        % Assegniamo il normale corretto alla matrice
        normals(i, :) = normal;
    end
end
