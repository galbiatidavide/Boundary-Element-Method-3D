function [G, H_hat, H, alpha] = assembly_quadratic(triangles, face_normals, ...
    points, n_elements, TR)

% ASSEMBLY_LINEAR Assembles the global matrices using linear elements
%                 (vectorization enabled, ONLY TRIANGULAR MESHES WITH FLAT 
%                  TRIANGLES!)
%
% INPUTS
% triangles: 3D coordinates of the triangles in the triangulation
% face_normals: outward normals of each triangle
% points: 3D coordinates of the collocation points
% n_elements: total number of elements in the triangulation
% TR: the triangulation
%
% OUTPUTS
% G: the matrix which corresponds to the single layer potential
% H_hat: the matrix which corresponds to the double layer potential
% H: H_hat with modified diagonal terms 
%    (to obtain H * u_boundary = G * dn_u_boundary)
% alpha: the interior solid angle at each node

%% INITIALIZATION OF STRUCTURES
n_points = size(points,1);
nodes = TR.Points;

n_nodes = size(nodes,1);
is_boundary_system = (n_nodes == n_points); % checks if the nodes are
                                            % the same as the points
if is_boundary_system
    is_boundary_system = (nodes == points);
end
H_hat = zeros(n_points, n_nodes);
G = zeros(n_points, n_nodes);
%basis_functions = @(csi) [1 - csi(1) - csi(2); csi(1); csi(2)];
basis_functions = @(csi) ...
    [(1-csi(1)-csi(2)).*(1-2*csi(1)-2*csi(2)); ...
    csi(1).*(2*csi(1)-1); ...
    csi(2).*(2*csi(2)-1); ...
    4*csi(1).*(1-csi(1)-csi(2)); ...
    4*csi(1).*csi(2); ...
    4*csi(2).*(1-csi(1)-csi(2))];


% quadrature_nodes = [0.5 0;
%                      0.5 0.5;
%                      0 0.5];

%% IO LAVOREREI SUL TRIANGOLO DI RIFERIMENTO DI AREA 0.5 E NON SU QUELLO DI QUADRATURA
%% perche' la nostra TR.ConnectivityList ha 6 colonne ed e' tutto piu' coerente.
quadrature_nodes = [0 0; 1 0; 0 1; 0.5 0; 0.5 0.5; 0 0.5];
basis_functions_eval = zeros(size(quadrature_nodes,1), 6);
%% MATRICE DIAGONALE: ok
for q = 1:size(quadrature_nodes,1)
   basis_functions_eval(q,:) = basis_functions(quadrature_nodes(q,:));
end



%% CONSTRUCT DATA STRUCTURES TO ENABLE VECTORIZATION

% creating a block diagonal matrix (each block is equal to
% basis_functions_eval)
basis_functions_eval_helper = repmat({basis_functions_eval}, 1, n_points);
basis_functions_eval_rep = blkdiag(basis_functions_eval_helper{:});
% each *_points_rep (where * can be x, y or z) is a column vector containing
% {[*_i; *_i; *_i]} for i=1:n_points. For instance:
% x_nodes_rep = [x_1; x_1; x_1; x_2; x_2; x_2; ...; x_{n_points};
%                x_{n_points}; x_{n_points}]
x_points_rep = repmat(points(:,1), 1, 6);
y_points_rep = repmat(points(:,2), 1, 6);
z_points_rep = repmat(points(:,3), 1, 6);
x_points_rep = x_points_rep';
y_points_rep = y_points_rep';
z_points_rep = z_points_rep';
x_points_rep = reshape(x_points_rep, 6*n_points, 1);
y_points_rep = reshape(y_points_rep, 6*n_points, 1);
z_points_rep = reshape(z_points_rep, 6*n_points, 1);

%% LOOP OVER THE ELEMENTS
for j = 1:n_elements
    % extracting information about the considered element
    iglo = TR.ConnectivityList(j,:);
    triangle = squeeze(triangles(j,:,:));
    normal = face_normals(j,:)';
    x = triangle(:,1);
    y = triangle(:,2);
    z = triangle(:,3);

    x_vertex = x(1:3);
    y_vertex = y(1:3);
    z_vertex = z(1:3);
    
    % constructing a data structure to enable vectorization: *_tria_rep
    % is a column vector which contains [*_{1}, *_{2}, *_{3}] 
    % concatenated n_points times; the index (which goes from 1 to 3) 
    % represent the local index of the vertexes of the considered element
    % For instance:
    % x_tria_rep = [x_1; x_2; x_3; x_1; x_2; x_3; ...; x_1; x_2; x_3];
    x_tria_rep = repmat(x, n_points, 1);
    y_tria_rep = repmat(y, n_points, 1);
    z_tria_rep = repmat(z, n_points, 1);
    
    % values for mid-point quadrature formula in a triangle
    %% la sua formula e' giusta: area = area del triangolo fisico (Jacobiano).
    ons = [1 1 1];

    area = 0.5 * sqrt(det([x_vertex';y_vertex';ons])^2 + ...
        det([y_vertex';z_vertex';ons])^2 + ...
        det([z_vertex';x_vertex';ons])^2);

 
    % parametrized coordinates for midpoint quadrature formula in a
    % triangle
    %% siccome basis_function_eval e' diagonale la trasformazione e' come moltiplicare per I.
    %% applico solo uno shift spostandomi sulle coordinate di quel triangolo rispetto a tutti gli altri.
    x_param = basis_functions_eval_rep * (x_tria_rep - x_points_rep);
    y_param = basis_functions_eval_rep * (y_tria_rep - y_points_rep); 
    z_param = basis_functions_eval_rep * (z_tria_rep - z_points_rep);
    r_param = sqrt(x_param.^2 + y_param.^2 + z_param.^2); 
    
    % computing the green function and its normal derivative 
   
    green_fun = 1/(4*pi) ./ r_param;
    green_der = - 1/(4*pi) * ([x_param y_param z_param] * normal) ...
                ./(r_param.^3);
  
    %% stampa queste due variabili: green_der contiene solo dei NaN che levi mandanoli a 0 (week singularity)
    %% green_fun invece contiene 6 valori inf vale a dire 6 strong singularitues. Sono le integrazioni quando valuti
    %% nodo p_i sullo stesso nodo p_i.
    % reshaping to obtain a matrix of shape (n_points x 3)
    green_fun = reshape(green_fun, 6, n_points)';
    green_der = reshape(green_der, 6, n_points)';
    
    % setting to 0 the values at the singular nodes, if solving the boundary
    % system
    %metti a 0 i Nan e metti a 0 gli inf che dopo vanno modificati (perdi
    %anche tutti gli altri valori sulla riga con inf ma tanto li
    %ricalcoli con single layer potential)
    
    if is_boundary_system
        green_fun(iglo,:) = 0.0;
        green_der(iglo,:) = 0.0;
    end



    %%Questo potrebbe essere diag(0..., 1/3, 1/3, 1/3) a quanto ho capito se
    %%i pesi sono equivalenti la loro somma deve fare area del triangolo 
    %% => A = 1/3A + 1/3A + 1/3A. Sul triangolo di riferimento avresti A = 0.5
    %% i pesi sono 1/6 e sono tre. Quindi 0.5 = 0.5/? + 0.5/? + 0.5/?. Nel nostro caso A
    %% e' fisica e non vale 0.5. Tuttavia i pesi devono comunque sommarsi e fare A_triangolo fisico.
    %% Io ho capito cosi.

    
    %%I pesi dei vertici sono nulli per avere precisione DI ORDINE TRE: via
    %%le prime tre colonne di green fun. Se interpreto bene green_fun
    %%(50x6) in ogni riga si vede come ogni nodo "influenza tutti i nodi
    %%del j-esimo triangolo). E' come se ogni colonna rappresentasse un
    %%punto di quel triangolo. Infatti quando il ciclo esterno scorre all'elemento
    %%successivo avrai come ogni nodo influenza i nodi del j+1-esimo
    %%triangolo.
    W = 2*diag([0, 0, 0, 1/6, 1/6, 1/6]);
    %W = 1/6*diag([1,1,1,1,1,1]);
    %W = diag([1/30, 1/30, 1/30, 9/30, 9/30, 9/30]); %%DO NOT USE


    %%Se stampi questa variabile vedrai che le prime tre colonne sono
    %%nulle, perche' il loro peso e' nullo.    
    %(green_fun * basis_functions_eval * W);

    
    % updating the global matrixes
    G(:,iglo) = G(:,iglo) + green_fun * basis_functions_eval * W * area; 
    H_hat(:,iglo) = H_hat(:,iglo) + green_der * basis_functions_eval * W * area;

    % correcting G with the singular values, if solving the boundary
    % system (singular values of H are equal to 0).

    %% Ogni nodo del triangolo viene valutato ora:
    if is_boundary_system
        for local_index=1:6
            if local_index <=3
            G(iglo(local_index),iglo) = G(iglo(local_index),iglo) + ...
                single_layer_potential_singular_quadratic_vertices(area, triangle, local_index);
            else
             G(iglo(local_index),iglo) = G(iglo(local_index),iglo) + ...
                single_layer_potential_singular_quadratic_midpoints(area, triangle, local_index);
            end
        end
    end
end



%% COMPUTATION OF H
% computing the solid angle and H matrix if necessary
if nargout >= 3
    alpha = solid_angle_quadratic(nodes, triangles, face_normals, n_elements);
    H = H_hat + diag(alpha/(4*pi));
end

end