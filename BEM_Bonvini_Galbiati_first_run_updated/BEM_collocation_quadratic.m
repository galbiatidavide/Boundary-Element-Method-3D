function [u_boundary_nodal, dn_u_boundary_nodal, TR, u_boundary, ...
    dn_u_boundary] = BEM_collocation_quadratic(data)

% BEM_COLLOCATION_QUADRATIC Main code of the BEM with quadratic elements and
%                      collocation
%
% INPUTS
% data: the struct which is the output of the function create_data_quadratic() in 
%       create_data_quadratic.m
%
% OUTPUTS
% u_boundary_nodal: the solution at the boundary (nodal values)
% u_n_boundary_nodal: the normal derivative of the solution at the boundary
%                     (nodal values)
% TR: the triangulation
% u_boundary: the solution at the boundary (element values)
% u_n_boundary: the normal derivative of the solution at the boundary
%               (element values)

%%
fprintf('\n')
fprintf('-- LAPLACE EQUATION: 3D BEM WITH COLLOCATION AND QUADRATIC ELEMENTS --\n')
fprintf('\n')

%% PREPROCESSING

% This file is used to load the mesh generated with GMSH in a .msh file.
% Since MATLAB does not support quadratic elements this mesh is generated
% in python with GMSH and then the connectivity, points and triangles data
% are loaded separately as .mat files.

% To do so you need to run the files 'mesh_generator.ipynb' and
% 'mesh_data_struct_generator.ipynb' that are in the folder 'mesh'.

% READING PROVIDED MESH FILE
fprintf('-- IMPORTING MESH FILE --\n')
tic
% Load the .mat file into MATLAB
TR = load('mesh_data_struct.mat');
toc

% COMPUTING BARICENTERS COORDINATES
fprintf('-- COMPUTING COORDINATES --\n')
tic
% for i = 1:n_elements
%     % extracting 3D coordinates of the vertices of each flat triangle
%     triangles(i,:,:) = [TR.Points(TR.ConnectivityList(i,1),:); ...
%                         TR.Points(TR.ConnectivityList(i,2),:); ...
%                         TR.Points(TR.ConnectivityList(i,3),:)];
% end
triangles = TR.triangles;
nodes = TR.Points;
n_elements = size(TR.ConnectivityList,1);
% computing normals to each face of the triangulation
face_normals = calculateNormals(TR);
toc
fprintf('-- END PREPROCESSING --\n\n')

%% CORE

% ASSEMBLY PHASE: COMPUTING AND STORING THE INTEGRALS OF THE 
% GREEN FUNCTION AND ITS DERIVATIVE
fprintf('-- ASSEMBLING GLOBAL MATRICES --\n')
tic
[G, H_hat, H, alpha] = assembly_quadratic(triangles, face_normals, nodes, ...
    n_elements, TR);
toc

% ASSIGNING BOUNDARY CONDITIONS (FULLY DIRICHLET PROBLEM)
fprintf('-- ASSIGNING BOUNDARY CONDITIONS --\n')
dirichlet = data.dirichlet;
neumann = data.neumann;
tic
[u_boundary_nodal, dn_u_boundary_nodal, bc_type] = assign_bcs('mixed', ...
    dirichlet, neumann, nodes);
toc

% GETTING THE SOLUTION OF THE LINEAR SYSTEM
fprintf('-- SOLUTION OF THE LINEAR SYSTEM --\n')
tic
[A, rhs] = recombine_matrices(G, H, u_boundary_nodal, dn_u_boundary_nodal, bc_type);
P = diag(diag(A));
A = P\A;
rhs = P\rhs;
if data.enable_iterative_solver(1) ~= 0
    tol = data.enable_iterative_solver(2);
    max_iterations = data.enable_iterative_solver(3);
    sol = gmres(A,rhs,[],tol,max_iterations);
else
    sol = A\rhs;
end

[u_boundary_nodal, dn_u_boundary_nodal] = get_solution_boundary(u_boundary_nodal, ...
    dn_u_boundary_nodal, sol, bc_type);
u_boundary = compute_elements_values_linear(TR, u_boundary_nodal);
dn_u_boundary = compute_elements_values_linear(TR, dn_u_boundary_nodal);
toc
fprintf('-- END OF CORE COMPUTATIONS --\n')

%% POST PROCESSING

if data.enable_postprocessing(1) ~=0
    fprintf('\n-- POSTPROCESSING --\n')
    tic

    % Linearize the triangulation (discard midpoints) to pass it to trisurf
    quad_elements = TR.ConnectivityList;  % Connectivity list of quadratic elements
    quad_points = TR.Points;  % Node coordinates
    
    % Extract only the first three nodes (1, 2, 3) of each quadratic element to create linear elements
    linear_elements = quad_elements(:, 1:3);  % Keep only the vertex nodes
    
    % Create a new triangulation object with the linearized elements
    TR_linearized = triangulation(linear_elements, quad_points);

    % PLOT THE SOLUTION OF THE NORMAL DERIVATIVE OF U AT THE BOUNDARY
    figure
    trisurf(TR_linearized, dn_u_boundary);
    colorbar
    axis equal
    title('Normal derivative of the solution at the boundary')
    hold off

    % PLOT THE SOLUTION OF U AT THE BOUNDARY
    figure
    trisurf(TR_linearized, u_boundary)
    colorbar
    axis equal
    title('Solution at the boundary')
    hold off

    % COMPUTE AND PLOT THE SOLUTION OF U IN THE INTERNAL DOMAIN ON A X-Y
    % PLANE
    postprocessing_domain(u_boundary_nodal, dn_u_boundary_nodal, ...
        triangles, face_normals, n_elements, TR, ...
        data.enable_postprocessing(2), 'quadratic');

    toc
    fprintf('-- END POSTPROCESSING --\n')
end

end

