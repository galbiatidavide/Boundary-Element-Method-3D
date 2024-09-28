function [u_boundary_nodal, dn_u_boundary_nodal, TR, u_boundary, ...
    dn_u_boundary] = BEM_collocation_quadratic(data)

% BEM_COLLOCATION_LINEAR Main code of the BEM with linear elements and
%                      collocation
%
% INPUTS
% data: the struct which is the output of the function create_data() in 
%       create_data.m
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
close all
clc
addpath('core')
addpath('mesh')
addpath('int_green3d-1.1')
fprintf('\n')
fprintf('-- LAPLACE EQUATION: 3D BEM WITH COLLOCATION AND LINEAR ELEMENTS --\n')
fprintf('\n')

%% PREPROCESSING
% 
% % INITIAL CHECK
% if length(data.filename) ~= 1
%     error('Only one STL file at a time must be provided');
% else
%     data.filename = data.filename{1};
% end

% READING PROVIDED MESH FILE
% fprintf('-- IMPORTING MESH FILE --\n')
% tic
% [TR,fileformat,attributes,solidID] = stlread(data.filename);
% toc
% fprintf(data.filename);
% fprintf(' read!\n');
% if data.enable_postprocessing ~= 0
%     figure
%     trimesh(TR, 'EdgeColor', 'b')
%     axis equal
%     title('Imported mesh')
% end

% This file is used to load the mesh generated with GMSH in a .msh file.
% Since MATLAB does not support quadratic elements this mesh is generated
% in python with GMSH and then the connectivity, points and triangles data
% are loaded separately as .mat files.

% Load the .mat file into MATLAB
TR = load('mesh_data_struct.mat');

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
face_normals = faceNormal(TR);
toc
fprintf('-- END PREPROCESSING --\n\n')

%% CORE

% ASSEMBLY PHASE: COMPUTING AND STORING THE INTEGRALS OF THE 
% GREEN FUNCTION AND ITS DERIVATIVE
fprintf('-- ASSEMBLING GLOBAL MATRICES --\n')
tic
[G, H_hat, H, alpha] = assembly_linear(triangles, face_normals, nodes, ...
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

    % PLOT THE SOLUTION OF THE NORMAL DERIVATIVE OF U AT THE BOUNDARY
    figure
    trisurf(TR, dn_u_boundary);
    colorbar
    axis equal
    title('Normal derivative of the solution at the boundary')
    hold off

    % PLOT THE SOLUTION OF U AT THE BOUNDARY
    figure
    trisurf(TR, u_boundary)
    colorbar
    axis equal
    title('Solution at the boundary')
    hold off

    % COMPUTE AND PLOT THE SOLUTION OF U IN THE INTERNAL DOMAIN ON A X-Y
    % PLANE
    postprocessing_domain(u_boundary_nodal, dn_u_boundary_nodal, ...
        triangles, face_normals, n_elements, TR, ...
        data.enable_postprocessing(2), 'linear');

    toc
    fprintf('-- END POSTPROCESSING --\n')
end

end

