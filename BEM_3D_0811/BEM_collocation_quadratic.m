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
TR = data.TR;
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
face_normals = compute_normals(TR);
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
[u_boundary_nodal, dn_u_boundary_nodal, bc_type] = assign_bcs(data.bc_type, ...
    dirichlet, neumann, nodes);
toc

% GETTING THE SOLUTION OF THE LINEAR SYSTEM
fprintf('-- SOLUTION OF THE LINEAR SYSTEM --\n')
tic
[A, rhs] = recombine_matrices(G, H, u_boundary_nodal, dn_u_boundary_nodal, bc_type);
% P = diag(diag(A));
% A = P\A;
% rhs = P\rhs;
% if data.enable_iterative_solver(1) ~= 0
%     tol = data.enable_iterative_solver(2);
%     max_iterations = data.enable_iterative_solver(3);
%     sol = gmres(A,rhs,[],tol,max_iterations);
% else
%     sol = A\rhs;
% end
% Block Jacobi Preconditioner
% Assuming A is block-structured, divide A into blocks and invert each block

% Assume you have a vector of indices for Neumann boundary conditions
% Recombine matrices (unchanged)
% block_size = size(G, 1);  % Assume block size is provided in data structure
% num_blocks = size(A, 1) / block_size;  % Number of blocks
% 
% % Create a preconditioner function
% function x = block_gauss_seidel(A, b, block_size, max_iterations, neumann_indices)
%     x = zeros(size(b));  % Initial guess
%     for iter = 1:max_iterations
%         for i = 1:num_blocks
% 
%             idx = (i-1)*block_size + 1:i*block_size;  % Indices for the current block            
%             % Skip Neumann indices
%             if (neumann_indices(idx)=='N')
%                 % Calculate the residual
%                 r = b(idx) - A(idx, :) * x;  % Compute residual using current x
%                 % Update the current block using the lower triangular part of A
%                 x(idx) = x(idx) + A(idx, idx) \ r;  % Solve for x_i
%             end
%         end
%     end
% end
% 
% % Preconditioning with Block Gauss-Seidel
% max_iterations = 200;  % Number of iterations for the preconditioner
% x0 = zeros(size(rhs));  % Initial guess for preconditioner
% sol_preconditioned = block_gauss_seidel(A, rhs, block_size, max_iterations, bc_type);
% 
% % Apply the preconditioner to the rhs
% % Skip Neumann indices in the rhs preconditioning
% rhs_preconditioned = rhs;
% for i = 1:length(rhs)
%     if (bc_type(i) == 'D')
%         rhs_preconditioned(i) = rhs(i) - A(i, :) * sol_preconditioned;  % Adjust rhs only for non-Neumann conditions
%     end
% end

% Solver configuration
if data.enable_iterative_solver(1) ~= 0
    tol = data.enable_iterative_solver(2);
    max_iterations = data.enable_iterative_solver(3);

    % Solve the preconditioned system using BiCG
    sol = gmres(A, rhs,[],tol, max_iterations);
    % x0 = 1e6*ones(size(A,1),1);
    % sol = bicg(A,rhs,[],max_iterations,[],[],x0);



else
    % Direct solver if iterative solver is not enabled
%     P = diag(diag(A));
%     A = P\A;
%     rhs = P\rhs;
%     sol = A\rhs;
    % sol = lsqr(A,rhs,1e-6,100000);
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
    
    % Find the unique vertex nodes (1, 2, 3) and their coordinates
    unique_vertex_indices = unique(linear_elements(:));  % Unique vertex indices from linearized elements
    linear_points = quad_points(unique_vertex_indices, :);  % Coordinates of vertex nodes
    
    % Adjust the element connectivity to reference the reduced set of points
    [~, new_indices] = ismember(linear_elements, unique_vertex_indices);  % Map old indices to new ones
    TR_linearized = triangulation(new_indices, linear_points);  % Create linearized triangulation
    
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
%     postprocessing_domain(u_boundary_nodal, dn_u_boundary_nodal, ...
%         triangles, face_normals, n_elements, TR, ...
%         data.enable_postprocessing(2), 'quadratic');
    toc
    fprintf('-- END POSTPROCESSING --\n')
end

end

