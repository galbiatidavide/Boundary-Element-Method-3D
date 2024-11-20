function [u_boundary, dn_u_boundary, TR] = BEM_collocation_constant(data)

% BEM_COLLOCATION_CONSTANT Main code of the BEM with constant elements and
%                      collocation
%
% INPUTS
% data: the struct which is the output of the function create_data() in 
%       create_data.m
%
% OUTPUTS
% u_boundary: the solution at the boundary
% u_n_boundary: the normal derivative of the solution at the boundary
% TR: the triangulation

%%
close all
clc
addpath('core')
addpath('mesh')
addpath('int_green3d-1.1')
fprintf('\n')
fprintf('-- LAPLACE EQUATION: 3D BEM WITH COLLOCATION AND CONSTANT ELEMENTS --\n')
fprintf('\n')

%% PREPROCESSING

% INITIAL CHECK
if length(data.filename) ~= 1
    error('Only one STL file at a time must be provided');
else
    data.filename = data.filename{1};
end

% READING PROVIDED MESH FILE
fprintf('-- IMPORTING MESH FILE --\n')
tic
[TR,fileformat,attributes,solidID] = stlread(data.filename);
toc
fprintf(data.filename);
fprintf(' read!\n');
if data.enable_postprocessing ~= 0
    figure
    trimesh(TR, 'EdgeColor', 'b')
    axis equal
    title('Imported mesh')
end

% COMPUTING BARICENTERS COORDINATES
fprintf('-- COMPUTING BARICENTERS COORDINATES --\n')
n_elements = length(TR.ConnectivityList);
triangles = zeros(n_elements, 3, 3);
tic
for i = 1:n_elements
    % extracting 3D coordinates of the vertices of each flat triangle
    triangles(i,:,:) = [TR.Points(TR.ConnectivityList(i,1),:); ...
                        TR.Points(TR.ConnectivityList(i,2),:); ...
                        TR.Points(TR.ConnectivityList(i,3),:)];
end
% extract baricenter of each element
baricenters = mean(triangles, 2);
baricenters = squeeze(baricenters);
% computing normals to each face of the triangulation
face_normals = faceNormal(TR);
toc
fprintf('-- END PREPROCESSING --\n\n')


%% CORE

% ASSEMBLY PHASE: COMPUTING AND STORING THE INTEGRALS OF THE 
% GREEN FUNCTION AND ITS DERIVATIVE
fprintf('-- ASSEMBLING GLOBAL MATRICES --\n')
tic
[G, H_hat, H] = assembly_constant(triangles, face_normals, ...
    baricenters, n_elements);
toc

% ASSIGNING BOUNDARY CONDITIONS (FULLY DIRICHLET PROBLEM)
fprintf('-- ASSIGNING BOUNDARY CONDITIONS --\n')
dirichlet = data.dirichlet;
neumann = data.neumann;
tic
[u_boundary, dn_u_boundary, bc_type] = assign_bcs('mixed', ...
    dirichlet, neumann, baricenters);
toc

% GETTING THE SOLUTION OF THE LINEAR SYSTEM
fprintf('-- SOLUTION OF THE LINEAR SYSTEM --\n')
tic
[A, rhs] = recombine_matrices(G, H, u_boundary, dn_u_boundary, bc_type);
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
[u_boundary, dn_u_boundary] = get_solution_boundary(u_boundary, ...
    dn_u_boundary, sol, bc_type);
toc
fprintf('-- END OF CORE COMPUTATIONS --\n')

%% POST PROCESSING

if data.enable_postprocessing(1) ~=0
    fprintf('\n-- POSTPROCESSING --\n')
    tic
    if (size(dn_u_boundary,1) < size(dn_u_boundary,2))
        dn_u_boundary = dn_u_boundary';
    end

    %PLOT THE SOLUTION OF THE NORMAL DERIVATIVE OF U AT THE BOUNDARY
    incenters = incenter(TR);
    figure
    subplot(1,2,1)
    sgtitle('Normal derivatives of the solution at the boundary')
    trisurf(TR, 'FaceColor','cyan','FaceAlpha',0.3);
    dn_u_vector = dn_u_boundary .* abs(face_normals);
    axis equal
    hold on  
    quiver3(incenters(:,1), incenters(:,2), incenters(:,3), ...
         dn_u_vector(:,1), dn_u_vector(:,2), dn_u_vector(:,3), ...
         0.5,'color','r');
    legend('Surface','d_{n}u')
    subplot(1,2,2)
    trisurf(TR, dn_u_boundary);
    colorbar
    axis equal
    hold off

    % PLOT THE SOLUTION OF U AT THE BOUNDARY
    figure
    trisurf(TR, u_boundary)
    colorbar
    axis equal
    hold on  
    title('Solution at the boundary')
    hold off

    % COMPUTE AND PLOT THE SOLUTION OF U IN THE INTERNAL DOMAIN ON A X-Y 
    % PLANE
    postprocessing_domain(u_boundary, dn_u_boundary, ...
        triangles, face_normals, n_elements, TR, ...
        data.enable_postprocessing(2), 'constant')

    toc
    fprintf('-- END POSTPROCESSING --\n\n\n')
end

end



