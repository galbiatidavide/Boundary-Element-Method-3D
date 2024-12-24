%% CONVERGENCE TEST FUNCTION
function [relative_errors, mesh_sizes] = convergence_test_linear_constant(data, norm_type, n_sample_points, which_bem)

% CONVERGENCE_TEST Used to perform a convergence test (with either constant
%                  and linear elements)
%
% INPUTS
% data: the struct which is the output of the function create_data() in 
%       create_data.m
% norm_type: a value which allows the user to choose the norm to evaluate
%            the convergence with
% n_sample_points: number of points per dimension to generate a grid which
%                  contains the triangulation; will be used to keep the 
%                  inside points only
% which_bem: 'constant' or 'linear'

for k=1:length(data)
    current_data = data(k);
    if (strcmp(which_bem, 'constant'))
        [u_boundary_nodal, dn_u_boundary_nodal, TR] = ...
            BEM_collocation_constant(current_data);
    elseif (strcmp(which_bem, 'linear'))
        [u_boundary_nodal, dn_u_boundary_nodal, TR] = ...
            BEM_collocation_linear(current_data);
    else
        error('Inexistent kind of BEM chosen')
    end
    % storing the relative error and the mesh size
    for n=1:length(norm_type)
        [relative_errors(n,k,:), mesh_sizes(k)] = ...
            compute_relative_error(norm_type(n), u_boundary_nodal, ...
            dn_u_boundary_nodal, current_data.u_exact, current_data.dn_u_exact, ...
            TR, n_sample_points, which_bem);
    end
end

% display mesh sizes
fprintf('\nMESH SIZES');
disp(mesh_sizes);
fprintf('RELATIVE ERRORS\n\n');

for n=1:length(norm_type)
    % display relative errors w.r.t. norm
    fprintf('Norm: %d\n', norm_type(n));
    disp(squeeze(relative_errors(n,:,:)));
    
    % set figure to plot
    figure
    sgtitle(['CONVERGENCE TEST: NORM ', num2str(norm_type(n))], 'Fontsize', 32)
    
    % plotting the convergence results (normal derivative at the boundary)
    subplot(1,3,1)
    mesh_sizes05 = mesh_sizes.^0.5;
    mesh_sizes2 = mesh_sizes.^2;
    loglog(mesh_sizes, relative_errors(n,:,1) / relative_errors(n,1,1), 'r-', ...
        'LineWidth', 2);
    grid on
    hold on
    loglog(mesh_sizes, mesh_sizes05 / mesh_sizes05(1), 'k--');
    loglog(mesh_sizes, mesh_sizes / mesh_sizes(1), 'k-');
    loglog(mesh_sizes, mesh_sizes2 / mesh_sizes2(1), 'k-*');
    legend('relative errors', '1/2 order', '1 order', ...
        '2 order', 'Location', 'Southeast', 'Fontsize', 16)
    title('Normal derivative of the solution at the boundary', 'Fontsize', 18)

    % plotting the convergence results (solution at the boundary)
    subplot(1,3,2)
    loglog(mesh_sizes, relative_errors(n,:,2) / relative_errors(n,1,2), 'r-', ...
        'LineWidth', 2);
    grid on
    hold on
    loglog(mesh_sizes, mesh_sizes05 / mesh_sizes05(1), 'k--');
    loglog(mesh_sizes, mesh_sizes / mesh_sizes(1), 'k-');
    loglog(mesh_sizes, mesh_sizes2 / mesh_sizes2(1), 'k-*');
    legend('relative errors', '1/2 order', '1 order', ...
        '2 order',  'Location', 'Southeast', 'Fontsize', 16)
    title('Solution at the boundary', 'Fontsize', 18)

    % plotting the convergence results (solution at the internal points)
    subplot(1,3,3)
    loglog(mesh_sizes, relative_errors(n,:,3) / relative_errors(n,1,3), ...
        'r-', 'LineWidth', 2);
    grid on
    hold on
    loglog(mesh_sizes, mesh_sizes05 / mesh_sizes05(1), 'k--');
    loglog(mesh_sizes, mesh_sizes / mesh_sizes(1), 'k-');
    loglog(mesh_sizes, mesh_sizes2 / mesh_sizes2(1), 'k-*');
    legend('relative errors', '1/2 order', '1 order', ...
        '2 order',  'Location', 'Southeast', 'Fontsize', 16)
    title('Solution at the internal points', 'Fontsize', 18)
end

warning('The reason why some of the relative errors are null may be the fact that you imposed the values via bcs')

end

%% FUNCTION TO COMPUTE THE RELATIVE ERROR WITH
function [relative_errors, mesh_size] = compute_relative_error(norm_type, ...
    u_boundary_nodal, dn_u_boundary_nodal, u_exact, dn_u_exact, TR, ...
    n_sample_points, which_bem)
% COMPUTE_RELATIVE_ERROR Used to compute the relative error 
%
% INPUTS
% norm_type: a value which allows the user to choose the norm to evaluate
%            the convergence with
% u_boundary_nodal: the solution on the boundary
% dn_u_boundary_nodal: the normal derivative of the solution on the boundary
% u_exact: the exact solution (lambda function)
% dn_u_exact: the normal derivative of the exact solution (lambda function)
% TR: the triangulation
% n_sample_points: number of points per dimension to generate a grid which
%                  contains the triangulation; will be used to keep the 
%                  inside points only
% which_bem: 'constant' or 'linear'
%
% OUTPUTS:
% relative_errors: the relative errors (1: of the normal derivative at the 
%                  boundary, 2: of the solution at the boundary, 3: of the 
%                  solution at internal points), computed with the desired
%                  norm (specified by norm_type)
% mesh_size: the biggest edge of the triangulation

% extracting information from the triangulation
points = TR.Points;
n_elements = length(TR.ConnectivityList);
triangles = zeros(n_elements, 3, 3);
mesh_size = -Inf;
for i = 1:n_elements
    % extracting 3D coordinates of the vertices of each flat triangle
    v1 = TR.Points(TR.ConnectivityList(i,1),:);
    v2 = TR.Points(TR.ConnectivityList(i,2),:);
    v3 = TR.Points(TR.ConnectivityList(i,3),:);
    d1 = sqrt(sum((v2 - v3).^2));
    d2 = sqrt(sum((v1 - v3).^2));
    d3 = sqrt(sum((v2 - v1).^2));
    triangles(i,:,:) = [v1; ...
                        v2; ...
                        v3];
    mesh_size = max(mesh_size, max(d1, max(d2, d3)));
end
% extract baricenter of each element
baricenters = mean(triangles, 2);
baricenters = squeeze(baricenters);
%computing face normals
face_normals = faceNormal(TR);

% generating an alphaShape
SH = alphaShape(points(:,1), points(:,2), points(:,3));

% generating a grid which surely contains all the points inside the
% triangulation
epsilon = 1e-1*max(abs(points)); % to guarantee that you don't keep points 
                                 % on the boundary
coordinate_extremes = [min(points,[],1) + epsilon; ...
                       max(points,[],1) - epsilon];
x_coord = linspace(coordinate_extremes(1,1), coordinate_extremes(2,1),...
 n_sample_points(1));
y_coord = linspace(coordinate_extremes(1,2), coordinate_extremes(2,2),...
 n_sample_points(2));
z_coord = linspace(coordinate_extremes(1,3), coordinate_extremes(2,3),...
 n_sample_points(3));
[XX, YY, ZZ] = meshgrid(x_coord, y_coord, z_coord);
XX = reshape(XX, ...
 n_sample_points(1) * n_sample_points(2) * n_sample_points(3), 1, 1);
YY = reshape(YY, ...
 n_sample_points(1) * n_sample_points(2) * n_sample_points(3), 1, 1);
ZZ = reshape(ZZ, ...
 n_sample_points(1) * n_sample_points(2) * n_sample_points(3), 1, 1);

% keep only the points which are inside the triangulation
is_in_shape = inShape(SH, XX, YY, ZZ);
XX_in = XX(is_in_shape);
YY_in = YY(is_in_shape);
ZZ_in = ZZ(is_in_shape);
points_in = [XX_in, YY_in, ZZ_in];

% evaluating the exact solution
u_exact_eval = u_exact(XX_in, YY_in, ZZ_in);
u_exact_eval_boundary = u_exact(baricenters(:,1), baricenters(:,2),...
    baricenters(:,3));
dn_u_exact_eval_boundary = dn_u_exact(baricenters(:,1), baricenters(:,2),...
    baricenters(:,3));

% computing the solution with BEM
if (strcmp(which_bem, 'constant'))
    u_exact_eval_boundary = u_exact(baricenters(:,1), baricenters(:,2),...
        baricenters(:,3));
    dn_u_exact_eval_boundary = dn_u_exact(baricenters(:,1), baricenters(:,2),...
        baricenters(:,3));
    solution = get_solution_domain_constant(u_boundary_nodal, ...
        dn_u_boundary_nodal, triangles, face_normals, points_in, n_elements);
elseif (strcmp(which_bem, 'linear'))
     u_exact_eval_boundary = u_exact(points(:,1), points(:,2),...
        points(:,3));
    dn_u_exact_eval_boundary = dn_u_exact(points(:,1), points(:,2),...
        points(:,3));
    solution = get_solution_domain_linear(u_boundary_nodal, ...
        dn_u_boundary_nodal, triangles, face_normals, points_in, ...
        n_elements, TR);
else
        error('Inexistent kind of BEM chosen')
end

% storing the relative errors (1: of the normal derivative at the boundary
%                              2: of the solution at the boundary
%                              3: of the solution at internal points)
relative_errors(1) = norm(dn_u_exact_eval_boundary - dn_u_boundary_nodal,...
    norm_type) / norm(dn_u_exact_eval_boundary, norm_type);
relative_errors(2) = norm(u_exact_eval_boundary - u_boundary_nodal, ...
    norm_type) / norm(u_exact_eval_boundary, norm_type);
relative_errors(3) = norm(u_exact_eval - solution, norm_type) /...
    norm(u_exact_eval, norm_type);

end

