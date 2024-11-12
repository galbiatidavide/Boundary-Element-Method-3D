%% DATA FUNCTION
function [data] = create_data_quadratic(filename)

% CREATE_DATA Used to retrieve information about the mesh, the BCs, 
%             and eventually the exact solutions
%
% INPUTS
% filename: (optional) name(s) of the STL file(s) where to retrieve the mesh
%
% OUTPUTS 
% data: a struct which should be the input of BEM_collocation_constant
%       and BEM_collocation_linear

data = struct();
data.TR = load(filename);

% a cell of strings with the name(s) of the STL file(s) used to 
% generate the mesh; if received as function input, it overwrites the 
% manually inserted one
% data.filename = {'cube_1.stl'};
% if nargin == 1
%     data.filename = filename;
% end

% a lambda function with Dirichlet conditions
data.dirichlet = @dir;

% a lambda function with Neumann conditions
data.neumann = @neu;

% enables the visualization: setting the first element to 0 might reduce 
% computational time, especially with fine meshes. If the first element is
% different from 0, the second element represents the fixed z coordinate
% for the plot inside the domain on a x-y plane with z fixed
data.enable_postprocessing = [1, 0.5];

% the exact solution
data.u_exact = @(x,y,z) x;
data.dn_u_exact = @(x,y,z) - 1.0 .* (x == 0) + ...
                           1.0 .* (x == 1);

% indicates if an iterative solver is used to solve the linear system;
% if the first element is different from 0, the second indicates the 
% the tolerance and the third indicates the number of gmres iterations
data.enable_iterative_solver = [0, 1e-6, 10];
data.bc_type = 'fully Dirichlet';
end

%% BOUNDARY CONDITIONS FUNCTIONS
% Warning: the Dirichlet and Neumann conditions must be set
% complementarily: where a Dirichlet condition is defined you need to set
% the corresponding values of Neumann condition to nan and viceversa

function value = dir(x,y,z)
value = x;
% value(x==0) = nan;
% value(x==1) = nan;
% value(x==0) = -1.0;
% value(x==1) = 1.0;
% value(y==0) = -1.0;
% value(y==1) = 1.0;
% value(z==0) = -1.0;
% value(z==1) = 1.0;
end

function value = neu(x,y,z)
value = nan(size(y));
% value(x==0) = -1.0;
% value(x==1) = 1.0;
% value(y==0) = 0.0;
% value(y==1) = 0.0;
% value(z==0) = 0.0;
% value(z==1) = 0.0;
end
