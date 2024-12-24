%% DATA FUNCTION
function [data] = create_data_mat(filename)

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

% a lambda function with Dirichlet conditions
data.dirichlet = @dir;

% a lambda function with Neumann conditions
data.neumann = @neu;

% enables the visualization: setting the first element to 0 might reduce 
% computational time, especially with fine meshes. If the first element is
% different from 0, the second element represents the fixed z coordinate
% for the plot inside the domain on a x-y plane with z fixed
data.enable_postprocessing = [1, 0.5];

% TEST CASE 1
% data.u_exact = @(x,y,z) x;
% data.dn_u_exact = @(x,y,z) - 1.0 .* (x == 0) + 1.0 .* (x == 1);

% TEST CASE 2
% data.u_exact = @(x,y,z) z.*cos(y).*exp(x);
% data.dn_u_exact = @(x,y,z) exp(x).*cos(y) .* (z == 1) - exp(x).*cos(y) .* (z == 0) + ...
%     -z.*cos(y).*(x == 0) + z.*cos(y).*exp(1).*(x == 1) + -z.*exp(x).*(y == 1);

R = 2;
data.u_exact = @(x,y,z) exp(x).*sin(z) + exp(z).*cos(y);
data.dn_u_exact = @(x,y,z) (1./(2*sqrt(R^2 - 2*R*sqrt(x.^2+y.^2) + x.^2 + y.^2 + z.^2))).*(exp(x).*sin(z).*(2*x.*(-R+sqrt(x.^2+y.^2)))./(sqrt(x.^2+y.^2)) - ...
    exp(z).*sin(y).*(2*y.*(-R+sqrt(x.^2+y.^2)))./(sqrt(x.^2+y.^2)) + (exp(z).*cos(y) + exp(x).*cos(z)).*z*2);

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

%TEST CASE 1
% function value = dir(x,y,z)
% value = x;
% % value(x==0) = nan;
% % value(x==1) = nan;
% % value(x==0) = -1.0;
% % value(x==1) = 1.0;
% % value(y==0) = -1.0;
% % value(y==1) = 1.0;
% % value(z==0) = -1.0;
% % value(z==1) = 1.0;
% end
% 
% function value = neu(x,y,z)
% value = nan(size(y));
% % value(x==0) = -1.0;
% % value(x==1) = 1.0;
% % value(y==0) = 0.0;
% % value(y==1) = 0.0;
% % value(z==0) = 0.0;
% % value(z==1) = 0.0;
% end

%TEST CASE 3
function value = dir(x,y,z)
value = exp(x).*sin(z) + exp(z).*cos(y);
end

function value = neu(x,y,z)
value = nan(size(y));
end