function [u_boundary,dn_u_boundary,bc_type] = assign_bcs(type, dirichlet,...
    neumann, points)

% ASSIGN_BCS Assigns Dirichlet or Neumann boundary conditions
% 
% INPUTS
% type: type of the boundary conditions to be assigned (either 
%       'fully Dirichlet', 'mixed', or 'fully Deumann')
% dirichlet: handle function to evaluate the Dirichlet bcs with
% neumann: handle function to evaluate the Neumann bcs with
% points: 3D coordinates to evaluate the bcs at
%
% OUTPUTS
% u_boundary: evaluation of u at the boundary
% dn_u_boundary: evaluation of the normal derivative of u at the boundary
% bc_type: type of boundary conditions in each node ('D' or 'N')

if strcmp(type,'fully Dirichlet')
    u_boundary = dirichlet(points(:,1), points(:,2), points(:,3));
    dn_u_boundary = nan(size(u_boundary));
    bc_type = repmat('D', size(points,1), 1);
elseif strcmp(type,'fully Neumann')
    dn_u_boundary = neumann(points(:,1), points(:,2), points(:,3));
    u_boundary = nan(size(dn_u_boundary));
    bc_type = repmat('N', size(points,1), 1);
elseif strcmp(type, 'mixed')
    u_boundary = dirichlet(points(:,1), points(:,2), points(:,3));
    dn_u_boundary = neumann(points(:,1), points(:,2), points(:,3));
    bc_type = repmat('D', size(points,1), 1);
    bc_type(isnan(u_boundary)) = 'N';
else
    error('Inexistent boundary condition type');
end

end

