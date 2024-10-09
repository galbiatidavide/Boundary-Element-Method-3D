function [u_boundary, dn_u_boundary] = get_solution_boundary(u_boundary, ...
    dn_u_boundary, sol, bc_type)

% GET_SOLUTION_BOUNDARY combines the previous assigned boundary conditions
%                       with the solution of the BEM system
% 
% INPUTS
% u_boundary: evaluation of u at the boundary before assigning the BEM
%             system solution at unknown nodes' values
% dn_u_boundary: evaluation of the normal derivative of u at the boundary 
%                before assigning the BEM system solution at unknown nodes'
%                values
% sol: solution of the BEM linear system
% bc_type: array of chars: 'D' stands for Dirichlet and 'N' for Neumann
%          conditions at the element indicated by the index
%
% OUTPUTS
% u_boundary: evaluation of u at the boundary after assigning the BEM
%             system solution at unknown nodes' values
% dn_u_boundary: evaluation of the normal derivative of u at the boundary
%                after assigning the BEM system solution at unknown nodes'
%                values

is_dirichlet = (bc_type=='D');
is_neumann = (bc_type=='N');
if (is_dirichlet + is_neumann ~= ones(length(bc_type),1))
    error('Error in bc_type array')
end


u_boundary(is_neumann) = sol(is_neumann);
dn_u_boundary(is_dirichlet) = -sol(is_dirichlet);


end

