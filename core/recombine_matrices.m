function [A, rhs] = recombine_matrices(G,  H, u_boundary, dn_u_boundary,...
    bc_type)

% RECOMBINE_MATRICES defines the system matrix and the rhs
% 
% INPUTS
% G: the matrix which corresponds to the single layer potential
% H: the matrix which corresponds to the double layer potential
% u_boundary: evaluation of u at the boundary
% dn_u_boundary: evaluation of the normal derivative of u at the boundary
% bc_type: array of chars: 'D' stands for Dirichlet and 'N' for Neumann
%          conditions at the element indicated by the index
%
% OUTPUTS
% A: the system matrix
% rhs: the right hand side

is_dirichlet = (bc_type=='D');
is_neumann = (bc_type=='N');
if (is_dirichlet + is_neumann ~= ones(length(bc_type),1))
    error('Error in bc_type array')
end
% n = size(H,1);
% m = size(H, 2);
% if(n == m)
%     for i = 1:n
%         % Calculate the sum of the row excluding the diagonal element
%         rowSum = sum(G(i, :));
% 
%         % Get the diagonal element
% 
%         % Print the sum and diagonal element for the current row
%         fprintf('Row %d: Sum = %d\n', i, rowSum);
%     end
% end
A = zeros(size(H));
A(:,is_dirichlet) = G(:,is_dirichlet);
A(:,is_neumann) = H(:,is_neumann);
assigned = zeros(size(bc_type));
assigned(is_dirichlet) = u_boundary(is_dirichlet);
rhs = -H * assigned;
assigned = zeros(size(bc_type));
assigned(is_neumann) = dn_u_boundary(is_neumann);
rhs = rhs + G * assigned;
% A = G;
% rhs = -H*u_boundary(is_dirichlet);
% A = H;
% rhs = G * dn_u_boundary(is_neumann);


end

