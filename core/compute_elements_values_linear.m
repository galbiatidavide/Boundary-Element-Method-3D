function [element_values] = compute_elements_values_linear(TR, ...
    nodal_values)

% COMPUTE_ELEMENTS_VALUES_LINEAR computes the values of the solution for
%                                each element averaging the nodal values at
%                                the vertexes which define the triangle
%
% INPUT
% TR: the triangulation
% nodal_values: the nodal values of the solution
%
% OUTPUT
% element_values: the values of the solution at each element

element_values = mean(nodal_values(TR.ConnectivityList),2);

end

