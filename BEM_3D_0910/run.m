close all
clear
clc

% Specifiy if you want constant (0), linear (1) or quadratic (2) elements
element_order = 2;

addpath('core')
addpath('mesh')
addpath('int_green3d-1.1')
addpath('normals')

if element_order == 0
    data = create_data();
    [~, ~, ~] = BEM_collocation_constant(data);
elseif element_order == 1
    data = create_data();
    [~, ~, ~] = BEM_collocation_linear(data);
elseif element_order == 2
    data = create_data_quadratic;
    [~, ~, ~] = BEM_collocation_quadratic(data);
else
    fprintf('Error: Insert a valid element order in the run file: constant (0), linear (1) or quadratic (2) elements\n')
end
