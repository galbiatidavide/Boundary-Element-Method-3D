function [errors_quadratic_to_squeeze, sizes_quadratic] = run_convergence_test_quadratic()

addpath('core')
addpath('mesh')
addpath('int_green3d-1.1')
addpath('normals')

data1 = create_data_mat('torus_1.mat');
data2 = create_data_mat('torus_2.mat');
data = [data1, data2];

[errors_quadratic_to_squeeze, sizes_quadratic] = convergence_test_quadratic(data,2, [10,10,10], 'quadratic');

end