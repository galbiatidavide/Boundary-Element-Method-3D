function [errors_linear_to_squeeze, sizes] = run_convergence_test_linear_constant()



addpath('core')
addpath('mesh')
addpath('int_green3d-1.1')
addpath('normals')

data1 = create_data_stl('torus_linear_1.stl');
data2 = create_data_stl('torus_linear_2.stl');

data = [data1, data2];

[errors_linear_to_squeeze, sizes] = convergence_test_linear_constant(data,2, [10,10,10], 'linear');

end