% Run file for convergence test
close all
clear
clc

addpath('core')
addpath('mesh')
addpath('int_green3d-1.1')
addpath('normals')

data1 = create_data_quadratic('mesh_data_struct_1_giusta.mat');
data05 = create_data_quadratic('mesh_data_struct_0.5_giusta.mat');
data025 = create_data_quadratic('mesh_data_struct_0.25_giusta.mat');
%data0125 = create_data_quadratic('mesh_data_struct_0.125_giusta.mat');

data = [data1, data05, data025];%, data0125];

convergence_test_quadratic(data,2, [20,20,20], 'quadratic')