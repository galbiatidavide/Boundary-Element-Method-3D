% This file is used to load the mesh generated with GMSH in a .msh file.
% Since MATLAB does not support quadratic elements this mesh is generated
% in python with GMSH and then the connectivity, points and triangles data
% are loaded separately as .mat files.

% Load the .mat file into MATLAB
TR = load('mesh_data_struct.mat');