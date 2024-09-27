# COLLOCATION BEM FOR 3D LAPLACE EQUATION

This code contains the implementation of the BEM with constant and linear 
elements and collocation for 3D Laplace equation. Only flat triangular meshes 
are considered.

All the meshes must be inside the folder *mesh* for the correct functioning 
of the code.

This code depends on the the kernel functions provided in the folder
*int_green3d-1.1* (external library).

In the folder *core* you can find the core functions developed for the 
implementation of the Boundary Element Method. 

In the main folder are located the functions you need to call from command line: this is where you need to navigate to, if you want to follow the following instructions.

### BASIC USAGE (From command line)

1. `clear all`

2. Modify *create_data.m*

3. `data = create_data()`

4. `[~, ~, ~] = BEM_collocation_*(data)` (where * can be `constant` or `linear`)

### USAGE FOR CONVERGENCE TEST (From command line)

1. `clear all`

2. Modify *create_data.m*

3. `data = create_data()`

4. `norm_type = ...`

5. `n_sample_points = [..., ..., ...]`

6. `convergence_test(data, norm_type, n_sample_points, *)` (where * can be the string `'constant'` or `'linear'`)