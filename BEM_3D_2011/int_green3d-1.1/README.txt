Description
int_green3d_tri and int_green3d_poly are designed to analytically integrate the Green’s function and its gradient over plane faces in 3d. These functions implements the formulas in:

R. Graglia, "Numerical Integration of the linear shape functions times the 3-d Green's function or its gradient on a plane triangle", IEEE Transactions on Antennas and Propagation, vol. 41, no. 10, Oct 1993, pp. 1448--1455

They are intended to be used as kernel functions for the implementation of the boundary element method and other integral methods.

Changelog
Version 1.1: this version includes the function versor.m, missing in the previous version (thanks to Justin Solomon for pointing out!)

Function Reference
* int_green3d_tri.m: integrates 1/r and grad(1/r) analytically over a triangle in 3d
* int_green3d_poly.m: generalizes the previous function to arbitrary plane polygons in 3d
* demo_tri: simple test to demonstrate the syntax of the functions. Results of int_green3d_tri  and int_green3d_poly are compared with the built-in quad2d. The integration domain is a triangle.
* demo_quad: simple test to demonstrate the syntax of the functions. Results of int_green3d_tri  and int_green3d_poly are compared with the built-in quad2d. The integration domain is a square.
* private/eulerangle.m & private/rotation.m: service functions to obtain the rotation matrices

Copyright (c) 2016
Fabio Freschi
Politecnico di Torino   | The University of Queensland
fabio.freschi@polito.it | f.freschi@uq.edu.au