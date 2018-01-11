addpath('../../utils_grid');
X=[7.3 10; 0 2;-pi pi];
eta=[0.2;0.2;2*pi/35];

[xrv, xl, xtbl] = grid_uniform_generate(X, eta);

T=[9 9.6;0 0.6;-pi pi];
targetSet = grid_uniform_subset(T, xtbl, false, false);