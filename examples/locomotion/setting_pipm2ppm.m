%% problem data
% % constants
% g = 9.81;
% z1= 1
% z2= 0.7
omega1= 3.1321; % sqrt(g/z1);
omega2= 3.7436; % sqrt(g/z2);

% % keyframe states (determined by contact locations)
x1= [0; 0.5];
x2= [0.6, 1.7];

% % zeta boundaries for R_inter
zeta1_Rinter= [0, 10];
zeta2_Rinter= [-1000000, -0.1];

% % constants for mapping from Euclidean to manifold
Dx= 0.0002;
zeta_0 = 1e-5; 
x1_0= [x1(1)+Dx; sqrt(x1(2)^2+omega1^2*Dx^2)];
x2_0= [x2(1)+Dx; sqrt(x2(2)^2-omega2^2*Dx^2)];

% % discretization parameters
inc_t= 0.02; % [s]
X= [-0.15 0.8; 0.3 1.85]; % state space
eta= [0.004;0.004];
mu= 0.02;
u= 2:mu:4;

% % robust margins [sigma, zeta]
ds1= 0.002;
dz1= 0.05;
q_pipm= [-ds1 ds1; -dz1 dz1];

ds2= 0.06;
dz2= 0.005;
q_ppm= [-ds2 ds2; -dz2 dz2];

% % the margin within robust set 1 (simulation use only)
d_rb1= [0.0005;0.02];