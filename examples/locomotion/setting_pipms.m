%%%
% Problem data: PIPM to PIPM
%%%

% % constants
% g = 9.81;
% z1= 1;
% z2= 1.2;
omega1= 3.1321; % sqrt(g/z1);
omega2= 2.8592; % sqrt(g/z2);

% % keyframe states (determined by contact locations)
x1= [0; 0.5];
x2= [0.5; 0.6];

% % asymptotes y= (+/-)k(x-xcontact): they set the safety constraints
k1= omega1;
k2= omega2;

% % constants for mapping from Euclidean to manifold
zeta_0 = 1e-5; 
x1_0= [0.00005; 0.5];
x2_0= [0.50006; 0.6];


% % discretization parameters
inc_t = 0.02; % [s]

X= [-0.15 0.7; 0.2 1.1]; % state space
eta= [0.005;0.005];

mu= 0.02; % control space
u= 2:mu:4;

rbset1= [-0.002 0.002; -0.05, 0.05]; % robust sets: [sigma, zeta]
rbset2= [-0.006 0.006; -0.05, 0.05];

d_rb1= [0.0005;0.001];  % the margin within robust set 1 (simulation use only)