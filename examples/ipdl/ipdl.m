function dx=ipdl(t, x, u)
% the vector field of inverted pendulum
%   @states:     x: 2-dimensional, x(1)-angular position; x(2)-rate of x(2)
%   @controls:   u: normalized torque
%   

% pendulum system parameters
M= 0.5; % mass of the cart
m= 0.2; % mass of the pendulum
b= 0.1; % coefficient of friction of cart
J= 0.006; % inertia
g= 9.8; % gravitational constant
l= 0.3; % pendulum length

% the vector field
dx= [x(2);((m*g*l)*sin(x(1)) -b*x(2) +l*cos(x(1))*u)/(J+m*l^2)];