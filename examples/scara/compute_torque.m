function torque = compute_torque(u, w, theta)
% Compute the torques for joints by the angular velocities
%
% SYNTEX
% ---------
% torque = compute_torque(w, theta)
% 
% INPUTS
% ---------
%   u(n x 1)     : the accelerates of angular velocities
%	w(n x 1)	 : the vector of angular velocities of joints
%   theta(n x 1) : the vector of joint angles
%
% RETURNS
% ---------
%	torque(n x 1)      : the torque corresponding to the w, theta

m1 = 0.1;
m2 = 0.1;
l1 = 0.15;
l2 = 0.15;
r1 = 0.5*l1;
r2 = 0.5*l2;
I1 = 1.33e-5;
I2 = 1.33e-5;

z1= I1 + I2 + m1*r1^2 + m2*(l1^2+r2^2);
z2= m2*l1*r1;
z3= I2 + m2*r2^2;

M= [z1+2*z2*cos(theta(2)), z3+z2*cos(theta(2));
    z3+z2*cos(theta(2)), z3];
C= [-z2*sin(theta(2))*w(2), -z2*sin(theta(2))*(w(1)+w(2));
    z2*sin(theta(2))*w(1), 0];

torque= M*u + C*w;

end