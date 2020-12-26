function dx = scara(t, x, u)
% Dynamics of 2 link SCARA robot
% SYNTEX
% ---------
% x_derivative = scara(t, x, u)
% 
% INPUTS
% ---------
%	x		    : states [theta1; theta2; theta1_dot; theta2_dot]
%   u			: control torques [tau1; tau2]
%
% RETURNS
% ---------
%	dx          : the derivative of states.

m1 = 0.1;
m2 = 0.1;
l1 = 0.15;
l2 = 0.15;
r1 = 0.5*l1;
r2 = 0.5*l2;
I1 = 1.33e-5;
I2 = 1.33e-5;

z1= I1 + I2 + m1*r1^2 + m2*(l1^2+r2^2);
z2= m2*l1*r2;
z3= I2 + m2*r2^2;
% z4= I1 + m1*r1^2 + m2*l1^2; % z4= z1-z3

detM= z3*(z1-z3) - z2^2*cos(x(2))^2; % a1 > 0
a= z2 * sin(x(2)) * (x(4) + 2*x(3)) * x(4);
c= z2 * x(3).^2 * sin(x(2)) - u(2);
b= z2*cos(x(2));

dx = [x(3); x(4); ...
    (z3*u(1) + z3*a + (z3+b)*c)/detM;...
    ((z1+2*b)*(-c)-(z3+b)*(u(1)+a))/detM];
end