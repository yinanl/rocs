function dx= aircraftlong(t,x,u)
% AIRCRAFTLONG - the longitudinal dynamics of an aircraft
% SYNTEX
% ---------
% dx/dt = aircraft(t, x, u)
%
% INPUTS
% ---------
%   x:  x(1)-v velocity; x(2)-gamma flight path angle; x(3)-h altitude
%   u:  u(1)-thrust; u(2)-angle of attack
%
% RETURNS
% ---------
%   dx/dt: the derivative of the state x w.r.t. time.
mg= 60000.0*9.81;
mi = 1.0/60000;

c= 1.25+4.2*u(2);
dx= zeros(3,1);

dx(1)= mi*(u(1)*cos(u(2))-(2.7+3.08*c^2)*x(1)^2-mg*sin(x(2)));
dx(2)= (1.0/(60000*x(1)))*(u(1)*sin(u(2))+68.6*c*x(1)^2-mg*cos(x(2)));
dx(3)= x(1)*sin(x(2));

end