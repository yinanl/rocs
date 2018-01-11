function dx = car(t,x,u)
% ROBOT2D - the dynamics of 2-d robot
% SYNTEX
% ---------
% x_next = robot2d(t, x, u)
%
% INPUTS
% ---------
%   x:  x(1)-x position; x(2)-y position; x(3)-orientation
%   u:  u(1)-velocity; u(2)-steering angle
%
% RETURNS
% ---------
%   dx: the velocity vector at (x,u)

alpha = atan(tan(u(2))/2);

dx = zeros(3, 1);

dx(1) = u(1)*cos(alpha+x(3))/cos(alpha);
dx(2) = u(1)*sin(alpha+x(3))/cos(alpha);
dx(3) = u(1)*tan(u(2));