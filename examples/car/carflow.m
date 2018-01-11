function x_next=carflow(t,x,u)
% SOLVE_ROBOT2D - the solution of 2-d robot.
% SYNTEX
% ---------
% x_next = solve_robot2d(t, x, u)
%
% INPUTS
% ---------
%   x:  x(:,1)-x position; x(:,2)-y position; x(:,3)-orientation
%   u:  u(1)-velocity; u(2)-steering angle
%   t:  time duration
%
% RETURNS
% ---------
%   x_next: the next state after t under control u

x_next= x;

alpha = atan(tan(u(2))/2);
if(u(1)==0)
    x_next = x;
elseif(u(2)==0)
    x_next(:,3) = x(:,3);
    x_next(:,1) = u(1)*cos(x(:,3))*t + x(:,1);
    x_next(:,2) = u(1)*sin(x(:,3))*t + x(:,2);
else
%     x_next(:,3) = u(1)*tan(u(2))*t + x(:,3);
%     x_next(:,1) = (sin(alpha+x_next(:,3))-sin(alpha+x(:,3)))/cos(alpha)/tan(u(2)) + x(:,1);
%     x_next(:,2) = -(cos(alpha+x_next(:,3))-cos(alpha+x(:,3)))/cos(alpha)/tan(u(2)) + x(:,2);
    at= alpha + u(1)*tan(u(2))*t/2;
    r = 2*sin(u(1)*tan(u(2))*t/2)/cos(alpha)/tan(u(2));
    
    x_next(:,1) = r*cos(at+x(:,3)) + x(:,1);
    x_next(:,2) = r*sin(at+x(:,3)) + x(:,2);
    x_next(:,3) = u(1)*tan(u(2))*t + x(:,3);
end