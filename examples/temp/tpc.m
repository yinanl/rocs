function dx= tpc(t, x, u)
% TPC - the dynamics of a temperature controller
% SYNTEX
% ---------
% x_next = tpc(t, x, u)
%
% INPUTS
% ---------
%   t: time
%   x:  x(1)-room temperature; x(2)-heater temperature
%   u:  four modes (u= off(1), heating(2), cooling(3), on(4))
%
% RETURNS
% ---------
%   dx: the next step state
r1= 0.002;
r2= 0.1;
c= 16;

switch u
    case 1
        dx= [-r1*(x(1)-c); 0];
    case 2
        dx= [-r1*(x(1)-x(2)); r2];
    case 3
        dx= [-r1*(x(1)-x(2)); -r2];
    case 4
        dx= [-r1*(x(1)-x(2)); 0];
    otherwise
        error('Wrong mode.')
end
