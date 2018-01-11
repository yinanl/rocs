function dx= dcdc(t, x, u)
% DCDC - the dynamics of a DCDC converter
% SYNTEX
% ---------
% x_next = dcdc(t, x, u)
%
% INPUTS
% ---------
%   ts: samping time
%   x:  x(1)-capacitor voltage; x(2)-inductor current
%   u:  two modes (u= 1 or 2)
%
% RETURNS
% ---------
%   dx: the next step state
xc= 70;
xl= 3;
rc= 0.005;
rl= 0.05;
r0= 1;
vs= 1;

A1= [-rl/xl 0;
    0 -1/(xc*(rc+r0))];
b1= [vs/xl; 0];

A2= [(-1/xl)*(rl + r0*rc/(r0+rc)) (-1/xl)*(r0/(r0+rc));
    (1/xc)*(r0/(r0+rc)) (-1/xc)*(1/(r0+rc))];
b2= b1;


I= eye(2);

if(u==1)
    dx= expm(A1*t)*x + A1\I*(expm(A1*t)-I)*b1;
elseif(u==2)
    dx= expm(A2*t)*x + A2\I*(expm(A2*t)-I)*b2;
else
    error('Wrong mode.')
end
% 
% O= [1.15 1.55; 5.45/5 5.85/5];
% x0= [1.2; 1.12];

% O= [1.3 1.7; 1.1 1.3];
% x0= [1.6; 1.25];

% \cite{FribourgGPM16}
% O= [1.55 2.15; 1.0 1.4];