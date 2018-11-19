function [tsim,xsim,usim,dsim]= simulate_trajectory(inc_t,x0,xgrid,ctlr,Q,R,D)
% Simulate the controlled trajectory
%
% INPUT
% ----------
% inc_t: simulation time step.
% x0: the initial state (an 1 x 2 array).
% xgrid: a grid of the local state space.
% ctlr(1 or 2): a structure with fields
%   m: the walking mode.
%   q: the nominal state (or keyframe state).
%   w: the nominal control input.
%   id_goal: indices of the target set.
%   id_init: indices of the initial set.
%   least: the least restrictive controller
%   dt: the extra time stay in the goal set.
% Q,R: coefficients for selecting control values
%   min_{u}(Q*(u-w)^2+R*(u-u_pre)^2)
% D: the disturbance bound.
%
% OUTPUT
% ----------
% tsim, xsim, usim, dsim: arrays of simulation time, state, control and
% disturbance.


% % initial condition
x= x0;
w= ctlr.u1;
t= 0;

dsim= []; xsim= []; usim= []; tsim= [];

% % first semi-step: pipm
ind= get_stateid(x, xgrid);
while (~any(ismember(ind, ctlr.id_goal1)))
    % % compute control input
    uall= find(ctlr.least1(ind,2:end));
    [~,uid]= min(Q*(ctlr.ugrid1(uall)-ctlr.u1).^2 + R*(ctlr.ugrid1(uall)-w).^2);
    w= ctlr.ugrid1(uall(uid));
    if(isempty(w))
        disp('Not within the winning set 1.')
        return;
    end
    
    % % generate a random disturbance
    d= D.*([-1;-1] + 2*rand(2,1));
    % % record
    dsim= [dsim; d'];
    xsim= [xsim; x];
    usim= [usim; w];
    tsim= [tsim; t];
    % % update
    x= vectorfield(ctlr.m1, ctlr.q1, inc_t, x, w, d);
    t= t+inc_t;
    % % check update state
    ind= get_stateid(x, xgrid);
end
t1= t;

% % continue mode 1 until t>tp1 or get out of guardset
while (any(ismember(ind, ctlr.id_goal1)) && t<t1+ctlr.dt1)
    % % generate a random disturbance
    d= D.*([-1;-1] + 2*rand(2,1));
    dsim= [dsim; d'];
    xsim= [xsim; x];
    usim= [usim; w];
    tsim= [tsim; t];
    % % update
    x= vectorfield(ctlr.m1, ctlr.q1, inc_t, x, w, d);
    t= t+inc_t;
    ind= get_stateid(x, xgrid);
end
if (~any(ismember(ind, ctlr.id_goal1)))
    x= xsim(end,:);
    t= t-inc_t;
end

% % second semi-step
w= ctlr.u2;
ind= get_stateid(x, xgrid);
while (~any(ismember(ind, ctlr.id_goal2)))
    % % compute control input
    ind= get_stateid(x, xgrid);
    uall= find(ctlr.least2(ind,2:end));
    [~,uid]= min(Q*(ctlr.ugrid2(uall)-ctlr.u2).^2 + R*(ctlr.ugrid2(uall)-w).^2);
    w= ctlr.ugrid2(uall(uid));
    if(isempty(w))
        disp('Not within the winning set 2.')
        return;
    end
    
    % % generate a random disturbance
    d= D.*([-1;-1] + 2*rand(2,1));
    dsim= [dsim; d'];
    xsim= [xsim; x];
    usim= [usim; w];
    tsim= [tsim; t];
    
    % % update
    x= vectorfield(ctlr.m2, ctlr.q2, inc_t, x, w, d);
    t= t+inc_t;
    %         [sigma, zeta]= manifold('ppm', x, omega2, x2, x2_0, zeta_0);
    ind= get_stateid(x, xgrid);
end
t2= t;

% % continue the last control config for ~0.05s
while (any(ismember(ind, ctlr.id_goal2)) && t<t2+ctlr.dt2)
    % % generate a random disturbance
    d= D.*([-1;-1] + 2*rand(2,1));
    dsim= [dsim; d'];
    xsim= [xsim; x];
    usim= [usim; w];
    tsim= [tsim; t];
    % % update
    x= vectorfield(ctlr.m2, ctlr.q2, inc_t, x, w, d);
    t= t+inc_t;
    %         [sigma, zeta]= manifold('ppm', x, omega2, x2, x2_0, zeta_0);
    ind= get_stateid(x, xgrid);
end