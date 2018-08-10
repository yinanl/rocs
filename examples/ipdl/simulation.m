
addpath('../../matlab')
%% load spec & controller
% data saved in .mat:
% - Tree-structrued controller: ctree, cindex, cvalue.
% - All input values: U.
% - Workarea: X.
% - Sampling time: ts.
% - Target area: G.
load('data_ipdlCobuchi.mat')


%% simulation
Tsim = 1;
H= 0:9;
h= ts/10;
tfill= h*H';

x0 = [1; 1];

t= 0;
x= x0;
xsim= [];
usim= [];
tsim= [];
while( t <= Tsim)
    
    % compute control
    xid= find(x(1)>=pavings(:,1) & x(1)<=pavings(:,2) & ...
        x(2)>=pavings(:,3) & x(2)<=pavings(:,4)); % direct search
    uid= find(ctlr(xid(1),:));
    u= U(uid(end),:);
    
    % compute next state
    [tt,xx]= ode45(@(t,x) ipdl(t,x,u), [0 ts], x);
    
    % append state for simulation
    usim= cat(1, usim, repmat(u,size(tt,1),1)); % a col
    tsim= cat(1, tsim, t+tt); % a col
    xsim= cat(1, xsim, xx);
        
    % update x, t
    x= xx(end,:)';  % a col
    t = t + tt(end);
end


%% display

% % plot setting
cr= [0.6350 0.0780 0.1840];
cb= [0 0.4470 0.7410];
cy= [0.9290 0.6940 0.1250];
FS= 16; % fontsize
LW= 1.5; % lineweight


% % state and control signals
figure
subplot(2,1,1)

% xsim= interp1(tts, xts, tsim);
plot(tsim, xsim, 'LineWidth',LW)

axis([0, Tsim, min(X(:,1)), max(X(:,2))])

ylabel({'$\theta,\; \dot{\theta}$'}, 'interpreter','latex',...
    'FontUnits','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')

hl= legend({'$\theta(t)$','$\dot{\theta}(t)$'}, 'Interpreter', 'latex',...
    'FontUnits','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold');

subplot(2,1,2)
plot(tsim, usim,'LineWidth',LW)

axis([0, Tsim, U(1), U(end)])

xlabel({'$t(s)$'}, 'interpreter','latex',...
    'FontUnits','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')

ylabel({'$u$'}, 'interpreter','latex',...
    'FontUnits','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')

% print figure
% print -depsc ipdlsim.eps
% % trajectory in 2-d plane
