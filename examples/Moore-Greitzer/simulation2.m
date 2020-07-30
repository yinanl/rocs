
addpath('../../matlab')
%% load spec & controller
% data saved in .mat:
% - Tree-structrued controller: ctree, cindex, cvalue.
% - All input values: U.
% - Workarea: X.
% - Sampling time: ts.
% - Target area: G.
global a B H W lc cx cy aH H2 W2
a = 1/3.5;
B = 2.0;
H = 0.18;
W = 0.25;
lc = 8.0;
cx = 1.0/lc;
cy = 1.0/(4*lc*B*B);
aH = a+H;
H2 = H/(2.0*W*W*W);
W2 = 3*W*W;

% load('data_caseII_Cobuchi2.mat')
load('data_caseIIReachstay2.mat')
% load('data_engineReach2.mat')

winid= find(any(ctlr,2));
winset= pavings(winid,:);
wc= [(winset(:,1)+winset(:,2))/2,...
    (winset(:,3)+winset(:,4))/2];
loseid= find(~any(ctlr,2));
loseset= pavings(loseid,:);


%% simulation
Tsim = 10;
% ts= 1
% H= 0:9;
% h= ts/10;
% tfill= h*H';

x0= [0.5343; 0.6553];
muold= 0.66;

t= 0;
x= x0;
xsim= [];
xts= x0';
usim= [];
tsim= [];

while( t <= Tsim)
    % compute control
    xid= find(x(1)>=pavings(:,1) & x(1)<=pavings(:,2) & ...
        x(2)>=pavings(:,3) & x(2)<=pavings(:,4)); % direct search
    uid= find(ctlr(xid(1),:));
    
%     % select a random u from all valid control values
%     pick=randperm(numel(uid));
%     u= U(uid(pick(1)),:);
%     % select the minimum value
%     uall= U(uid,:);
%     [val, ind]= min(abs(uall(:,1)));
%     u= uall(ind,:);
%     % select the first/last value
%     u= U(uid(end),:);

    % select the minimum mu (u(2)) change
    uall= U(uid,:);
    [val, ind]= min(abs(uall(:,2)-muold));
    u= uall(ind,:);

    % compute next state
    [tt,xx]= ode45(@(t,x) mg(t,x,u), [0 ts], x);
    
    % append state for simulation
    usim= cat(1, usim, repmat(u,size(tt,1),1)); % a col
    tsim= cat(1, tsim, t+tt); % a col
    xsim= cat(1, xsim, xx);
    xts= cat(1, xts, xx(end,:));
        
    % update x, t
    x= xx(end,:)';  % a col
    t= t + tt(end);
    muold= u(2);
end


%% display
% % plot setting
cr= [0.6350 0.0780 0.1840];
cb= [0 0.4470 0.7410];
% cy= [0.9290 0.6940 0.1250];
cg= [0.7,0.7,0.7];
cge= [0.4660 0.6740 0.1880];
FS= 16; % fontsize
LW= 1.5; % lineweight


% % state and control signals
hf1= figure;
subplot(3,1,1)
% xsim= interp1(tts, xts, tsim);
plot(tsim, xsim, 'LineWidth',LW)
axis([0, Tsim, 0.3, 0.7])
ylabel({'$x,\; y$'}, 'interpreter','latex',...
    'FontSize',FS, 'FontName','Times', 'FontWeight','bold')
hl= legend({'$x(t)$','$y(t)$'}, 'Interpreter', 'latex',...
    'FontSize',FS, 'FontName','Times', 'FontWeight','bold');
subplot(3,1,2)
plot(tsim, usim(:,1),'LineWidth',LW)
axis([0, Tsim, min(U(:,1)), max(U(:,1))])
xlabel({'$t(s)$'}, 'interpreter','latex',...
    'FontSize',FS, 'FontName','Times', 'FontWeight','bold')
ylabel({'$u$'}, 'interpreter','latex',...
    'FontSize',FS, 'FontName','Times', 'FontWeight','bold')
subplot(3,1,3)
plot(tsim, usim(:,2),'LineWidth',LW)
axis([0, Tsim, min(U(:,2)), max(U(:,2))])
xlabel({'$t(s)$'}, 'interpreter','latex',...
    'FontSize',FS, 'FontName','Times', 'FontWeight','bold')
ylabel({'$\mu$'}, 'interpreter','latex',...
    'FontSize',FS, 'FontName','Times', 'FontWeight','bold')

% plot 2-d state space trajectory
hf2= figure;
ax= gca;
% plot(wc(:,1), wc(:,2), '.')
% plot2_boxes(winset(:,1:4), [0.5,0.5,0.5], 'k', 1);
% plot2_boxes(pavings((tag==1)&any(ctlr,2),1:4), cg, cg, 1);
hold on
rectangle('Position',[X(1,1),X(2,1),X(1,2)-X(1,1),X(2,2)-X(2,1)],...
    'LineWidth',LW, 'LineStyle', '-') % plot the state space
rectangle('Position',[G(1,1),G(2,1),G(1,2)-G(1,1),G(2,2)-G(2,1)],...
    'LineWidth',LW, 'LineStyle', '-', 'EdgeColor', cge) % plot the target space
rectangle('Position',[xobs(1,1),xobs(2,1),xobs(1,2)-xobs(1,1),xobs(2,2)-xobs(2,1)],...
    'LineWidth',LW, 'LineStyle', '-', 'EdgeColor', 'k')
% p= plot(xsim(:,1), xsim(:,2),'-+','LineWidth',LW);
p= plot(xts(:,1), xts(:,2),'-+','LineWidth',LW);
p.Color= [39,64,139]/255;
axis([X(1,1) X(1,2) X(2,1) X(2,2)])
plot(x0(1), x0(2), '.','MarkerFaceColor',cr, 'MarkerEdgeColor',cr)%[176 23 31]/255)[0.8500, 0.3250, 0.0980])
xlabel({'$x$'}, 'interpreter','latex',...
    'FontSize',FS, 'FontName','Times', 'FontWeight','bold')
ylabel({'$y$'}, 'interpreter','latex',...
    'FontSize',FS, 'FontName','Times', 'FontWeight','bold')

% print figure
% print -depsc ipdlsim.eps
% % trajectory in 2-d plane
