clear
clc
addpath('../../matlab')
%% load spec & controller
vf= @MG;

%%% Controller data %%%
% - ts: Sampling time.
% - U : All input values.
% - X : Workarea.
% - A(can be empty): obstacles.
% - G(can be empty): Target area.
% - pavings: Tree-structrued controller.
% - tag: indicating if a cell is inside the winning set.
% - ctlr: all valid control inputs for each cell in pavings.

%%% Load from .mat file %%%
% load('controller_2d_caseIReachstay.mat')

%%% Load from .h5 file %%%
ctlrfile= 'controller_2d_caseIReachstay.h5';
ts= h5read(ctlrfile, '/ts');
X= h5read(ctlrfile, '/X')';
U= h5read(ctlrfile, '/U')';
A= permute(h5read(ctlrfile, '/xobs'), [3,2,1]);
G= permute(h5read(ctlrfile, '/G'), [3,2,1]);
pavings= h5read(ctlrfile, '/pavings')';
tag= h5read(ctlrfile, '/tag');
ctlr= h5read(ctlrfile, '/ctlr')';


%% Winning set
% winid= find(any(ctlr,2));
% winset= pavings(winid,:);
winset= pavings(tag==1,:);
wc= [(winset(:,1)+winset(:,2))/2,...
    (winset(:,3)+winset(:,4))/2];
loseid= find(~any(ctlr,2));
loseset= pavings(loseid,:);

figure
% plot(wc(:,1), wc(:,2), '.')
plot2_boxes(winset(:,1:4), [0.5,0.5,0.5], 'k', 1);


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
    if (isempty(uid))
        error("Invalid controller.");
    else
        % select a random u from all valid control values
        pick=randperm(numel(uid));
        u= U(uid(pick(1)),:);
%         % select the minimum value
%         uall= U(uid,:);
%         [val, ind]= min(abs(uall(:,1)));
%         u= uall(ind,:);
%         % select the first/last value
%         u= U(uid(end),:);
    end
    
    % append state for simulation
    usim= cat(1, usim, repmat(u,size(tt,1),1)); % a col
    tsim= cat(1, tsim, t+tt); % a col
    xts= cat(1, xts, x');
    
    % compute next state
    [tt,xx]= ode45(@(t,x) vf(t,x,u), [0 ts], x);
    xsim= cat(1, xsim, xx);
        
    % update x, t
    x= xx(end,:)';  % a col
    t= t + tt(end);
end


%% display
% % plot setting
cr= [0.6350 0.0780 0.1840];
cb= [0 0.4470 0.7410];
% cy= [0.9290 0.6940 0.1250];
cg= [0.7,0.7,0.7];
cge= [0.4660 0.6740 0.1880];
cdg= [34,139,34]/255;
FS= 16; % fontsize
LW= 1.5; % lineweight


% % state and control signals
hf1= figure;
subplot(3,1,1)
% xsim= interp1(tts, xts, tsim);
plot(tsim, xsim, 'LineWidth',LW)
axis([0, Tsim, 0.3, 0.7])
ylabel({'$\Phi,\; \Psi$'}, 'interpreter','latex',...
    'FontSize',FS, 'FontName','Times', 'FontWeight','bold')
hl= legend({'$\Phi(t)$','$\Psi(t)$'}, 'Interpreter', 'latex',...
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
plot2_boxes(winset(:,1:4), [0.7,0.7,0.7], [0.7,0.7,0.7], 1);
% plot2_boxes(pavings((tag==1)&any(ctlr,2),1:4), cg, cg, 1);
hold on

% % plot the state space
% rectangle('Position',[X(1,1),X(2,1),X(1,2)-X(1,1),X(2,2)-X(2,1)],...
%     'LineWidth',LW, 'LineStyle', '-') 
rectangle('Position',[0.40,0.6,0.16,0.08],...
    'LineWidth',LW, 'LineStyle', '-')

% % plot the target set
rectangle('Position',[G(1,1),G(2,1),G(1,2)-G(1,1),G(2,2)-G(2,1)],...
    'Curvature',[1,1],...
    'LineWidth',LW, 'LineStyle', '-', 'EdgeColor', cdg)

% % plot the avoid area
rectangle('Position',[A(1,1),A(2,1),...
    A(1,2)-A(1,1), A(2,2)-A(2,1)],...
    'LineWidth',LW, 'LineStyle', '-', 'EdgeColor', 'k')

% % plot the simulation
pline= plot(xsim(:,1), xsim(:,2),'-','LineWidth',LW);
pline.Color= [39,64,139]/255;
ppts= plot(xts(:,1), xts(:,2), '.', 'LineWidth',LW);
ppts.Color= [189,183,107]/255;

plot(x0(1), x0(2), '.', 'MarkerSize', 20,...
    'MarkerFaceColor',cb, 'MarkerEdgeColor',cb) %[176 23 31]/255)[0.8500, 0.3250, 0.0980])
xlabel({'$\Phi$'}, 'interpreter','latex',...
    'FontSize',FS, 'FontName','Times', 'FontWeight','bold')
ylabel({'$\Psi$'}, 'interpreter','latex',...
    'FontSize',FS, 'FontName','Times', 'FontWeight','bold')

% axis([X(1,1) X(1,2) X(2,1) X(2,2)])
% axis([X(1,1) X(1,2) 0.62 0.68])
axis([0.40 0.56 0.6 0.68])

% print figure
% print -depsc ipdlsim.eps
% % trajectory in 2-d plane
