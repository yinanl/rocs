
addpath('../../matlab')
%% load problem & controller
% data saved in .mat:
% - ctree, cindex, cvalue: Tree-structrued controller.
% - U : Two modes.
% - X : Workarea.
% - G: Target area.
% - ts: Sampling time.
load('data_ddeInv.mat')

func= @(x) [0.1*x(1)+0.1*x(3)-0.2*x(4); 0.4*x(1)+0.1*x(2)+0.4*x(3)+0.5*x(4);x(1);x(2)];


%% simulation
x0= [-0.8; -0.8];

Tsim= 50;
x= x0;

H= 0:9;
h= ts/10;
tfill= h*H';

uold= [];
tsim= [];
usim= uold;

xts= [];
tts= [];

t= 0;
while(t<Tsim)
    
    % compute control input
    xid= find(x(1)>=pavings(:,1) & x(1)<=pavings(:,2) & ...
        x(2)>=pavings(:,3) & x(2)<=pavings(:,4)); % direct search
    
    xt= func(x);
    
    % append state for simulation
    tsim= cat(1, tsim, t+tfill);
    
    xts= cat(1, xts, x');
    tts= cat(1, tts, t);
    
    % move to the next step
    x= xt;
    t= t + ts;
    
end


%% display
% define color
cr= [0.6350 0.0780 0.1840];
cb= [0 0.4470 0.7410];
cy= [0.9290 0.6940 0.1250];

FS= 16; % fontsize
LW= 1.5; % lineweight

% plot time-state curves
hf1= figure;
xsim= interp1(tts, xts, tsim);
plot(tsim, xsim, 'LineWidth',LW)
axis([0, Tsim, min(X(1,1), X(2,1)), max(X(1,2), X(2,2))])
ylabel({'$x_1(t),\; x_2(t)$'}, 'interpreter','latex',...
    'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')

hl= legend({'$x_2(t)$','$x_2(t)$'}, 'Interpreter', 'latex',...
    'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold');


% plot 2-d state space trajectory
hf2= figure;
plot2_boxes(pavings(tag>0,:), [0.5,0.5,0.5], 'k', 1);
hold on
rectangle('Position',[X(1,1),X(2,1),X(1,2)-X(1,1),X(2,2)-X(2,1)],...
    'LineWidth',LW, 'LineStyle', '-')
p= plot(xts(:,1), xts(:,2),'-+','LineWidth',LW);
% p.Color= [39,64,139]/255;
plot(x0(1), x0(2), 'o','MarkerFaceColor',cr, 'MarkerEdgeColor',cr)%[176 23 31]/255)[0.8500, 0.3250, 0.0980])
axis([X(1,1) X(1,2) X(2,1) X(2,2)])
xlabel({'$x_1$'}, 'interpreter','latex',...
    'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')
ylabel({'$x_2$'}, 'interpreter','latex',...
    'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')
