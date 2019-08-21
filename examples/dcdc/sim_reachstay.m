
addpath('../../matlab')
%% load problem & controller
% data saved in .mat:
% - ctree, cindex, cvalue: Tree-structrued controller.
% - U : Two modes.
% - X : Workarea.
% - G: Target area.
% - ts: Sampling time.

% load('data_dcdcReachStay.mat')
load('data_dcdcCoBuchi.mat')
% load('data_dcdcReach.mat')

fm= @dcdc;


%% simulation
x0= [.7; 5.4/5];
% x0= [1.401; 1.181];

Tsim= 100;
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
    uid= find(ctlr(xid(1),:));
%     uid= ctlr_feasible(x, ctree, cindex, cvalue);
    if(isempty(uold))
        u= U(uid(1));
    elseif(find(uid==uold))
        u= uold;
    else
        u= U(uid(1));
    end
    uold= u;
    
    xt= fm(ts,x,u);
    
    % append state for simulation
    usim= cat(1, usim, repmat(u,size(H,2),1));
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
cg= [0.7,0.7,0.7];
cge= [0.4660 0.6740 0.1880];

FS= 16; % fontsize
LW= 1.5; % lineweight

% plot time-state curves
hf1= figure;
subplot(2,1,1)
xsim= interp1(tts, xts, tsim);
plot(tsim, xsim, 'LineWidth',LW)
axis([0, Tsim, min(X(1,1), X(2,1)), max(X(1,2), X(2,2))])
ylabel({'$x_1(t),\; x_2(t)$'}, 'interpreter','latex',...
    'Units','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')

hl= legend({'$x_2(t)$','$x_2(t)$'}, 'Interpreter', 'latex',...
    'Units','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold');

subplot(2,1,2)
plot(tsim, usim,'LineWidth',LW)
axis([0, Tsim, 0, 3])
xlabel({'$Time(s)$'}, 'interpreter','latex',...
    'Units','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')
ylabel({'$u(t)$'}, 'interpreter','latex',...
    'Units','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')

% plot 2-d state space trajectory
hf2= figure;
plot2_boxes(pavings((tag==1)&any(ctlr,2),:), cg, cg, 1);
hold on
rectangle('Position',[X(1,1),X(2,1),X(1,2)-X(1,1),X(2,2)-X(2,1)],...
    'LineWidth',LW, 'LineStyle', '-') % plot the state space
rectangle('Position',[G(1,1),G(2,1),G(1,2)-G(1,1),G(2,2)-G(2,1)],...
    'LineWidth',LW, 'LineStyle', '-', 'EdgeColor', cge) % plot the target space
p= plot(xts(:,1), xts(:,2),'-+','LineWidth',LW);
% p.Color= [39,64,139]/255;
axis([X(1,1) X(1,2) X(2,1) X(2,2)])
plot(x0(1), x0(2), 'o','MarkerFaceColor',cr, 'MarkerEdgeColor',cr)%[176 23 31]/255)[0.8500, 0.3250, 0.0980])
xlabel({'$i_l$'}, 'interpreter','latex',...
    'Units','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')
ylabel({'$V_c$'}, 'interpreter','latex',...
    'Units','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')