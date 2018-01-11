% example from \cite{GirardGM16,RunggerZ16}

addpath('../../matlab')
%% load problem & controller
% data saved in .mat:
% - ctree, cindex, cvalue: Tree-structrued controller.
% - U : Two modes.
% - X : Workarea.
% - G: Target area.
% - ts: Sampling time.
load('data_dcdc_spec.mat')
load('data_dcdc_cbox.mat')
load('data_dcdc_ctree.mat')

fm= @dcdc;


%% simulation
x0= [1.2; 1.12];

Tsim= 50;
x= [1.2; 1.12];

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
    uid= ctlr_feasible(x, ctree, cindex, cvalue);
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

FS= 16; % fontsize
LW= 1.5; % lineweight

% plot time-state curves
hf1= figure;
subplot(2,1,1)
xsim= interp1(tts, xts, tsim);
plot(tsim, xsim, 'LineWidth',LW)
axis([0, Tsim, min(G(1,1), G(2,1)), max(G(1,2), G(2,2))])
ylabel({'$x_1(t),\; x_2(t)$'}, 'interpreter','latex',...
    'FontUnits','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')

hl= legend({'$x_2(t)$','$x_2(t)$'}, 'Interpreter', 'latex',...
    'FontUnits','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold');

subplot(2,1,2)
plot(tsim, usim,'LineWidth',LW)
axis([0, Tsim, 0, 3])
xlabel({'$Time(s)$'}, 'interpreter','latex',...
    'FontUnits','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')
ylabel({'$u(t)$'}, 'interpreter','latex',...
    'FontUnits','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')

% plot 2-d state space trajectory
hf2= figure;
plot2_boxes(pavings(tag>0,:), [0.5,0.5,0.5], 'k', 1);
hold on
rectangle('Position',[G(1,1),G(2,1),G(1,2)-G(1,1),G(2,2)-G(2,1)],...
    'LineWidth',LW, 'LineStyle', '-')
p= plot(xts(:,1), xts(:,2),'-+','LineWidth',LW);
% p.Color= [39,64,139]/255;
plot(x0(1), x0(2), 'o','MarkerFaceColor',cr, 'MarkerEdgeColor',cr)%[176 23 31]/255)[0.8500, 0.3250, 0.0980])
axis([G(1,1) G(1,2) G(2,1) G(2,2)])
xlabel({'$x_1$'}, 'interpreter','latex',...
    'FontUnits','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')
ylabel({'$x_2$'}, 'interpreter','latex',...
    'FontUnits','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')