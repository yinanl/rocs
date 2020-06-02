% example from \cite{GirardGM16,RunggerZ16}

addpath('../../matlab')
%% load problem & controller
% data saved in .mat:
% - ctree, cindex, cvalue: Tree-structrued controller.
% - U : Four modes.
% - X : Workarea.
% - G: Target area.
% - ts: Sampling time.
load('data_tpc_reachstay.mat')

vf= @tpc;


%% simulation
x0= [25; 15];
x0= [12; 27];
% x0= [13, 16];
% x0= [26, 21];

Tsim= 1000;

% H= 0:9;
% h= ts/10;
% tfill= h*H';

uold= [];
usim= uold;

% xts= [];
% tts= [];
xsim= [];
tsim= [];
t= 0;
x= x0;
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
    
    [tt,xt]= ode45(@(t, x) vf(t, x, u), [0, ts], x);
    
    % append state for simulation
    usim= cat(1, usim, repmat(u, size(tt,1), 1));
    xsim= cat(1, xsim, xt);
    tsim= cat(1, tsim, t+tt);
%     usim= cat(1, usim, repmat(u,size(H,2),1));
%     tsim= cat(1, tsim, t+tfill);
%     
%     xts= cat(1, xts, x');
%     tts= cat(1, tts, t);
    
    % move to the next step
    x= xt(end,:)';
    t= t + ts;
    
end


%% display
% define color
cr= [0.6350 0.0780 0.1840];
cb= [0    0.4470    0.7410];
co= [0.8500    0.3250    0.0980];
cy= [0.9290    0.6940    0.1250];
cp= [0.4940    0.1840    0.5560];
cg= [0.4660    0.6740    0.1880];

FS= 16; % fontsize
LW= 1.5; % lineweight

% plot time-state curves
hf1= figure;
subplot(2,1,1)
% xsim= interp1(tts, xts, tsim);
plot(tsim, xsim, 'LineWidth',LW)
axis([0, Tsim, min(X(1,1), X(2,1)), max(X(1,2), X(2,2))])
ylabel({'$x_1(t),\; x_2(t)$'}, 'interpreter','latex',...
    'FontSize',FS, 'FontName','Times', 'FontWeight','bold')

hl= legend({'$x_1(t)$','$x_2(t)$'}, 'Interpreter', 'latex',...
    'FontSize',FS, 'FontName','Times', 'FontWeight','bold');

subplot(2,1,2)
plot(tsim, usim,'LineWidth',LW)
axis([0, Tsim, 0, max(U)])
xlabel({'$Time(s)$'}, 'interpreter','latex',...
    'FontSize',FS, 'FontName','Times', 'FontWeight','bold')
ylabel({'$u(t)$'}, 'interpreter','latex',...
    'FontSize',FS, 'FontName','Times', 'FontWeight','bold')

% plot 2-d state space trajectory
hf2= figure;
% plot2_boxes(pavings(tag>0,:), [0.5,0.5,0.5], 'k', 1);
hold on
rectangle('Position',[X(1,1),X(2,1),X(1,2)-X(1,1),X(2,2)-X(2,1)],...
    'LineWidth',LW, 'LineStyle', '-')  % workspace
patch([G(1,1),G(1,2),G(1,2),G(1,1)], [G(2,1),G(2,1),G(2,2),G(2,2)],...
    cg, 'FaceAlpha', 0.7)  % goal area
hobs= patch([xobs(1,1),xobs(1,2),xobs(1,2),xobs(1,1)],...
    [xobs(2,1),xobs(2,1),xobs(2,2),xobs(2,2)], [0.7 0.7 0.7],...
    'FaceAlpha', 0.7);  % avoiding area
hatchfill(hobs, 'cross', 45, 7);

p= plot(xsim(:,1), xsim(:,2),'LineWidth',LW);
% p.Color= [39,64,139]/255;
plot(x0(1), x0(2), 'o','MarkerFaceColor',cr, 'MarkerEdgeColor',cr)%[176 23 31]/255)[0.8500, 0.3250, 0.0980])
axis([X(1,1) X(1,2) X(2,1) X(2,2)])
xlabel({'$x_1$'}, 'interpreter','latex',...
    'FontSize',FS, 'FontName','Times', 'FontWeight','bold')
ylabel({'$x_2$'}, 'interpreter','latex',...
    'FontSize',FS, 'FontName','Times', 'FontWeight','bold')