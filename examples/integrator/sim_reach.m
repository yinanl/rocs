
addpath('../../matlab')
%% load problem & controller
% data saved in .mat:
% - cpavings[nBox x 2n]: a column of intervals that represent the state space
% - ctlr[nBox x nU]: a column of permissible controls (marked 1)
% - tag[nBox x 1]: a column of tags. The tag of the winning set is 1.
% - U : Two modes.
% - X : Workarea.
% - G: Target area.
% - ts: Sampling time.
load('data_reach.mat')

% define the discrete-time dynamics of a double integrator
Dintegrator= @(ts,x,u) [x(1) + ts*x(2) + 0.5*ts^2*u(1); x(2) + ts*u(1)];


%% simulation
x0= [1.2; 1.12];

Tsim= 10;
x= [1.2; 1.12];

H= 0:9;
h= ts/10;
tfill= h*H';

uold= [];
tsim= [];
usim= uold;

xts= x';
tts= 0;

t= 0;
while(x(1)>G(1,2) || x(1)<G(1,1) ||...
        x(2)>G(2,2) || x(2)<G(2,1))
    
    % compute control input
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
    
    xt= Dintegrator(ts,x,u);
    
    % append state for simulation
    usim= cat(1, usim, repmat(u,size(H,2),1));
    tsim= cat(1, tsim, t+tfill);
    xts= cat(1, xts, xt');
    tts= cat(1, tts, t+ts);
    
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
    'FontSize',FS, 'FontName','Times', 'FontWeight','bold')

hl= legend({'$x_2(t)$','$x_2(t)$'}, 'Interpreter', 'latex',...
    'FontSize',FS, 'FontName','Times', 'FontWeight','bold');

subplot(2,1,2)
plot(tsim, usim,'LineWidth',LW)
axis([0, Tsim, 0, 3])
xlabel({'$Time(s)$'}, 'interpreter','latex',...
    'FontSize',FS, 'FontName','Times', 'FontWeight','bold')
ylabel({'$u(t)$'}, 'interpreter','latex',...
    'FontSize',FS, 'FontName','Times', 'FontWeight','bold')

% plot 2-d state space trajectory
hf2= figure;
plot2_boxes(pavings(tag>0,:), [0.5,0.5,0.5], 'k', 1);
hold on
rectangle('Position',[G(1,1),G(2,1),G(1,2)-G(1,1),G(2,2)-G(2,1)],...
    'LineWidth',LW, 'LineStyle', '-')
p= plot(xts(:,1), xts(:,2),'-+','LineWidth',LW);
p.Color= [39,64,139]/255;
plot(x0(1), x0(2), 'o','MarkerFaceColor',cr, 'MarkerEdgeColor',cr)%[176 23 31]/255)[0.8500, 0.3250, 0.0980])
axis([X(1,1) X(1,2) X(2,1) X(2,2)])
xlabel({'$x_1$'}, 'interpreter','latex',...
    'FontSize',FS, 'FontName','Times', 'FontWeight','bold')
ylabel({'$x_2$'}, 'interpreter','latex',...
    'FontSize',FS, 'FontName','Times', 'FontWeight','bold')