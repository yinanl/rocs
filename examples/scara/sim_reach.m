clear
clc
addpath('../../matlab')
%% load spec & controller
vf= @scara;
fm= @twoIntegrators;

%%% Controller data %%%
% - U : All input values.
% - X : Workarea.
% - ts: Sampling time.
% - G(can be empty): Target area.
% - A(can be empty): obstacles.
% - pavings: Tree-structrued controller.
% - tag: indicating if a cell is inside the winning set.
% - ctlr: all valid control inputs for each cell in pavings.

%%% Load from *.mat file %%%
% load('data_2dbint_reach1.mat');

%%% Load from .h5 file %%%
ctlrfile= 'controller_2dbint_reach1.h5';
ts= h5read(ctlrfile, '/ts');
X= h5read(ctlrfile, '/X')';
U= h5read(ctlrfile, '/U')';
G= permute(h5read(ctlrfile, '/G'), [3,2,1]);
pavings= h5read(ctlrfile, '/pavings')';
tag= h5read(ctlrfile, '/tag');
ctlr= h5read(ctlrfile, '/ctlr')';


%% simulation
% x0= [0; 0.0; 0; 0];
% x0= [0; 1.2; 0; 0];
x0= [0.5; -0.88; 0.05; -0.05];

tspan= [0, ts];
x= x0;
t= 0;

tsim= [];
xsim= [];
usim= [];
xidsim= [];
torqsim= [];
while(x(1)>G(1,2) || x(1)<G(1,1) ||...
        x(2)>G(2,2) || x(2)<G(2,1) || ...
        x(3)>G(3,2) || x(3)<G(3,1) || ...
        x(4)>G(4,2) || x(4)<G(4,1))  % not reach goal
    
    % compute control input
    xid= find(x(1)>=pavings(:,1) & x(1)<=pavings(:,2) & ...
        x(2)>=pavings(:,3) & x(2)<=pavings(:,4) & ...
        x(3)>=pavings(:,5) & x(3)<=pavings(:,6) & ...
        x(4)>=pavings(:,7) & x(4)<=pavings(:,8));
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
    
    %%% Convert from the angular accelarations to torques %%%
    torque= compute_torque(u', x(3:4), x(1:2));
    
    xsim= [xsim; x'];
    tsim= [tsim; t];
    usim= [usim; u];
    xidsim= [xidsim; xid];
    torqsim= [torqsim; torque'];
    
    y= fm(ts, x, u);
    t= t + ts;
    x= y';
    
%     [tt, xx]= ode45(@(t,x) vf(t,x,torque), [0, ts], x);
%     t= t+tt(end);
%     x= xx(end,:)';
end


%% display
%%% define color %%%
pink = [255,182,193]/255;
gold = [1,0.84,0];
lightgold = [255,236,139]/255;
gray = [0.6,0.6,0.6];
lightblue = [176 226 255]/255;
orange= [0.8500 0.3250 0.0980];

FS= 16; % fontsize
LW= 1.5; % lineweight

%%% display the obstacles and target %%%
figure
hold on
obs= pavings(tag==-1,:);
oc= [(obs(:,1)+obs(:,2))/2,(obs(:,3)+obs(:,4))/2,...
    (obs(:,5)+obs(:,6))/2,(obs(:,7)+obs(:,8))/2];
plot2_boxes([obs(:,1:2),obs(:,3:4)], [0.5,0.5,0.5], 'k', 1)
%%% reference obstacles %%%
l1 = 0.15;
l2 = 0.15;
h = 0.8*l1;
r = 0.5*l1;
a1 = atan2(h,r);
a2= asin(h/l1);
N1= 50;
N2= 30;
theta11= 0:a2/N1:a2;
theta21= pi - theta11 - atan2(h-l1.*sin(theta11), l1.*cos(theta11)-r);
theta12= a2:(a1-a2)/N2:a1;
theta22= pi - theta12 + atan2(l1.*sin(theta12)-h, l1.*cos(theta12));
theta1= [theta11, theta12];
theta2= [theta21, theta22];
plot(theta1, theta2, 'LineWidth', 2)

%%% Display target set %%%
rectangle('Position', [G(1,1), G(2,1), ...
    G(1,2)-G(1,1), G(2,2)-G(2,1)],...
    'EdgeColor','g','FaceColor','g')

%%% Display winning set %%%
win= pavings(tag==1,:);
wc= [(win(:,1)+win(:,2))/2,(win(:,3)+win(:,4))/2,...
    (win(:,5)+win(:,6))/2,(win(:,7)+win(:,8))/2];
% plot(wc(:,1), wc(:,2),'.')
% figure(1)
% plot2_boxes(win(:,1:4), [0.5,0.5,0.5], 'k', 1) % theta1-theta2
% figure(2)
% plot2_boxes([win(:,1:2),win(:,5:6)], [0.5,0.5,0.5], 'k', 1) % theta1-w1
% figure(3)
% plot2_boxes([win(:,3:4),win(:,7:8)], [0.5,0.5,0.5], 'k', 1) % theta2-w2

%%% display phase trajectories %%%
% plot(xsim(:,1),xsim(:,2),...
%     'Marker','.','MarkerEdgeColor','r','MarkerFaceColor','r')
plot(xsim(:,1),xsim(:,2), 'LineWidth', LW)

axis([X(1,:) X(2,:)])
xlabel({'$\theta_1$ (rad)'}, 'interpreter','latex',...
    'Units','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')
ylabel({'$\theta_2$ (rad)'}, 'interpreter','latex',...
    'Units','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')

%%% Display computed torque %%%
wsim= xsim(1:end-1,3:4)';
thesim= xsim(1:end-1,1:2)';
tq= [0:0.01:tsim(end)]';
torq= interp1(tsim(1:end),torqsim,tq,'previous');

figure
plot(tq, torq, 'LineWidth', LW)
axis([0, tsim(end), min(min(torqsim)), max(max(torqsim))])
xlabel({'$t(s)$'}, 'interpreter','latex',...
    'Units','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')
ylabel({'$\tau_1,\;\tau_2$'}, 'interpreter','latex',...
    'Units','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')

hl= legend({'$\tau_1(t)$','$\tau_2(t)$'}, 'Interpreter', 'latex',...
    'Units','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold');