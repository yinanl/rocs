
addpath('../../matlab')
%% load problem & controller
% data saved in .mat:
% - ctree, cindex, cvalue: Tree-structrued controller.
% - U : Two modes.
% - X : Workarea.
% - G: Target area.
% - ts: Sampling time.
load('data_U_grid.mat')
load('data_abstbased.mat')
X= [0 10; 0 10; -3.4 3.4];
G= [9 9.5; 0 0.5; -3.4 3.4];
ts= 0.3;
fm= @car;
vf= @car;

gid= goalset+ones(numel(goalset),1);
oid= avoidset+ones(numel(avoidset),1);
obs=xgrid(oid,:);
goal= xgrid(gid,:);
win= xgrid(any(leastctlr,2),:);

% figure
% hold on
% plot(obs(:,1), obs(:,2), '.','color','k')
% % plot3(win(:,1),win(:,2),win(:,3),'o','color','b')
% plot(win(:,1),win(:,2),'o','color','b')



%% simulation
% x0= [7.6; 0.4; pi/2];
x0= [0.6; 0.6; pi/2];
% x0= [9.8; 3.2; 0];

tspan= [0, ts];
x= x0;
t= 0;

tsim= t;
xsim= x';
usim= [0, 0];

xxx= repmat([x,x],1,1,size(G,3));
dGx= x-G;
flaginG= dGx(1,1,:)>0 & dGx(2,1,:)>0 & dGx(3,1,:)>0 ...
    & dGx(1,2,:)<0 & dGx(2,2,:)<0 & dGx(3,2,:)<0;
while(isempty(find(flaginG,1)))  % not reach goal
    
    % compute control input
    xid= get_stateid(xsim(end,:), xgrid);
    uid= find(leastctlr(xid(1),2:end));
    
    if (isempty(uid))
        error("Invalid controller.");
    else
%         uall= U(uid,:);
%         [val, ind]= min(abs(uall(:,1)));
%         u= uall(ind,:);
        u= U(uid(1),:);
    end
    
%     xx= fm(taus,x',u);
%     t= t+taus;
%     x= xx';
%     xsim= [xsim; xx];
%     tsim= [tsim; t];
%     usim= [usim; u];
    
    [tt, xx]= ode45(@(t,x) vf(t,x,u), tspan, x);
    x= xx(end,:)';
    t= t + tt(end,:);
    xsim= [xsim; xx(end,:)];
    tsim= [tsim; t];
    usim= [usim; u];
    
    xxx= repmat([x,x],1,1,size(G,3));
    dGx= x-G;
    flaginG= dGx(1,1,:)>0 & dGx(2,1,:)>0 & dGx(3,1,:)>0 ...
        & dGx(1,2,:)<0 & dGx(2,2,:)<0 & dGx(3,2,:)<0;
end

%% display
% define color
pink = [255,182,193]/255;
gold = [1,0.84,0];
lightgold = [255,236,139]/255;
gray = [0.6,0.6,0.6];
lightblue = [176 226 255]/255;
orange= [0.8500 0.3250 0.0980];

FS= 16; % fontsize
LW= 1.5; % lineweight

% car size
w=0.16;
h=0.08;

figure
hold on

% % whole area
rectangle('Position', [X(1,1), X(2,1), X(1,2)-X(1,1), X(2,2)-X(2,1)],...
    'EdgeColor','k', 'LineWidth',2)

% % goal area
rectangle('Position', [G(1,1), G(2,1), ...
    G(1,2)-G(1,1), G(2,2)-G(2,1)],...
    'EdgeColor',gold,'FaceColor',gold)
plot(goal(:,1),goal(:,2),'*','color','r')

% % avoid area
walls= [1,1.2, 0,9; 
    2.2,2.4, 0,5;
    2.2,2.4, 6,10;
    3.4,3.6, 0,9;
    4.6,4.8, 1,10;
    5.8,6, 0,6;
    5.8,6, 7,10;
    7,7.2, 1,10;
    8.2,8.4, 0,8.5;
    8.4,9.3, 8.3,8.5;
    9.3,10, 7.1,7.3;
    8.4,9.3, 5.9,6.1;
    9.3,10, 4.7,4.9;
    8.4,9.3, 3.5,3.7;
    9.3,10, 2.3,2.5];
for i=1:size(walls,1)
    rectangle('Position', [walls(i,1), walls(i,3), ...
        walls(i,2)-walls(i,1), walls(i,4)-walls(i,3)],...
        'EdgeColor',0.6*ones(3,1), 'FaceColor',0.6*ones(3,1))
end
plot(obs(:,1), obs(:,2), '.','color','k')
axis([X(1,:) X(2,:)])


% simulated path
for i= 1:size(xsim,1)
    plot_rectangle_angle(xsim(i,1),xsim(i,2),w,h,xsim(i,3))
end
plot(xsim(:,1),xsim(:,2),...
    'Marker','.','MarkerEdgeColor','r','MarkerFaceColor','r')
plot(xsim(:,1),xsim(:,2))

% axis equal
axis([X(1,:) X(2,:)])

xlabel('$x$', 'interpreter','latex',...
    'Units','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')
ylabel('$y$', 'interpreter','latex',...
    'Units','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')


% % % time-control curves
% usim= [usim(2:end, :); usim(end,:)];
% tq= [0:0.01:tsim(end)]';
% uq= interp1(tsim,usim,tq,'previous');
% 
% hf2=figure;
% plot(tq, uq, 'LineWidth', 2);
% axis([0, tsim(end), -2, 2])
% xlabel({'$t(s)$'}, 'interpreter','latex',...
%     'Units','points', 'FontSize',FS,...
%     'FontName','Times', 'FontWeight','bold')
% ylabel({'$u_1,\;u_2$'}, 'interpreter','latex',...
%     'FontUnits','points', 'FontSize',FS,...
%     'FontName','Times', 'FontWeight','bold')
% 
% hl= legend({'$u_1(t)$','$u_2(t)$'}, 'Interpreter', 'latex',...
%     'Units','points', 'FontSize',FS,...
%     'FontName','Times', 'FontWeight','bold');