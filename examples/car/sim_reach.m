clear
clc
addpath('../../matlab/')
%% load spec & controller
%%% Define ODEs or DEs for car kinematics %%%
vf= @car; % ODEs
fm= @carflow; % DEs

%%% Controller data %%%
% - ts: Sampling time.
% - U : All input values.
% - X : Workarea.
% - A(can be empty): Obstacles.
% - G(can be empty): Target area.
% - pavings: Tree-structrued controller.
% - tag: indicating if a cell is inside the winning set.
% - ctlr: all valid control inputs for each cell in pavings.

%%% Load from .mat file %%%
% load('data_carReach1.mat')
% load('data_carReach2.mat')
% load('data_carBuchi.mat')

%%% Load from .h5 file %%%
ctlrfile= 'controller_carBuchi.h5';
ts= h5read(ctlrfile, '/ts');
X= h5read(ctlrfile, '/X')';
U= h5read(ctlrfile, '/U')';
A= permute(h5read(ctlrfile, '/xobs'), [3,2,1]);
G= permute(h5read(ctlrfile, '/G'), [3,2,1]);
pavings= h5read(ctlrfile, '/pavings')';
tag= h5read(ctlrfile, '/tag');
ctlr= h5read(ctlrfile, '/ctlr')';


%% simulation
x0= [7.6; 0.4; pi/2];
Tsim= 50;

x= x0;
t= 0;
tsim= [];
xsim= [];
usim= [];
while(t<Tsim && ...
        x(1)>G(1,2) || x(1)<G(1,1) ||...
        x(2)>G(2,2) || x(2)<G(2,1) || ...
        x(3)>G(3,2) || x(3)<G(3,1))  % not reach goal
    %%% compute control input %%%
    xid= find(x(1)>=pavings(:,1) & x(1)<=pavings(:,2) & ...
        x(2)>=pavings(:,3) & x(2)<=pavings(:,4) & ...
        x(3)>=pavings(:,5) & x(3)<=pavings(:,6));
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
    
    %%% Store simulation data %%%
    xsim= [xsim; x'];
    tsim= [tsim; t];
    usim= [usim; u];
    
    %%% Update the state of the dyanmical system %%%
    %%% use difference equations
%     xx= fm(taus,x',u);
%     t= t+taus;
%     x= xx';
    %%% use ode
    [tt, xx]= ode45(@(t,x) vf(t,x,u), [0, ts], x);
    
    %%% update x, t %%%
    x= xx(end,:)';
    t= t + tt(end,:);
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

% % avoid area
if (~isempty(A))
    for i= 1:size(A, 3)
        rectangle('Position', [A(1,1,i), A(2,1,i), ...
            A(1,2,i)-A(1,1,i), A(2,2,i)-A(2,1,i)],...
            'EdgeColor',gray, 'FaceColor',gray)
    end
end


% % simulated path
for i= 1:size(xsim,1)
    plot_rectangle_angle(xsim(i,1),xsim(i,2),w,h,xsim(i,3))
end
plot_rectangle_angle(x(1),x(2),w,h,x(3))
plot([xsim(:,1); x(1)],[xsim(:,2); x(2)],...
    'Marker','.','MarkerEdgeColor','r','MarkerFaceColor','r')
plot([xsim(:,1); x(1)],[xsim(:,2); x(2)], 'LineWidth', LW)

axis equal
axis([X(1,:) X(2,:)])

xlabel('x position')
ylabel('y position')


% % time-control curves
tq= [0:0.01:tsim(end)]';
uq= interp1(tsim,usim,tq,'previous');

hf2=figure;
plot(tq, uq, 'LineWidth', 2);
axis([0, tsim(end), -2, 2])
xlabel({'$t(s)$'}, 'interpreter','latex',...
    'Units','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')
ylabel({'$u_1,\;u_2$'}, 'interpreter','latex',...
    'Units','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')

hl= legend({'$u_1(t)$','$u_2(t)$'}, 'Interpreter', 'latex',...
    'Units','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold');

% legend('wheel velocity', 'steering angle')