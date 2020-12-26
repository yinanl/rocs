clear
clc
addpath('../../matlab')
%% load spec & controller
fm= @twoIntegrators;
vf= @scara;

G1= [0.4980, 0.5772; 1.5739, 1.7055; -0.1 0.1; -0.1, 0.1];
G2= [0.4903, 0.6069; -0.9363, -0.8363; -0.1 0.1; -0.1, 0.1];

%%% Read the specification %%%
% - n_dba: # of DBA nodes.
% - n_props: # of all propositions (2^AP)
% - q0: the initial node of the DBA.
% - acc: accepting nodes.
% - q_prime: the DBA transition matrix.
spec= 'gb2';
[n_dba,n_props,q0,acc,q_prime]=read_spec(strcat(spec,'.txt'));

%%% Controller data %%%
% - ts: Sampling time.
% - U : All input values.
% - X : Workarea.
% - A(can be empty): obstacles.
% - G(can be empty): Target area.
% - pavings: Tree-structrued controller.
% - tag: indicating if a cell is inside the winning set.
% - ctlr: all valid control inputs for each cell in pavings.
controller.ctlr= cell(n_dba,1);
controller.partitions= cell(n_dba,1);
controller.tags= cell(n_dba,1);
controller.obs= cell(n_dba,1);

%%% If controllers are stored in *.mat %%%
% for i=1:n_dba
%     ctlrfile = strcat('data_', spec, '_w', num2str(i-1),'.mat');
%     load(ctlrfile);
%     controller.ctlr{i}= ctlr;
%     controller.partitions{i}= pavings;
%     controller.tags{i}= tag;
%     controller.obs{i}= controller.partitions{i}(controller.tags{i}==-1,:);
%     figure
%     plot2_boxes([controller.obs{i}(:,1:2),controller.obs{i}(:,3:4)],...
%         [0.5,0.5,0.5], 'k', 1)
% end

%%% If controllers are stored in *.h5 %%%
for i=1:n_dba
    ctlrfile = strcat('controller_', spec, '_w', num2str(i-1),'.h5');
    controller.ctlr{i}= h5read(ctlrfile, '/ctlr')';
    controller.partitions{i}= h5read(ctlrfile, '/pavings')';
    controller.tags{i}= h5read(ctlrfile, '/tag');
    controller.obs{i}= controller.partitions{i}(controller.tags{i}==-1,:);
end
X= h5read(ctlrfile, '/X')';
U= h5read(ctlrfile, '/U')';
ts= h5read(ctlrfile, '/ts');


%% simulation
x0= [0;0;0;0];
T= 50;

x= x0;  % a column
t= 0;
% q= [0,1,2];
qid= q0+1;
nacc= 0;

tsim= [];
xsim= [];  % a row
usim= [];
torqsim= [];
qsim= [];
rng('shuffle')
while(t<T)
    %%% compute control input %%%
    if (qid < n_dba+1)
        par= controller.partitions{qid};
        ctlr= controller.ctlr{qid};
        xid= find(x(1)>=par(:,1) & x(1)<=par(:,2) & ...
            x(2)>=par(:,3) & x(2)<=par(:,4) & ...
            x(3)>=par(:,5) & x(3)<=par(:,6) & ...
            x(4)>=par(:,7) & x(4)<=par(:,8));
        uid= find(ctlr(xid(1),:));
        if (isempty(uid))
            s1= 'Invalid controller for q';
            s2= num2str(qid-1);
            s= strcat(s1,s2);
            error(s);
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
    end
    
    %%% Convert from the angular accelarations to torques %%%
    torqsim_inter= compute_torque(u', x(3:4), x(1:2));
    
    tsim= [tsim; t];
    xsim= [xsim; x'];
    usim= [usim; u];
    torqsim= [torqsim; torqsim_inter'];
    qsim= [qsim; qid-1];
    
    %%% Update the state of the DBA %%%
    p= labeling(x, G1, G2);
    qid= q_prime(qid, p+1);
    if(qid-1==acc)
        nacc= nacc+1;
    end
    
    %%% Update the state of the dyanmical system %%%
    %%% use discrete-time two double integrators
    xx= fm(ts,x',u);
    x= xx';
    t= t+ts;
    %%% use the full scara ode
%     [tt, xx]= ode45(@(t,x) vf(t,x,torq), [0, ts], x);
%     x= xx(end,:)';
%     t= t + tt(end,:);
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

% % display the obstacles and target
figure
hold on
plot2_boxes([controller.obs{q0+1}(:,1:2),controller.obs{q0+1}(:,3:4)], ...
    [0.5,0.5,0.5], 'k', 1)
% % reference obstacles
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
% % display target area
rectangle('Position', [G1(1,1), G1(2,1), ...
    G1(1,2)-G1(1,1), G1(2,2)-G1(2,1)],...
    'EdgeColor','g','FaceColor','g')
rectangle('Position', [G2(1,1), G2(2,1), ...
    G2(1,2)-G2(1,1), G2(2,2)-G2(2,1)],...
    'EdgeColor','g','FaceColor','g')

% % display phase trajectories
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

% % display computed torque
tq= [0:0.01:tsim(end)]';
torqsim_inter= interp1(tsim,torqsim,tq,'previous');

figure
plot(tq, torqsim_inter, 'LineWidth', LW)
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