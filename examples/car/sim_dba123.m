clear
clc
addpath('../../matlab/')
%% Work space setup
G1= [1,2; 0.5,2;-pi,pi];
G2= [0.5,2.5;7.5,8.5;-pi,pi];
G3= [7.1,9.1;4.6,6.4;-pi,pi];
G= {G1, G2, G3};

A= zeros(3,2,4);
A(:,:,1)= [1.6,5.7;4.0,5.0;-pi,pi];
A(:,:,2)= [3,5;5,8;-pi,pi];
A(:,:,3)= [4.3,5.7;1.8,4.0;-pi,pi];
A(:,:,4)= [5.7,8.5;1.8,2.5;-pi,pi];


%% load spec & controller
%%% Define ODEs or DEs for car kinematics %%%
vf= @car; % ODEs
fm= @carflow; % DEs

%%% Read the specification %%%
% - n_dba: # of DBA nodes.
% - n_props: # of all propositions (2^AP)
% - q0: the initial node of the DBA.
% - acc: accepting nodes.
% - q_prime: the DBA transition matrix.
spec= 'dba3';
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
controller.ctlrs= cell(n_dba,1);
controller.partitions= cell(n_dba,1);
controller.tags= cell(n_dba,1);

%%% If controllers are stored in *.mat %%%
% for i=1:n_dba
%     ctlrfile = strcat('data_', spec, '_w', num2str(i-1),'.mat');
%     load(ctlrfile);
%     controller.ctlrs{i,1}= ctlr;
%     controller.partitions{i,1}= pavings;
%     controller.tags{i,1}= tag;
% end

%%% If controllers are stored in *.h5 %%%
for i=1:n_dba
    ctlrfile = strcat('controller_', spec, '_w', num2str(i-1),'.h5');
    controller.ctlrs{i}= h5read(ctlrfile, '/ctlr')';
    controller.partitions{i}= h5read(ctlrfile, '/pavings')';
    controller.tags{i}= h5read(ctlrfile, '/tag');
end
ts= h5read(ctlrfile, '/ts');
X= h5read(ctlrfile, '/X')';
U= h5read(ctlrfile, '/U')';


%% simulation
rng('shuffle')
T= 50;
num_acc= 3;

%%% Set up an initial condition %%%
% x0= [1.3; 5; 3*pi/4];
% x0= [3; 2; pi/2];
x0= [6; 1; pi/2];
% x0= [9; 5; pi/4];

x= x0;  % a column
t= 0;
qid= q0+1;
nacc= 0;

tsim= [];
xsim= [];
usim= [];
qsim= [];
while(t<T || nacc<num_acc)
    %%% compute control input %%%
    if (qid < n_dba+1)
        par= controller.partitions{qid};
        ctlr= controller.ctlrs{qid};
        xid= find(x(1)>=par(:,1) & x(1)<=par(:,2) & ...
            x(2)>=par(:,3) & x(2)<=par(:,4) & ...
            x(3)>=par(:,5) & x(3)<=par(:,6));
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
    
    %%% Store simulation data %%%
    if(x(3)>pi)
        x(3) = x(3)-2*pi;
    end
    if(x(3)<-pi)
        x(3) = x(3)+2*pi;
    end
    xsim= [xsim; x'];
    tsim= [tsim; t];
    usim= [usim; u];
    qsim= [qsim; qid-1];
    
    %%% Update the state of the DBA %%%
    p= labeling1(x, G1, G2, G3);
    qid= q_prime(qid, p+1);
    if(qid-1==acc)
        nacc= nacc+1;
    end
    
    %%% Update the state of the dyanmical system %%%
    %%% use difference equations
    %     xx= fm(taus,x',u);
    %     t= t+taus;
    %     x= xx';
    %     xsim= [xsim; xx];
    %     tsim= [tsim; t];
    %     usim= [usim; u];
    %%% use ode
    [tt, xx]= ode45(@(t,x) vf(t,x,u), [0, ts], x);
    
    %%% update x, t %%%
    x= xx(end,:)';
    t= t + tt(end,:);
end


%% display
%%% define color %%%
pink = [255,182,193]/255;
gold = [1,0.84,0];
lightgold = [255,236,139]/255;
gray = [0.6,0.6,0.6];
lightblue = [176 226 255]/255;
orange= [0.8500 0.3250 0.0980];
colors=get(groot,'DefaultAxesColorOrder');

FS= 16; % fontsize
LW= 2; % lineweight

%%% car size %%%
w=0.16;
h=0.08;

%%% Display workspace %%%
figure
hold on
%%% whole area
rectangle('Position', [X(1,1), X(2,1), ...
    X(1,2)-X(1,1), X(2,2)-X(2,1)],...
    'EdgeColor','k', 'LineWidth',2)
%%% goal area
for k=1:size(G,2)
    goal= G{k};
    for i=1:size(goal,3)
        rectangle('Position', [goal(1,1,i), goal(2,1,i), ...
            goal(1,2,i)-goal(1,1,i), goal(2,2,i)-goal(2,1,i)],...
            'EdgeColor',gold,'FaceColor',gold)
    end
end
%%% avoid area
if (~isempty(A))
    for i= 1:size(A, 3)
        rectangle('Position', [A(1,1,i), A(2,1,i), ...
            A(1,2,i)-A(1,1,i), A(2,2,i)-A(2,1,i)],...
            'EdgeColor',gray, 'FaceColor',gray)
    end
end
%%% Labels for areas
text(1.4,1.2,'$a_1$', 'interpreter','latex', 'FontName','Times', 'FontSize',FS)
text(1.3,8.0,'$a_2$', 'interpreter','latex', 'FontName','Times', 'FontSize',FS)
text(8.0,5.5,'$a_3$', 'interpreter','latex', 'FontName','Times', 'FontSize',FS)
text(4,5,'$o$', 'interpreter','latex', 'FontName','Times', 'FontSize',FS)


%%% Simulated path %%%
for i= 1:size(xsim,1)
    plot_rectangle_angle(xsim(i,1),xsim(i,2),w,h,xsim(i,3))
end
% plot(xsim(:,1),xsim(:,2),...
%     'Marker','.','MarkerEdgeColor','r','MarkerFaceColor','r')
plot(x0(1), x0(2), 'LineWidth', LW,...
    'Marker', '^', 'MarkerSize', 10,...
    'MarkerEdgeColor', colors(7,:), 'MarkerFaceColor', colors(7,:));
plot(xsim(end,1), xsim(end,2), 'LineWidth', LW,...
    'Marker', 'v', 'MarkerSize', 10,...
    'MarkerEdgeColor', colors(5,:), 'MarkerFaceColor', colors(5,:));

%%% plot trajectory with different colors for different DBA state %%%
qdiff = diff([q0; qsim]);
qchks= [1; find(qdiff); numel(qsim)];
q = q0;
for k= 1:numel(qchks)-1
    q = q + qdiff(qchks(k));
    xs= xsim(qchks(k):qchks(k+1), :);
    plot(xs(:,1),xs(:,2), 'color', colors(q+1,:), 'LineWidth', LW)
end

% axis equal
axis([X(1,:) X(2,:)])
xlabel({'$x$'}, 'interpreter','latex',...
    'Units','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')
ylabel({'$y$'}, 'interpreter','latex',...
    'Units','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')


%%% Refine time scale %%%
tq= [0:0.01:tsim(end)]';
uq= interp1(tsim,usim,tq,'previous');
kq= interp1(tsim,qsim,tq,'previous');


hf2=figure;
%%% time-control curves %%%
subplot(2,1,1)
plot(tq, uq, 'LineWidth', LW);
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

%%% time-q curves %%%
ax2= subplot(2,1,2);
hold(ax2,'on');
kdiff = diff([q0; kq]);
kchks= [1; find(kdiff); numel(kq)];
q = q0;
for k= 1:numel(kchks)-1
    q = q + kdiff(kchks(k));
    plot(ax2, tq(kchks(k):kchks(k+1)), kq(kchks(k):kchks(k+1)),...
        'color', colors(q+1,:), 'LineWidth', LW)
end
axis([0, tsim(end), -1, n_dba])
xlabel({'$t(s)$'}, 'interpreter','latex',...
    'Units','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')
ylabel({'$q_k$'}, 'interpreter','latex',...
    'Units','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')