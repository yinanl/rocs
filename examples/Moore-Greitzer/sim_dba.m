clear
clc
addpath('../../matlab/')
%% Set up work space
e= 0.003;
% G= [0.4519-e, 0.4519+e;
%     0.6513-e, 0.6513+e];
G= [0.5039-e, 0.5039+e; 0.6605-e, 0.6605+e];

% xobs= [0.497, 0.503; 
%     0.650, 0.656];
A= [0.520, 0.526; 0.658, 0.664];

vf= @MG;


%% load spec & controller

%%% Read the specification %%%
% - n_dba: # of DBA nodes.
% - n_props: # of all propositions (2^AP)
% - q0: the initial node of the DBA.
% - acc: accepting nodes.
% - q_prime: the DBA transition matrix.
spec= 'FGb';
[n_dba,n_props,q0,acc,q_prime]=read_spec(strcat(spec,'.txt'));

%%% Load controllers %%%
% data saved in .mat:
% - U : All input values.
% - X : Workarea.
% - ts: Sampling time.
% - G(can be empty): Target area.
% - A(can be empty): obstacles.
% - pavings: Tree-structrued controller.
% - tag: indicating if a cell is inside the winning set.
% - ctlr: all valid control inputs for each cell in pavings.
controller.ctlr= cell(n_dba,1);
controller.partitions= cell(n_dba,1);
controller.tags= cell(n_dba,1);
% % Using .mat
for i=1:n_dba
    ctlrfile = strcat('data_', spec, '_w', num2str(i-1),'.mat');
    load(ctlrfile);
    controller.ctlr{i,1}= ctlr;
    controller.partitions{i,1}= pavings;
    controller.tags{i,1}= tag;
end

% % Using .h5
% for i=1:n_dba
%     ctlrfile = strcat('controller_', spec, '_w', num2str(i-1),'.h5');
%     controller.ctlr{i}= h5read(ctlrfile, '/ctlr')';
%     controller.partitions{i}= h5read(ctlrfile, '/pavings')';
%     controller.tags{i}= h5read(ctlrfile, '/tag');
% end
% ts= h5read(ctlrfile, '/ts');
% X= h5read(ctlrfile, '/X')';
% U= h5read(ctlrfile, '/U')';


%% winning set
winset= controller.partitions{q0+1}(controller.tags{q0+1}==1,:);
% winid= find(any(controller.ctlr{i},2)); % q=0
% winset= controller.partition{i}(winid,:);
wc= [(winset(:,1)+winset(:,2))/2,...
    (winset(:,3)+winset(:,4))/2];
% figure
% plot(wc(:,1), wc(:,2), '.')
% plot2_boxes(winset(:,1:4), [0.5,0.5,0.5], 'k', 1);
% axis([X(1,:) X(2,:)])


%% simulation
x0= [0.5343; 0.6553];
Tsim= 10;
num_acc= 10;

x= x0;  % a column
t= 0;
qid= q0+1;
nacc= 0;

rng('shuffle')
tsim= [];
xsim= [];  % a row
usim= [];
qsim= [];
while(t<Tsim || nacc<num_acc)
    %%% compute control input %%%
    if (qid <= n_dba)
        par= controller.partitions{qid};
        ctlr= controller.ctlr{qid};
        xid= find(x(1)>=par(:,1) & x(1)<=par(:,2) & ...
            x(2)>=par(:,3) & x(2)<=par(:,4));
        uid= find(ctlr(xid(1),:));
        if (isempty(uid))
            s1= 'Invalid controller for q';
            s2= num2str(q(qid));
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
    
    xsim= [xsim; x'];
    tsim= [tsim; t];
    usim= [usim; u];
    qsim= [qsim; qid-1];
    
    %%% Update the state of the DBA %%%
    p= labeling(x, G, A);
    qid= q_prime(qid, p+1);
    if(qid-1==acc)
        nacc= nacc+1;
    end
    
    %%% Update the state of the dyanmical system %%%
    [tt, xx]= ode45(@(t,x) vf(t,x,u), [0, ts], x);
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

colors=get(groot,'DefaultAxesColorOrder');

FS= 16; % fontsize
LW= 2; % lineweight


figure
hold on

% % whole area
rectangle('Position', [X(1,1), X(2,1), X(1,2)-X(1,1), X(2,2)-X(2,1)],...
    'EdgeColor','k', 'LineWidth', LW)

% % goal area
rectangle('Position', [G(1,1), G(2,1), G(1,2)-G(1,1), G(2,2)-G(2,1)],...
        'EdgeColor',gold,'FaceColor',gold)

% % avoid area
rectangle('Position', [A(1,1), A(2,1), A(1,2)-A(1,1), A(2,2)-A(2,1)],...
    'EdgeColor',gray, 'FaceColor',gray)

% % simulated path
plot(xsim(:,1), xsim(:,2), 'LineWidth', 1.5)
plot(x0(1), x0(2), 'LineWidth', LW,...
    'Marker', '^', 'MarkerSize', 5,...
    'MarkerEdgeColor', colors(7,:), 'MarkerFaceColor', colors(7,:));
plot(xsim(end,1), xsim(end,2), 'LineWidth', LW,...
    'Marker', 'v', 'MarkerSize', 5,...
    'MarkerEdgeColor', colors(5,:), 'MarkerFaceColor', colors(5,:));

axis([X(1,:) X(2,:)])
% axis equal
xlabel({'$x$'}, 'interpreter','latex',...
    'Units','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')
ylabel({'$y$'}, 'interpreter','latex',...
    'Units','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')


%%% Interpolate data %%%
tq= [0:0.01:tsim(end)]';
uq= interp1(tsim,usim,tq,'previous');
kq= interp1(tsim,qsim,tq,'previous');
hf2=figure;

%%% time-control curves %%%
subplot(2,1,1)
plot(tq, uq, 'LineWidth', LW);
axis([0, tsim(end), min(min(usim)), max(max(usim))])
xlabel({'$t(s)$'}, 'interpreter','latex',...
    'Units','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')
ylabel({'$u_1,\;u_2$'}, 'interpreter','latex',...
    'Units','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')
hl= legend({'$u_1(t)$','$u_2(t)$'}, 'Interpreter', 'latex',...
    'Units','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold');

%%% time-node switching %%%
ax2= subplot(2,1,2);
hold(ax2,'on');
plot(ax2, tq, kq, 'color', colors(5,:), 'LineWidth', LW)
axis([0, tsim(end), -1, n_dba])
xlabel({'$t(s)$'}, 'interpreter','latex',...
    'Units','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')
ylabel({'$q_k$'}, 'interpreter','latex',...
    'Units','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')