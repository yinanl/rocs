
addpath('../../matlab/')
%% load spec & controller
% data saved in .mat:
% - ctree, cindex, cvalue: Tree-structrued controller.
% - U : All input values.
% - X : Workarea.
% - ts: Sampling time.
% - G: Target area.
% - A: obstacles.

G1= [1,2; 0.5,2;-pi,pi];
G2= [0.5,2.5;7.5,8.5;-pi,pi];
G3= [7.1,9.1;4.6,6.4;-pi,pi];
xobs= zeros(3,2,4);
xobs(:,:,1)= [1.6,5.7;4.0,5.0;-pi,pi];
xobs(:,:,2)= [3,5;5,8;-pi,pi];
xobs(:,:,3)= [4.3,5.7;1.8,4.0;-pi,pi];
xobs(:,:,4)= [5.7,8.5;1.8,2.5;-pi,pi];

% spec= 'p2';
spec= 'dba2';
nNodes= 4;
controller.ctlr= cell(nNodes,1);
controller.pavings= cell(nNodes,1);
controller.tags= cell(nNodes,1);
for i=1:nNodes
    ctlrfile = strcat('data_', spec, '_w', num2str(i-1),'.mat');
    load(ctlrfile);
    controller.ctlr{i,1}= ctlr;
    controller.partition{i,1}= pavings;
    controller.tags{i,1}= tag;
end

controller.X= X;
controller.U= U;
controller.ts= ts;
controller.G= {G1, G2, G3};
controller.A= xobs;

vf= @car;
fm= @carflow;


%% simulation
% x0= [1.3; 5; 3*pi/4];
% x0= [3; 2; pi/2];
% x0= [6; 1; pi/2];
x0= [9; 5; pi/4];

T= 60;

tspan= [0, ts];
x= x0;  % a column
t= 0;
q= [0, 1, 2, 3];
qid= 1;

tsim= t;
xsim= x';  % a row
usim= [0, 0];
qsim= q(qid);
while(t<T)
    %%% compute control input %%%
    par= controller.partition{ q(qid)+1 };
    ctlr= controller.ctlr{ q(qid)+1 };
    xid= find(x(1)>=par(:,1) & x(1)<=par(:,2) & ...
        x(2)>=par(:,3) & x(2)<=par(:,4) & ...
        x(3)>=par(:,5) & x(3)<=par(:,6));
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
    
    %%%
    % Update the state of the DBA
    %%%
    switch q(qid)
        case 3
            if( x(1)>=G1(1,1) && x(1)<=G1(1,2) &&...
                    x(2)>=G1(2,1) && x(2)<=G1(2,2) &&...
                    x(3)>=G1(3,1) && x(3)<=G1(3,2) ) % a1 is reached
                % jump to q0
                qid= 1; % q(1)= q0;
            else
                % stay in q3
                qid= 4;
            end
        case 2
            if( x(1)>=G2(1,1) && x(1)<=G2(1,2) &&...
                    x(2)>=G2(2,1) && x(2)<=G2(2,2) &&...
                    x(3)>=G2(3,1) && x(3)<=G2(3,2) ) % a2 is reached
                % jump to q3
                qid= 4; % q(4)= q3;
            else
                % stay in q2
                qid= 3;
            end
        case 1
            if( x(1)>=G3(1,1) && x(1)<=G3(1,2) &&...
                    x(2)>=G3(2,1) && x(2)<=G3(2,2) &&...
                    x(3)>=G3(3,1) && x(3)<=G3(3,2) ) % a3 is reached
                % jump to q2
                qid= 3;
            else
                % stay in q1
                qid= 2;
            end
        case 0
            if( x(1)>=G3(1,1) && x(1)<=G3(1,2) &&...
                    x(2)>=G3(2,1) && x(2)<=G3(2,2) &&...
                    x(3)>=G3(3,1) && x(3)<=G3(3,2) ) % a3 is reached
                % jump to q2
                qid= 3;
            else
                % jump in q1
                qid= 2;
            end
    end
    
    %%%
    % Update the state of the dyanmical system
    %%%
    %%% use difference equations %%%
    %     xx= fm(taus,x',u);
    %     t= t+taus;
    %     x= xx';
    %     xsim= [xsim; xx];
    %     tsim= [tsim; t];
    %     usim= [usim; u];
    %%% use ode %%%
    [tt, xx]= ode45(@(t,x) vf(t,x,u), tspan, x);
    x= xx(end,:)';
    t= t + tt(end,:);
    xsim= [xsim; xx(end,:)];
    tsim= [tsim; t];
    usim= [usim; u];
    qsim= [qsim; q(qid)];
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

% car size
w=0.16;
h=0.08;


figure
hold on

% % whole area
rectangle('Position', [controller.X(1,1), controller.X(2,1), ...
    controller.X(1,2)-controller.X(1,1), controller.X(2,2)-controller.X(2,1)],...
    'EdgeColor','k', 'LineWidth',2)

% % goal area
for k=1:size(controller.G,2)
    goal= controller.G{k};
    for i=1:size(goal,3)
        rectangle('Position', [goal(1,1,i), goal(2,1,i), ...
            goal(1,2,i)-goal(1,1,i), goal(2,2,i)-goal(2,1,i)],...
            'EdgeColor',gold,'FaceColor',gold)
    end
end

% % avoid area
if (~isempty(controller.A))
    for i= 1:size(controller.A, 3)
        rectangle('Position', [controller.A(1,1,i), controller.A(2,1,i), ...
            controller.A(1,2,i)-controller.A(1,1,i), controller.A(2,2,i)-controller.A(2,1,i)],...
            'EdgeColor',gray, 'FaceColor',gray)
    end
end

text(1.4,1.2,'$a_1$', 'interpreter','latex', 'FontName','Times', 'FontSize',FS)
text(1.3,8.0,'$a_2$', 'interpreter','latex', 'FontName','Times', 'FontSize',FS)
text(8.0,5.5,'$a_3$', 'interpreter','latex', 'FontName','Times', 'FontSize',FS)
text(4,5,'$o$', 'interpreter','latex', 'FontName','Times', 'FontSize',FS)

% % simulated path
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

for k= 1:numel(q)
%     qk= qsim==q(k);
%     qkchunk= [0;diff(qsim==q(k))];
    qkchunk= diff([0;qsim]==q(k));
    qkstart= find(qkchunk>0);
    kqe= find(qkchunk<0);
    if(numel(kqe)<numel(qkstart))
        kqe= [kqe; numel(qsim)];
    elseif(numel(kqe)>numel(qkstart))
        qkstart= [1; qkstart];
    end
    qkid= [qkstart kqe];
    for j= 1:size(qkid,1)
        xs= xsim(qkid(j,1):qkid(j,2),:);
        plot(xs(:,1),xs(:,2), 'color', colors(k,:), 'LineWidth', LW)
    end
%     xs= xsim(qsim==q(k) & qset<=qkchunk(j),:);
%     plot(xs(:,1),xs(:,2), 'color', colors(q(k),:), 'LineWidth', LW)
end
% axis equal
axis([X(1,:) X(2,:)])
xlabel({'$x$'}, 'interpreter','latex',...
    'Units','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')
ylabel({'$y$'}, 'interpreter','latex',...
    'Units','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')


% % time-control curves
usim= [usim(2:end, :); usim(end,:)];
tq= [0:0.01:tsim(end)]';
uq= interp1(tsim,usim,tq,'previous');

hf2=figure;
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

% % time-buchiset switching
qsim= [qsim(2:end, :); qsim(end,:)];
kq= interp1(tsim,qsim,tq,'previous');

ax2= subplot(2,1,2);
hold(ax2,'on');
for k= 1:numel(q)
%     kqk= [0;kq]==q(k);
    kqsec= diff([0;kq]==q(k));
    kqs= find(kqsec>0);
    kqe= find(kqsec<0);
    if(isempty(kqs))
        kqs= 1;
    end
    if(isempty(kqe))
        kqe= numel(kq);
    end
    if(numel(kqe)<numel(kqs))
        kqe= [kqe; numel(kq)];
    elseif(numel(kqe)>numel(kqs))
        kqs= [1; kqs];
    end
    for j= 1:numel(kqs)
        plot(ax2, tq(kqs(j):kqe(j)), kq(kqs(j):kqe(j)), 'color', colors(k,:), 'LineWidth', LW)
    end
end
% plot(tq, kq-1, 'LineWidth', 2);
axis([0, tsim(end), 0, 4])
xlabel({'$t(s)$'}, 'interpreter','latex',...
    'Units','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')
ylabel({'$q_k$'}, 'interpreter','latex',...
    'Units','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')