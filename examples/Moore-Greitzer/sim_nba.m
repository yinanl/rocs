
addpath('../../matlab/')
%% load spec & controller
% data saved in .mat:
% - ctree, cindex, cvalue: Tree-structrued controller.
% - U : All input values.
% - X : Workarea.
% - ts: Sampling time.
% - G: Target area.
% - A: obstacles.
global a B H W lc cx cy aH H2 W2
a = 1/3.5;
B = 2.0;
H = 0.18;
W = 0.25;
lc = 8.0;
cx = 1.0/lc;
cy = 1.0/(4*lc*B*B);
aH = a+H;
H2 = H/(2.0*W*W*W);
W2 = 3*W*W;

spec= 'nba';
nNodes= 3;
controller.ctlr= cell(nNodes,1);
controller.pavings= cell(nNodes,1);
controller.tags= cell(nNodes,1);
for i=1:nNodes
    ctlrfile = strcat('data_', spec, '_w', num2str(i-1),'.mat');
    load(ctlrfile);
    controller.ctlr{i,1}= ctlr;
    controller.partition{i,1}= pavings;
    controller.tags{i,1}= tag;
    winid= find(any(ctlr,2));
    winset= pavings(winid,:);
    % winset= pavings(tag==1,:);
%     wc= [(winset(:,1)+winset(:,2))/2,...
%         (winset(:,3)+winset(:,4))/2];
    plot2_boxes(winset(:,1:4), [0.5,0.5,0.5], 'k', 1);
end

controller.X= X;
controller.U= U;
controller.ts= ts;

e= 0.003;
G= [0.4519-e, 0.4519+e;
    0.6513-e, 0.6513+e];
controller.G= {G};

xobs= [0.497, 0.503; 
    0.650, 0.656];
controller.A= xobs;

vf= @mg;


%% winning set
i= 1;
winset= controller.partition{i}(controller.tags{i}==1,:);
% winid= find(any(controller.ctlr{i},2)); % q=0
% winset= controller.partition{i}(winid,:);
wc= [(winset(:,1)+winset(:,2))/2,...
    (winset(:,3)+winset(:,4))/2];
% figure
% plot(wc(:,1), wc(:,2), '.')
% plot2_boxes(winset(:,1:4), [0.5,0.5,0.5], 'k', 1);


%% simulation
x0= [0.5343; 0.6553];
muold= 0.66;

T= 10;

tspan= [0, ts];
x= x0;  % a column
t= 0;
q= [0, 1, 2];
qid= 1;

tsim= t;
xsim= x';  % a row
usim= [0, 0];
qsim= q(qid);
rng('shuffle')
while(t<T)
    %%% compute control input %%%
    if (qid < numel(q))
        par= controller.partition{ q(qid)+1 };
        ctlr= controller.ctlr{ q(qid)+1 };
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
    
    %%%
    % Update the state of the DBA
    %%%
    switch q(qid)
        case 0 % at q=0
            if( x(1)>=G(1,1) && x(1)<=G(1,2) &&...
                    x(2)>=G(2,1) && x(2)<=G(2,2) ) % b is reached
                % jump to q=1
                qid= 2;
            else
                % stay in q=0
                qid= 1;
            end
        case 1 % at q=1
            if( x(1)>=G(1,1) && x(1)<=G(1,2) &&...
                    x(2)>=G(2,1) && x(2)<=G(2,2) ) % b is reached
                % stay in q=1
                qid= 2;
            else
                % jump to q=2
                qid= 3;
            end
        case 2 % trapped in q=2
            qid= 2;
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

% % simulated path
plot(xsim(:,1),xsim(:,2),...
    'Marker','.','MarkerEdgeColor','r','MarkerFaceColor','r')
plot(x0(1), x0(2), 'LineWidth', LW,...
    'Marker', '^', 'MarkerSize', 10,...
    'MarkerEdgeColor', colors(7,:), 'MarkerFaceColor', colors(7,:));
plot(xsim(end,1), xsim(end,2), 'LineWidth', LW,...
    'Marker', 'v', 'MarkerSize', 10,...
    'MarkerEdgeColor', colors(5,:), 'MarkerFaceColor', colors(5,:));
% for k= 1:4
%     xs= xsim(qsim==q(k),:);
%     plot(xs(:,1),xs(:,2), 'color', colors(k,:), 'LineWidth', LW)
% end


% axis equal
axis([X(1,:) X(2,:)])
xlabel({'$x$'}, 'interpreter','latex',...
    'Units','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')
ylabel({'$y$'}, 'interpreter','latex',...
    'Units','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')


% % % time-control curves
% usim= [usim(2:end, :); usim(end,:)];
% tq= [0:0.01:tsim(end)]';
% uq= interp1(tsim,usim,tq,'previous');
% 
% hf2=figure;
% subplot(2,1,1)
% plot(tq, uq, 'LineWidth', LW);
% axis([0, tsim(end), -2, 2])
% xlabel({'$t(s)$'}, 'interpreter','latex',...
%     'Units','points', 'FontSize',FS,...
%     'FontName','Times', 'FontWeight','bold')
% ylabel({'$u_1,\;u_2$'}, 'interpreter','latex',...
%     'Units','points', 'FontSize',FS,...
%     'FontName','Times', 'FontWeight','bold')
% 
% hl= legend({'$u_1(t)$','$u_2(t)$'}, 'Interpreter', 'latex',...
%     'Units','points', 'FontSize',FS,...
%     'FontName','Times', 'FontWeight','bold');
% 
% % % time-node switching
% qsim= [qsim(2:end, :); qsim(end,:)];
% kq= interp1(tsim,qsim,tq,'previous');
% 
% ax2= subplot(2,1,2);
% hold(ax2,'on');
% for k= 1:numel(q)-1
% %     kqk= [0;kq]==q(k);
%     kqsec= diff([0;kq]==q(k));
%     kqs= find(kqsec>0); % starting index
%     kqe= find(kqsec<0); % ending index
%     if(isempty(kqs))
%         kqs= 1;
%     end
%     if(isempty(kqe))
%         kqe= numel(kq);
%     end
% %     if(numel(kqe)<numel(kqs))
% %         kqe= [kqe; numel(kq)];
% %     end
%     for j= 1:numel(kqs)
%         plot(ax2, tq(kqs(j):kqe(j)), kq(kqs(j):kqe(j)), 'color', colors(k,:), 'LineWidth', LW)
%     end
% end
% plot(ax2, [tq(end);tq(end)], [3;4], 'color', colors(5,:), 'LineWidth', LW)
% % subplot(2,1,2)
% % plot(tq, 4-kq, 'LineWidth', 2);
% axis([0, tsim(end), 0, 4])
% xlabel({'$t(s)$'}, 'interpreter','latex',...
%     'Units','points', 'FontSize',FS,...
%     'FontName','Times', 'FontWeight','bold')
% ylabel({'$q_k$'}, 'interpreter','latex',...
%     'Units','points', 'FontSize',FS,...
%     'FontName','Times', 'FontWeight','bold')