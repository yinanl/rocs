
addpath('../../matlab/')
%% load spec & controller
% data saved in .mat:
% - ctree, cindex, cvalue: Tree-structrued controller.
% - U : All input values.
% - X : Workarea.
% - ts: Sampling time.
% - G: Target area.
% - xobs: obstacles.
load('data_safe_landing.mat')
gbox= [63,75; -3*pi/180,0; 0,2.5];
vf= @aircraftlong;


%% Winning set, goal set
Gc= squeeze((G(:,1,:)+G(:,2,:))/2);
win= pavings(tag==1,:);
wc= [(win(:,1)+win(:,2))/2, (win(:,3)+win(:,4))/2, (win(:,5)+win(:,6))/2];
figure
hold on
plot3(Gc(1,:),Gc(2,:),Gc(3,:),'b.');
plot3(wc(:,1),wc(:,2),wc(:,3), 'g.');

% test if tag and ctlr is consistent: tag(i)=1 iff
% (any(ctlr(i))=1|pavings(i)\in G
T= xor(tag, any(ctlr,2)); % if tag(i)=any(ctlr(i)), then T(i)=0
P= pavings(find(T),:);


%% simulation
x0= [81; -pi/180; 55];

tspan= [0, ts];
x= x0;
t= 0;

tsim= t;
xsim= x';
usim= [0, 0];
while(x(1)>gbox(1,2) || x(1)<gbox(1,1) ||...
        x(2)>gbox(2,2) || x(2)<gbox(2,1) || ...
        x(3)>gbox(3,2) || x(3)<gbox(3,1) ||...
        x(1)*sin(x(2))<-0.91)  % not reach goal
    
    % compute control input
    xid= find(x(1)>=pavings(:,1) & x(1)<=pavings(:,2) & ...
        x(2)>=pavings(:,3) & x(2)<=pavings(:,4) & ...
        x(3)>=pavings(:,5) & x(3)<=pavings(:,6));
    uid= find(ctlr(xid(1),:));
    
    if (isempty(uid))
        error("Invalid controller.");
    else
%         % select a random u from all valid control values
%         pick=randperm(numel(uid));
%         u= U(uid(pick(1)),:);
        
        % select the minimum value
        uall= U(uid,:);
        [val, ind]= min(abs(uall(:,1)));
        u= uall(ind,:);
        
%         % select the first/last value
%         u= U(uid(end),:);
    end
     
    [tt, xx]= ode45(@(t,x) vf(t,x,u), tspan, x);
    x= xx(end,:)';
    t= t + tt(end,:);
    xsim= [xsim; xx(end,:)];
    tsim= [tsim; t];
    usim= [usim; u];
    
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

% % time-state curves
hf1= figure;
plot(tsim, xsim, 'LineWidth', 2);
xlabel({'$t(s)$'}, 'interpreter','latex',...
    'Units','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')
ylabel({'$x$'}, 'interpreter','latex',...
    'Units','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')
hl= legend({'$v(t)$', '$\gamma(t)$', '$h(t)$'}, 'Interpreter', 'latex',...
    'Units','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold');

% % time-control curves
usim= [usim(2:end, :); usim(end,:)];
tq= [0:0.01:tsim(end)]';
uq= interp1(tsim,usim,tq,'previous');

hf2=figure;
plot(tq, uq(:,2), 'LineWidth', 2);
axis([0, tsim(end), 0, 0.5])
xlabel({'$t(s)$'}, 'interpreter','latex',...
    'Units','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')
ylabel({'$u_2$'}, 'interpreter','latex',...
    'Units','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')
% ylabel({'$u_1,\;u_2$'}, 'interpreter','latex',...
%     'Units','points', 'FontSize',FS,...
%     'FontName','Times', 'FontWeight','bold')

% hl= legend({'$u_1(t)$','$u_2(t)$'}, 'Interpreter', 'latex',...
%     'Units','points', 'FontSize',FS,...
%     'FontName','Times', 'FontWeight','bold');
hl= legend({'angle of attack'}, 'Interpreter', 'latex',...
    'Units','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold');

% legend('wheel velocity', 'steering angle')