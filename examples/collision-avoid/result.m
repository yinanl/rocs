
addpath('../../matlab')
%% load problem & controller
% load('data_safeset_d0.8.mat')

% filename= 'data_invset_0.8.h5';
filename= "controller_safety_itvl_0.8-1.2-0.3-0.1.h5";
G= h5read(filename, '/G');
X= h5read(filename, '/X')';
U= h5read(filename, '/U')';
pavings= h5read(filename, '/pavings')';
tag= h5read(filename, '/tag');
ctlr= h5read(filename, '/ctlr')';
ts= h5read(filename, '/ts');
D= [-0.8 0.8; -0.8 0.8];
Dr= D(:,2)-D(:,1);
Dm= (D(:,1)+D(:,2))/2;

% safeset= zeros(size(G,3), 4);
% for i= 1:size(G,3)
%     safeset(i,:)= [G(1,:,i), G(2,:,i)];
% end
% winid= find(any(ctlr,2));
% winset= pavings(winid,:);
winset= pavings(tag==1, :);
wctr= [(winset(:,1)+winset(:,2))/2,...
    (winset(:,3)+winset(:,4))/2,...
    (winset(:,5)+winset(:,6))/2];

% plot(G(:,1), G(:,2), '.', 'MarkerSize', 10)
% axis('equal')

%% simulation
% x0= [-2; 2.4; -pi/4];
x0= [-1.3; 1.3; 0];
% x0= [0.6; -2.1; 4/pi];
Tsim= 50;

tsim= [];
usim= [];
dsim= [];
xsim= [];
xts= x0';
x= x0;
t= 0;

rng('shuffle')
while(t<Tsim)
    
    % compute control input
    if(x(3)>pi)
        x(3)= x(3) - 2*pi;
    end
    if(x(3)<-pi)
        x(3)= 2*pi+x(3);
    end
    xid= find(x(1)>=pavings(:,1) & x(1)<=pavings(:,2) & ...
        x(2)>=pavings(:,3) & x(2)<=pavings(:,4) & ...
        x(3)>=pavings(:,5) & x(3)<=pavings(:,6)); % direct search
    uid= find(ctlr(xid(1),:));
    % select a random u from all valid control values
    pick=randperm(numel(uid));
    u= U(uid(pick(1)),:);
    d= Dr.*rand(2, 1)-Dm;
    [tt,xx]= ode45(@(t,x) ca(t,x,u,d), [0 ts], x);
    
    % append state for simulation
    tsim= cat(1, tsim, t+tt);
    xsim= cat(1, xsim, xx);
    usim= cat(1, usim, repmat(u,size(tt,2),1));
    dsim= cat(1, dsim, repmat(d,size(tt,2),1));
    xts= cat(1, xts, xx(end,:));
    
    % move to the next step
    x= xx(end,:)';
    t= t + tt(end);
end

 
%% display
% define color
cr= [0.6350 0.0780 0.1840];
cb= [0 0.4470 0.7410];
cy= [0.9290 0.6940 0.1250];

FS= 16; % fontsize
LW= 1.5; % lineweight

% plot the initial safe set
% hf2= figure;
% plot2_boxes(safeset, [0.5,0.5,0.5], 'k', 1);
% hold on
% rectangle('Position',[X(1,1),X(2,1),X(1,2)-X(1,1),X(2,2)-X(2,1)],...
%     'LineWidth',LW, 'LineStyle', '-')
% axis([X(1,1) X(1,2) X(2,1) X(2,2)])
% xlabel({'$x_1$'}, 'interpreter','latex',...
%     'FontSize',FS, 'FontName','Times', 'FontWeight','bold')
% ylabel({'$x_2$'}, 'interpreter','latex',...
%     'FontSize',FS, 'FontName','Times', 'FontWeight','bold')


% plot the controlled invariant set
t=0: 0.01:2*pi;
r= 1.2;
x= r*sin(t);
y= r*cos(t);

figure
plot(x,y, 'LineWidth', LW)
hold on
plot2_boxes(winset(:,1:4), [0.5,0.5,0.5], 'k', 1);
rectangle('Position',[X(1,1),X(2,1),X(1,2)-X(1,1),X(2,2)-X(2,1)],...
    'LineWidth',LW, 'LineStyle', '-')
p1= plot(xsim(:,1), xsim(:,2), '-', 'LineWidth',LW);
p2= plot(xts(:,1), xts(:,2), '.', 'LineWidth',LW);
p1.Color= [39,64,139]/255; 
p2.Color= [189,183,107]/255;

axis([X(1,1) X(1,2) X(2,1) X(2,2)])
xlabel({'$x_r$'}, 'interpreter','latex',...
    'FontSize',FS, 'FontName','Times', 'FontWeight','bold')
ylabel({'$y_r$'}, 'interpreter','latex',...
    'FontSize',FS, 'FontName','Times', 'FontWeight','bold')
% plot3(wctr(:,1), wctr(:,2), wctr(:,3), '.', 'MarkerSize', 10)