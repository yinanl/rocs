addpath('../../matlab')

X= [0, pi/2; -pi, pi; -2, 2; -2, 2];
U= [-0.001, 0.001; -0.001, 0.001];


%% constraints
l1 = 0.15;
l2 = 0.15;
h = 0.8*l1;
r = 0.5*l1;
a1 = atan2(h,r);
a2= asin(h/l1);

% % initialized
load('data_initialization.mat')
obs= pavings(tag==-1, :);
plot2_boxes(obs(:,1:4), [0.5,0.5,0.5], 'k', 1);
hold on
% plot2_boxes(pavings((tag==1)&any(ctlr,2),1:4), cg, cg, 1);

% % reference
N1= 50;
N2= 30;
theta11= 0:a2/N1:a2;
theta21= pi - theta11 - atan2(h-l1.*sin(theta11), l1.*cos(theta11)-r);
theta12= a2:(a1-a2)/N2:a1;
theta22= pi - theta12 + atan2(l1.*sin(theta12)-h, l1.*cos(theta12));
theta1= [theta11, theta12];
theta2= [theta21, theta22];
plot(theta1, theta2, 'LineWidth', 2)


%% The exact SCARA dynamics
% T= 0.05;
% u= [0.001; 0.001];
% x0= [0.1; 0; 0; 0];
% [tt, yy]= ode45(@(t,x) scara(t,x,u), [0, T], x0);
% plot(tt, yy(:,1:end), 'LineWidth', 1.6)


%% A simplified two double integrator model
% global ts
% ts= 0.1;
% T= 50;
% 
% u= [0;-1];
% x0= [0; 1.6264; 0; 0];
% xt= [];
% tt= [];
% x= x0;
% for k= 0: T
%     y= twoIntegrators(x, u);
%     x= y;
%     xt= [xt; y];
%     tt= [tt; ts*k];
% end
% plot(tt, xt(:,1:end), 'LineWidth', 1.6)


%% interval conversion between (x,y) and (theta1, theta2)
% x1= [0.03, 0.07]; y1= [0.18, 0.22];
% x2= [0.25, 0.28]; y2= [0.02, 0.05];
% x1_itvl= interval(x1); y1_itvl= interval(y1);
% x2_itvl= interval(x2); y2_itvl= interval(y2);
% [th1_itvl1, th2_itvl1]= xy2theta(x1_itvl, y1_itvl, l1, l2);
% [th1_itvl2, th2_itvl2]= xy2theta(x2_itvl, y2_itvl, l1, l2);

th1= interval([0.4980, 0.5772]);
th2= interval([1.5739, 1.7055]);
[x1, y1]= theta2xy(th1, th2, l1, l2);

theta1= interval([0.4903, 0.6069]);
theta2= interval([-0.9363, -0.8363]);
[x2, y2]= theta2xy(theta1, theta2, l1, l2);

