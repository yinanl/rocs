
%% test continuous-time reachable set

% state and input spaces
% X= [-0.3 0.3;-0.1 0.1];
% U= -0.1:0.02:0.1;
X= [-8 8;-8 8];  % it includes [-2pi,2pi]
U= -10:0.5:10;
O= [-0.05 0.05; -0.01 0.01];

% sampling time
tau= 0.01;


%% plot vector field
addpath('../../matlab')
% % define color
cr= [0.6350 0.0780 0.1840];
cb= [0 0.4470 0.7410];
cy= [0.9290 0.6940 0.1250];


% % display vector field
figure
hold on

N= 20;
w= (X(:,2)-X(:,1))/N;

% plot2_vf(X, w, {@ipdlu}, 1, U(35))
plot2_vf(X, w, {@ipdlu}, 1)
rectangle('Position', [X(1,1) X(2,1) X(1,2)-X(1,1) X(2,2)-X(2,1)],...
    'LineWidth', 1, 'LineStyle', '-');
rectangle('Position', [O(1,1) O(2,1) O(1,2)-O(1,1) O(2,2)-O(2,1)],...
    'LineWidth', 1, 'LineStyle', '-');


%% test points
% some arbitrary test points
x01= interval([0.02, 0.03; -0.02, -0.01])';
x02= interval([0.1, 0.15; -0.01, 0.01])';
x03= interval([-0.27, -0.21; 0.03, 0.04])';
x04= interval([-0.1, -0.01; -0.06, -0.02])';
x05= interval([-0.05, -0.025; -0.01 -0.005])';

% the point to be compared with analytical bounds
eta= [0.001;0.001];
x0= interval([0.05; -0.01]-eta/2, [0.05; -0.01]+eta/2)';
% x0= interval([-0.2;-0.006]-eta/2, [-0.2;-0.006]+eta/2)';
% x0= x01;

nbox= size(x0,1);

% plot test points
plot2_boxes(x0, cb, cb, 0.7)


%% simulation + analytical bounds
% % display point simulation result
xtu = [];
m= 0.2; % mass of the pendulum
b= 0.1; % coefficient of friction of cart
J= 0.006; % inertia
g= 9.8; % gravitational constant
l= 0.3; % pendulum length

a1= (m*g*l)/(J+m*l^2);
a2= b/(J+m*l^2);
a3= l/(J+m*l^2);

for i= 1:nbox
    for j= 1:size(U,2)
        
        [~, xt]= ode45(@(t,x,u) ipdlu(t,x,U(j)),[0 tau], mid(x0(i,:)));
        xc= mid(x0(i,:));
        L= [0 1;sqrt(a1^2 + a3^2*U(j)^2) -a2];
        beta= expm(L*tau) * width(x0)'/2;
        xtab= interval([xt(end,:)'-beta, xt(end,:)'+beta])';
        
        xtu = cat(1, xtu, xtab);
        
        plot2_boxes(xtab, cy, cy, 0.7)
        
        dcu= xt(end,1) - xc(1);
        dcv= xt(end,2) - xc(2);
        q= quiver(xc(1), xc(2), dcu, dcv, 0);
        q.Color= 'k';
        q.LineWidth= 2;
        
    end
end

axis([X(1,1),X(1,2),X(2,1),X(2,2)])


%% validated solutions
% 
% % h0= tau * ones(nbox,1);  % the maximal time step
% % htol= 1e-3 * ones(nbox, 1);  % time step precision
% % hmin= htol;
% alpha= 2;
% tol_a= 1e-3;
% tol_r= 1e-1;
% 
% tc= {@tc1_ipdl, @tc2_ipdl};
% jactc= {@jac_tc1_ipdl};
% 
% for j= 8:8 %size(u0,2)
% %     [xh, xall, h]= vnode_solver_hoe2(x0, u0(j), tc, jactc, h0, htol, hmin, alpha);
%     [xh, xhis, time]= vnode_solver(tau, x0, u0(j), tc, jactc, tol_a, tol_r, alpha);
%     plot2_boxes(xh, cr, cr, 0.5)
% end
% 
% % axis([O(1,1),O(1,2),O(2,1),O(2,2)])
% xlabel({'$x_1$'}, 'interpreter','latex',...
%     'FontUnits','points', 'FontSize',16,...
%     'FontName','Times', 'FontWeight','bold')
% 
% ylabel({'$x_2$'}, 'interpreter','latex',...
%     'FontUnits','points', 'FontSize',16,...
%     'FontName','Times', 'FontWeight','bold')


%% save
% print('test_reach_vdp', '-depsc')