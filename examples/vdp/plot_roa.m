%%%
% Plot the region of attraction (ROA) of reversed VDP
%
% using interval solver
%%%

addpath('../../matlab')
clear all
clc

cr= [0.6350 0.0780 0.1840];
cb= [0 0.4470 0.7410];
cy= [0.9290 0.6940 0.1250];

gold= [255, 215, 0]/255;
turquo= [64,224,208]/255;
orchid= [218,112,214]/255;
slateblue= [123,104,238]/255;

FS= 12; % FontSize
LW= 1.5; % LineWidth
MS= 2; % MarkerSize


%% plot roa boundary
load('data_roavdp_real.mat')
load('data_roavdp.mat')

hf1= figure;
hold on

% % plot state space
rectangle('position',[-4,-4,8,8]);

% % plot vector field
N= 10;
w= (X(:,2)-X(:,1))/N;
plot2_vf(X, w, {@(t,x) [-x(2); x(1) + (x(1)^2-1)*x(2)]}, 2)

% % plot boundaries
[t0,r0]=cart2pol(bd0(:,1),bd0(:,2));
[tt0, order]= sort(t0);
[x0,y0]= pol2cart(tt0,r0(order));
h0= plot(x0,y0, 'LineWidth', LW, 'Color', 'k');

[t1,r1]=cart2pol(bd01(:,1),bd01(:,2));
[tt1, order]= sort(t1);
[x1,y1]= pol2cart(tt1,r1(order));
h1= plot(x1,y1,'-', 'LineWidth', LW, 'Color', gold);

[t2,r2]=cart2pol(bd02(:,1),bd02(:,2));
[tt2, order]= sort(t2);
[x2,y2]= pol2cart(tt2,r2(order));
h2= plot(x2,y2,'-', 'LineWidth', LW, 'Color', turquo);

[t3,r3]=cart2pol(bd03(:,1),bd03(:,2));
[tt3, order]= sort(t3);
[x3,y3]= pol2cart(tt3,r3(order));
h3= plot(x3,y3,'-', 'LineWidth', LW, 'Color', slateblue);


% % plot roa core
% P= [1.5 -0.5; -0.5 1];
% h4= plot2_ellipse(P, 1.43, cr, '--', LW);

% axis equal
axis([X(1,1) X(1,2) X(2,1) X(2,2)])

xlabel({'$x_1$'}, 'interpreter','latex',...
    'FontSize',FS, 'FontName','Times', 'FontWeight','bold')

ylabel({'$x_2$'}, 'interpreter','latex',...
    'FontSize',FS, 'FontName','Times', 'FontWeight','bold')

hl= legend([h0,h1,h2,h3],{'Limit cycle','$\varepsilon=0.03$',...
    '$\varepsilon=0.01$','$\varepsilon=0.005$'}, 'Interpreter', 'latex',...
    'FontSize',14, 'FontName','Times', 'FontWeight','bold');

% hl= legend([h0,h1,h2,h3,h4],{'Limit cycle','$\varepsilon=0.03$',...
%     '$\varepsilon=0.01$','$\varepsilon=0.005$', '$\Omega_{1.43}$'}, 'Interpreter', 'latex',...
%     'FontSize',14,...
%     'FontName','Times', 'FontWeight','bold');

print(hf1, 'compare_roas.eps', '-depsc2')


%% diplay winning set
% wbox= pavings(tag==1,:);
% w= wbox(1:5:end,:);
% % w= wbox;
% ct= [(w(:,1)+w(:,2))/2, (w(:,3)+w(:,4))/2];
% 
% figure
% hold on
% 
% for i= 1:size(ct,1)
%     plot(ct(i,1), ct(i,2), '.', 'markersize', 12) 
% end
% % plot2_boxes(w)
% 
% % axis([X(1,:) X(2,:)])


%% simulate test points
delta= 5;
tsim='100';

x01=[1.2; 1.6]; % >0.03
x02= [0.7; -1.2]; % (0.01, 0.03)
x03= [-1.7; 0.53]; % (0.005, 0.01)
x04= [-2.0; 0.1]; % <0.005

% % % simulate inside matlab using the perturbed system model 
% % % (slow, better to use simulink)
% for k= 1:1
%     [td,yd]= ode45(@(t,x) vdp_disturbed(t,x,delta), [0 tsim], [0.9;2.6]);
%     plot(yd(:,1),yd(:,2))
% end

% % call simulink
for l=1:10
    s= rng('shuffle');
    seed= double(s.Seed);
    x0= x02;
    simout=sim('sim_vdp_delta', 'StopTime', tsim,...
        'SaveState','on','StateSaveName','xoutNew',...
        'SaveOutput','on','OutputSaveName','youtNew');
end
                 
t=simout.get('tout');
xout=simout.get('xoutNew');
% plot(xout(:,1), xout(:,2))
