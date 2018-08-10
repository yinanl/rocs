%%%
% Plot reachable sets of test points for reversed VDP
%
% using interval solver
%%%

addpath('../../matlab')

% % define color
cr= [0.6350 0.0780 0.1840];
cb= [0 0.4470 0.7410];
cy= [0.9290 0.6940 0.1250];
gray = [0.6,0.6,0.6];
gold = [1,0.84,0];
lightgold = [255,236,139]/255;
orange= [0.8500, 0.3250, 0.0980];
pink = [255,182,193]/255;
orchid= [218,112,214]/255;
lightblue = [176, 226, 255]/255;
slateblue= [123,104,238]/255;
turquo= [64,224,208]/255;
green= [0.4660, 0.6740, 0.1880]*0.7;

FS= 12; % FontSize
LW= 1.5; % LineWidth
MS= 2; % MarkerSize


%% test points
O= [-3 3; -3 3]; % state space
h= 0.05; % time of integration
delta= 10; % bound of disturbance

vf= @(t,x) [-x(2); x(1) + (x(1)^2-1)*x(2)];

xc01= [0.5;0.3]; r1= 0.015/2; %r1= 0.0002/2;
xc02= [-1.2; -0.6]; r2= 0.0035/2; %r2= 0.0008/2;
xc03= [-2.3; 1.7]; r3= 0.0002/2; %r3= 0.0002/2
xc0= [xc01 xc02 xc03];
r= [r1;r2;r3];

x01= interval([xc01(1)-r1, xc01(1)+r1; xc01(2)-r1, xc01(2)+r1])';
x02= interval([xc02(1)-r2, xc02(1)+r2; xc02(2)-r2, xc02(2)+r2])';
x03= interval([xc03(1)-r3, xc03(1)+r3; xc03(2)-r3, xc03(2)+r3])';

x0= [x01; x02; x03];

% % initial intervals
% xc0= [0.5;0.3];
% 
% r1= 0.001/2;
% x01= interval([xc0(1)-r1, xc0(1)+r1; xc0(2)-r1, xc0(2)+r1])';
% 
% r2= 0.1/2;
% x02= interval([xc0(1)-r2, xc0(1)+r2; xc0(2)-r2, xc0(2)+r2])';
% 
% x0= [x01; x02];
% % x0= x01;

% conditions for k and eps
alpha= 0.5;
K= 5;
kbar= 4;
orderk= log(alpha*delta/K)+log(factorial(kbar))/log(h);
wmin= alpha*delta*h/(K*exp(h));


%% compare reachable sets
% % the approximated reachable set
xtTaylor1= [0.47677391, 0.49256354; 0.30469853, 0.32130625]; %[0.484563, 0.484774; 0.312894, 0.313115]
xtTaylor2= [-1.17003353, -1.16634472; -0.67403922, -0.66998864]; %[-1.16861090, -1.16776737; -0.67247739, -0.67154956];
xtTaylor3= [-2.39211797, -2.39189747; 1.99699834, 1.99746280]; %[-2.39211797, -2.39189747; 1.99699834, 1.99746280];

xtTaylor= [xtTaylor1;xtTaylor2;xtTaylor3];

for p= 1:3
    figure(p)
    hold on
    % the nominal reachable set
    w= 2*r(p)/10; 
    [X,Y]= meshgrid(xc0(1,p)-r(p):w:xc0(1,p)+r(p), xc0(2,p)-r(p):w:xc0(2,p)+r(p));
    Nx= size(X,1);
    Ny= size(X,2);

    xt= zeros(Nx*Ny,2);
    for i= 1:Nx
        for j= 1:Ny
            % the nominal system
            [tn,yn]= ode45(vf, [0 h], [X(i,j);Y(i,j)]);
            plot(yn(end,1),yn(end,2),'.','MarkerSize',5)
            xt(j+(i-1)*Nx,:)= yn(end,:);
    %         % the perturbed system
    %         for k= 1:30
    %             [td,yd]= ode45(@(t,x) vdp_disturbed(t,x,delta), [0 h], [X(i,j);Y(i,j)]);
    %             plot(yd(end,1),yd(end,2),'.')
    %         end
        end
    end
    
    vhat= (xtTaylor(2*(p-1)+1,2)-xtTaylor(2*(p-1)+1,1))*(xtTaylor(2*(p-1)+2,2)-xtTaylor(2*(p-1)+2,1));

    [verid,v0]= convhull(xt(:,1),xt(:,2));
    Pxt= Polyhedron('V',xt(verid,:));
    vhat/v0;

    % the perturbed reachable set
    % Pw= Polyhedron('V', delta*h*[-0.5 -0.5;0.5 -0.5;0.5 0.5;-0.5 0.5]);
    Pw= Polyhedron('V', delta*h*[-1 -1;1 -1;1 1;-1 1]);
    P= Pxt+Pw;
    v= P.volume;
    vhat/v;


    % % display
    xbl= xtTaylor(2*(p-1)+1,1);
    ybl= xtTaylor(2*(p-1)+2,1);
    xwid= xtTaylor(2*(p-1)+1,2)-xtTaylor(2*(p-1)+1,1);
    ywid= xtTaylor(2*(p-1)+2,2)-xtTaylor(2*(p-1)+2,1);
    rectangle('Position', [xbl ybl xwid ywid],...
        'LineWidth', 1, 'LineStyle', '-', 'FaceColor',gray);
    hold on

%     P.plot('Color', gray, 'Alpha', 0.3)
    Pxt.plot('Color', turquo)
    
    % % format display
    axis equal
    axis([xbl-xwid/2 xbl+xwid*1.5 ybl-ywid/2 ybl+ywid*1.5])
    xlabel({'$x_1$'}, 'interpreter','latex',...
        'FontUnits','points', 'FontSize',16,...
        'FontName','Times', 'FontWeight','bold')
    
    ylabel({'$x_2$'}, 'interpreter','latex',...
        'FontUnits','points', 'FontSize',16,...
        'FontName','Times', 'FontWeight','bold')

end


%% plot reachable sets in vector field
figure(4)
hold on

% % display vector field
N= 20;
w= (O(:,2)-O(:,1))/N;
fct= @(t,x) ([-x(2); x(1) + (x(1)^2-1)*x(2)]);
plot2_vf(O, w, {fct}, 3)
rectangle('Position', [O(1,1) O(2,1) O(1,2)-O(1,1) O(2,2)-O(2,1)],...
    'LineWidth', 1, 'LineStyle', '-');

% % display the limit cycle
load('data_roavdp_real.mat')
[t0,r0]=cart2pol(bd0(:,1),bd0(:,2));
[tt0, order]= sort(t0);
[xlc,ylc]= pol2cart(tt0,r0(order));
h0= plot(xlc,ylc, 'LineWidth', LW, 'Color', 'k');

% % display test points
% xaprior= [interval([0.48413, 0.500512; 0.298752, 0.314302])';...
%     interval([0.432022, 0.550026; 0.248998, 0.368803])'];
% plot2_boxes(xaprior, cy, cy, 0.5)

plot2_boxes(x0, cb, cb, 0.7)

xh= [interval(xtTaylor1)'; interval(xtTaylor2)'; interval(xtTaylor3)'];
plot2_boxes(xh, cr, cr, 0.7)

xhc= mid(xh);
xc= mid(x0);
dcu= xhc(:,1) - xc(:,1);
dcv= xhc(:,2) - xc(:,2);
q= quiver(xc(:,1), xc(:,2), dcu, dcv, 0);
q.Color= 'k';
q.LineWidth= 2;

% % % display solution by ode integration (reference)
% [tt,yy]= ode45(fct, [0 h], xc0);
% plot(yy(:,1),yy(:,2),'.-')

% % format display
axis([O(1,1),O(1,2),O(2,1),O(2,2)])
xlabel({'$x_1$'}, 'interpreter','latex',...
    'FontUnits','points', 'FontSize',16,...
    'FontName','Times', 'FontWeight','bold')

ylabel({'$x_2$'}, 'interpreter','latex',...
    'FontUnits','points', 'FontSize',16,...
    'FontName','Times', 'FontWeight','bold')

% % annotation
a1= annotation('ellipse',[.57 .55 .03 .03], 'Color', green);
t1= annotation('textbox',[.6 .55 .05 .05],'String','$x_1$','FitBoxToText','on');
t1.EdgeColor='none';
t1.Interpreter= 'latex';
t1.FontSize= 16;
t1.Color= green;

a2= annotation('ellipse',[.35 .41 .035 .035], 'Color', green);
t2= annotation('textbox',[.4 .41 .05 .05],'String','$x_2$','FitBoxToText','on');
t2.EdgeColor='none';
t2.Interpreter= 'latex';
t2.FontSize= 16;
t2.Color= green;

a3= annotation('ellipse',[.2 .73 .035 .07], 'Color', green);
t3= annotation('textbox',[.23 .7 .05 .05],'String','$x_3$','FitBoxToText','on');
t3.EdgeColor='none';
t3.Interpreter= 'latex';
t3.FontSize= 16;
t3.Color= green;

% % save
% print('test_reach_vdp', '-depsc')


%% compare the trajectories of perturbed and nominal systems
% [tn,yn]= ode45(vdp_nominal, [0 h], xc0);
% [td,yd]= ode45(@(t,x) vdp_disturbed(t,x,delta), [0 h], xc0);
% 
% tt= 0:0.001:h;
% yy=[interp1(tn,yn(:,1),tt'), interp1(td,yd(:,1),tt'), ...
%     interp1(tn,yn(:,2),tt'), interp1(td,yd(:,2),tt')];
% % yy=[interp1(tn,yn(:,1),tt'), interp1(td,yd(:,1),tt')];
% figure
% plot(tt,yy)
