%%%
% PIPM-->PIPM: Simulation of one walking step from one initial robust margin
% set to one final robust margin set.
%%%


%% display setting
FS= 16; % fontsize
LW= 1.5; % lineweight
LThick= 1.5; % linewidth
LThin= 1;
CP = get(0, 'DefaultAxesColorOrder'); % b,o,y,p,g,c,r
CM = parula(5); % p->b->g->o->y
filterSize= 1;


%% plot specification
% % read problem data from file
setting_pipms

% % bound of disturbance
% D0= [0; 0];
% D0= [0.01; 0.02];
D0= [0.05; 0.1];
% D0= [0.07; 0.12];
% D0= [0.09; 0.16];
% D0= [0.1; 0.2];

% % import data from c++ results
% xgrid: a grid of the state space
% ugrid: a set of control values
% goalset: target set indices
% initset: initial set indices
load('pipms/data_walk_pipms_D2_004.mat');

figure
hold on

% plot initial set
idx1= initset+ones(size(initset));
initial= xgrid(idx1,:);
% plot(initial(:,1), initial(:,2), '.')
bdinit= extract_boundary2D(initial, filterSize);
fill(bdinit(:,1),bdinit(:,2),CP(6,:),'LineStyle','none');

% plot goal set
idx2= goalset+ones(size(goalset));
goal= xgrid(idx2,:);
% plot(goal(:,1), goal(:,2), '.')
bdgoal= extract_boundary2D(goal, filterSize);
fill(bdgoal(:,1),bdgoal(:,2),CP(5,:),'LineStyle','none');

% contact locations
plot(x2(1), x2(2), 'o', 'MarkerFaceColor', CP(7,:))
plot(x1(1), x1(2), 'o', 'MarkerFaceColor', CP(7,:))

% % state space: bounded by asymptotes
xpos= X(1,1):0.01:X(1,2);
yasym11= k1*(xpos-x1(1));
yasym12= -k1*(xpos-x1(1));
yasym21= k2*(xpos-x2(1));
yasym22= -k2*(xpos-x2(1));
% plot(xasym,yasym11)
% plot(xasym,yasym12)
% plot(xasym,yasym21)
% plot(xasym,yasym22)

% % robust tubes
v1u= real(sqrt(omega1^2*(xpos-x1(1)).^2+x1(2)^2+omega1^2/x1(2)^2*rbset1(1,2)));
v1l= real(sqrt(omega1^2*(xpos-x1(1)).^2+x1(2)^2+omega1^2/x1(2)^2*rbset1(1,1)));
v2u= real(sqrt(omega2^2*(xpos-x2(1)).^2+x2(2)^2+omega2^2/x2(2)^2*rbset2(1,2)));
v2l= real(sqrt(omega2^2*(xpos-x2(1)).^2+x2(2)^2+omega2^2/x2(2)^2*rbset2(1,1)));
plot(xpos, v1u, '-', 'Color', CP(5,:), 'LineWidth', LThin)
plot(xpos, v1l, '-', 'Color', CP(5,:), 'LineWidth', LThin)
plot(xpos, v2u, '-', 'Color', CP(5,:), 'LineWidth', LThin)
plot(xpos, v2l, '-', 'Color', CP(5,:), 'LineWidth', LThin)


%% plot winning set
% % load control synthesis results
% leastctlr: N x (1+M); least restrictive controller; 1st col denotes
% winning set.
% optctlr: N x 3; (time) optimal controller; winning set, control index,
% optimal value.

% % first semi-step
optctlr1(:,2)= optctlr1(:,2) + ones(size(optctlr1,1),1); % 1-based index
win1= xgrid(any(leastctlr1,2),:);
% plot(win1(:,1), win1(:,2), '.') % plot the discrete winning set 1
bdwin1= extract_boundary2D(win1, filterSize);
hfill1=fill(bdwin1(:,1),bdwin1(:,2),...
    CM(4,:),'LineStyle','--','LineWidth',LW,'EdgeColor',CM(4,:));
alpha(hfill1,0.5)

% % second semi-step
optctlr2(:,2)= optctlr2(:,2) + ones(size(optctlr2,1),1); % 1-based index
prewin2= xgrid(any(leastctlr2,2),:);
[s, z]= manifold('pipm', prewin2, omega1, x1, x1_0, zeta_0); % trim winning set 2
win2= prewin2(s<=rbset1(1,2),:);
% plot(win2(:,1), win2(:,2), '.') % plot the discrete winning set 1
bdwin2= extract_boundary2D(win2, filterSize);
hfill2=fill(bdwin2(:,1),bdwin2(:,2),...
    CM(4,:),'LineStyle','--','LineWidth',LW,'EdgeColor',CM(4,:));
alpha(hfill2,0.5)

% % plot guard set
idx3= guardset+ones(size(guardset));
guard= xgrid(idx3,:);
bdguard= extract_boundary2D(guard, filterSize);
fill(bdguard(:,1),bdguard(:,2),CP(2,:),...
    'LineStyle','-','LineWidth',LThin,'EdgeColor',CP(2,:));
% plot(guard(:,1), guard(:,2), '.k','MarkerSize',3)


%% simulation
simNum= 5;
% % time limit of each semi step
% tp1= 0.3;
% tp2= 0.8;
dt1= 2*inc_t;
dt2= 5*inc_t;
% % Q, R: the weighting factors for selecting control values
Q= 0.5;
R= 0.5;

% % settings for abnormal disturbance
Ndiver= 3;
tdiver= 6*inc_t;

D= D0;

% % intersection of winning set and initial set
xid= 1:size(xgrid,1);
wxid= xid(any(leastctlr1,2));
xinit= xgrid(intersect(idx1,wxid),:);
% % initial states are selected from a inner set of xinit
[sigma, zeta]= manifold('pipm', xinit, omega1, x1, x1_0, zeta_0);
xinit_res= xinit(sigma>(rbset1(1,1)+d_rb1(1)) & sigma<(rbset1(1,2)-d_rb1(1))...
    & zeta>(rbset1(2,1)+d_rb1(2)) & zeta<(rbset1(2,2)-d_rb1(2)),:);
% plot(xinit_res(:,1), xinit_res(:,2), '.')
% bdinitwin= extract_boundary2D(xinit, filterSize);
% fill(bdinitwin(:,1),bdinitwin(:,2),CP(6,:),'LineStyle','none');

rng('shuffle');
for iter=1:simNum
    % % initial condition
    k= randi(size(xinit_res,1));
%     k= ceil(size(xinit,1)/2);
    x0= xinit(k,:);
    
    x= x0;
    w= omega1;
    t= 0;
    d= [0;0];
    
    dsim= d';
    xsim= x0;
    usim= [];
    tsim= t;
    msim= 1;
    
    % % first semi-step: a do while loop
    dx= xgrid - repmat(x,size(xgrid,1),1);
    [v, ind]= min(dx(:,1).^2+dx(:,2).^2);
    while (~any(ismember(ind, idx3)))
        uall= find(leastctlr1(ind,2:end));
        [uv,uid]= min(Q*(ugrid(uall)-omega1).^2 + R*(ugrid(uall)-w).^2);
        w= ugrid(uall(uid));
        if(isempty(w))
            error('Not within the winning set.')
        end
        %     u= optctlr1(ind,:);
        %     if(u(3)==inf)
        %         error('Not within the winning set.')
        %     end
        %     w= ugrid(u(2),:);

        % % An abnormal trajectory: a sudden jump
        if (iter==Ndiver && abs(t-tdiver)<eps)
            D= [5; 10];  % the abnormal large disturbance
        else
            D= D0;  % the normal disturbance
        end
        
        % % generate a random disturbance
        d= D.*([-1;-1] + 2*rand(2,1));
        y= vectorfield('pipm', x1, inc_t, x, w, d);
        
        dsim= [dsim; d'];
        xsim= [xsim; y];
        usim= [usim; w];
        tsim= [tsim; t+inc_t];
        msim= [msim; 1];
        % % update
        x= xsim(end,:);
        t= t+inc_t;
        dx= xgrid - repmat(x,size(xgrid,1),1);
        [v, ind]= min(dx(:,1).^2+dx(:,2).^2);
    end
    t1= t;
    % % continue mode 1 until t>t1+dt1 or get out of guardset
    while (any(ismember(ind, idx3)) && t<t1+dt1)
        % % generate a random disturbance
        d= D.*([-1;-1] + 2*rand(2,1));
        y= vectorfield('pipm', x1, inc_t, x, w, d);
        
        dsim= [dsim; d'];
        xsim= [xsim; y];
        usim= [usim; w];
        tsim= [tsim; t+inc_t];
        msim= [msim; 1];
        
        % % update
        x= xsim(end,:);
        t= t+inc_t;
        dx= xgrid - repmat(x,size(xgrid,1),1);
        [v, ind]= min(dx(:,1).^2+dx(:,2).^2);
    end
    if (~any(ismember(ind, idx3)))
        x= xsim(end-1,:);
        t= tsim(end-1,:);
        dsim= dsim(1:end-1,:);
        xsim= xsim(1:end-1,:);
        usim= usim(1:end-1,:);
        tsim= tsim(1:end-1,:);
        msim= msim(1:end-1,:);
        t= t-inc_t;
    end
    
    % % second semi-step
    [sigma, zeta]= manifold('pipm', x, omega2, x2, x2_0, zeta_0);
    while (sigma > rbset2(1,2) || sigma < rbset2(1,1) ||...
            zeta > rbset2(2,2) || zeta < rbset2(2,1))
        % % compute control input
        dx= xgrid - repmat(x,size(xgrid,1),1);
        [v, ind]= min(dx(:,1).^2+dx(:,2).^2);
        uall= find(leastctlr2(ind,2:end));
        [uv,uid]= min(Q*(ugrid(uall)-omega2).^2 + R*(ugrid(uall)-w).^2);
        w= ugrid(uall(uid));
        if(isempty(w))
            error('Not within the winning set.')
        end
        %     u= optctlr2(ind,:);
        %     if(u(3)==inf)
        %         error('Not within the winning set.')
        %     end
        %     w= ugrid(u(2),:);
        
        % % generate a random disturbance
        d= D.*([-1;-1] + 2*rand(2,1));
        y= vectorfield('pipm', x2, inc_t, x, w, d);

        dsim= [dsim; d'];
        xsim= [xsim; y];
        usim= [usim; w];
        tsim= [tsim; t+inc_t];
        msim= [msim; 2];
        
        % % update
        x= xsim(end,:);
        t= t+inc_t;
        [sigma, zeta]= manifold('pipm', x, omega2, x2, x2_0, zeta_0);
    end
    t2= t;
    % % continue the last control config for ~0.05s
    while (sigma <= rbset2(1,2) && sigma >= rbset2(1,1) &&...
            zeta <= rbset2(2,2) && zeta >= rbset2(2,1) &&...
            t<t2+dt2)
        
        % % generate a random disturbance
        d= D.*([-1;-1] + 2*rand(2,1));
        y= vectorfield('pipm', x2, inc_t, x, w, d);
        
        dsim= [dsim;d'];
        xsim= [xsim; y];
        usim= [usim; w];
        tsim= [tsim; t+inc_t];
        msim= [msim; 2];
        
        % % update
        x= xsim(end,:);
        t= t+inc_t;
        [sigma, zeta]= manifold('pipm', x, omega2, x2, x2_0, zeta_0);
        dx= xgrid - repmat(x,size(xgrid,1),1);
        [v, ind]= min(dx(:,1).^2+dx(:,2).^2);
    end
    if (sigma > rbset2(1,2) || sigma < rbset2(1,1) ||...
            zeta > rbset2(2,2) || zeta < rbset2(2,1))
        x= xsim(end-1,:);
        t= tsim(end-1,:);
        dsim= dsim(1:end-1,:);
        xsim= xsim(1:end-1,:);
        usim= usim(1:end-1,:);
        tsim= tsim(1:end-1,:);
        msim= msim(1:end-1,:);
    end
    if (iter == Ndiver)
        plot(xsim(:,1), xsim(:,2),'-b','LineWidth', LW)
    else
        plot(xsim(:,1), xsim(:,2),'-k','LineWidth', LW)
    end
end


% % format display
axis square
rectangle('Position', [X(1,1) X(2,1) X(1,2)-X(1,1) X(2,2)-X(2,1)])
axis([-0.2 0.8 0.1 1.2])
xlabel({'$x$[m]'}, 'interpreter','latex',...
    'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')

ylabel({'$\dot{x}$[m/s]'}, 'interpreter','latex',...
    'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')
