%%%
% PIPM-->PPM: Simulation of one walking step from one initial robust margin
% set to one final robust margin set.
%%%


%% display setting
FS= 16; % fontsize
LThick= 1.5; % linewidth
LThin= 1;
CP = get(0, 'DefaultAxesColorOrder'); % b,o,y,p,g,c,r
CM = parula(5); % p->b->g->o->y
filterSize= 1;


%% import data from c++ results
setting_pipm2ppm

% % bound of disturbance
% D= [0; 0];
% D= [0.01; 0.02];
D= [0.05; 0.1];
% D= [0.07; 0.12];
% D= [0.09; 0.16];
% D= [0.1; 0.2];
% D= [0.15; 0.3];
% D= [0.2; 0.4];

% % xgrid: a grid of the state space
% % ugrid: a set of control values
% % goalset: target set indices
% % initset: initial set indices
load('pipm2ppm/data_walk_pipm2ppm_D5_004.mat')


%% plot specification
figure
axBase= axes;
% axSub= axes('Position',[0.55,0.3,0.2,0.2]);
hold(axBase, 'on')
% hold(axSub, 'on')

% % plot initial set
idx1= initset+ones(size(initset));
initial= xgrid(idx1,:);
% plot(initial(:,1), initial(:,2), '.')
bdinit= extract_boundary2D(initial, filterSize);
fill(bdinit(:,1),bdinit(:,2),CP(6,:),...
    'LineStyle','none', 'Parent', axBase);

% % plot goal set
idx2= goalset+ones(size(goalset));
goal= xgrid(idx2,:);
% plot(goal(:,1), goal(:,2), '.')
bdgoal= extract_boundary2D(goal, filterSize);
fill(bdgoal(:,1),bdgoal(:,2),CP(5,:),...
    'LineStyle','none', 'Parent', axBase);

% % plot guard set
idx3= guardset+ones(size(guardset));
guard= xgrid(idx3,:);
bdguard= extract_boundary2D(guard, filterSize);
% fill(bdguard(:,1),bdguard(:,2),CP(2,:),...
%     'LineStyle','-','LineWidth',LThin,'EdgeColor',CP(2,:),...
%     'Parent', axBase);
% plot(axBase, guard(:,1), guard(:,2), '.k','MarkerSize',3)
[sg1, zg1]= manifold('pipm', guard, omega1, x1, x1_0, zeta_0);
[sg2, zg2]= manifold('ppm', guard, omega2, x2, x2_0, zeta_0);
zg= [zg1 zg2];

% % plot contact locations
plot(x2(1), x2(2), 'o', 'MarkerFaceColor', CP(7,:), 'Parent', axBase)
plot(x1(1), x1(2), 'o', 'MarkerFaceColor', CP(7,:), 'Parent', axBase)

% % plot robust tubes
xpos= (X(1,1):0.001:X(1,2))';
v1l= real(sqrt(omega1^2*(xpos-x1(1)).^2 + x1(2)^2 + omega1^2/x1(2)^2*q_pipm(1,1)));
v1u= real(sqrt(omega1^2*(xpos-x1(1)).^2 + x1(2)^2 + omega1^2/x1(2)^2*q_pipm(1,2)));
v2l= real(sqrt(-omega2^2*(xpos-x2(1)).^2 + x2(2)^2 - omega2^2/x2(2)^2*q_ppm(1,2)));
v2u= real(sqrt(-omega2^2*(xpos-x2(1)).^2 + x2(2)^2 - omega2^2/x2(2)^2*q_ppm(1,1)));
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
% optctlr1(:,2)= optctlr1(:,2) + ones(size(optctlr1,1),1); % 1-based index
win1= xgrid(any(leastctlr1,2),:);
% plot(win1(:,1), win1(:,2), '.') % plot the discrete winning set 1
bdwin1= extract_boundary2D(win1, filterSize);
hfill1=fill(bdwin1(:,1),bdwin1(:,2),...
    CM(4,:),'LineStyle','--','LineWidth',LThin,'EdgeColor',CM(4,:));
alpha(hfill1,0.5)

% % second semi-step
% optctlr2(:,2)= optctlr2(:,2) + ones(size(optctlr2,1),1); % 1-based index
% prewin2= xgrid(any(leastctlr2,2),:);
% [s, z]= manifold('ppm', prewin2, omega1, x1, x1_0, zeta_0); % trim winning set 2
% win2= prewin2(s<=rbset1(1,2),:);
win2= xgrid(any(leastctlr2,2),:);
% plot(win2(:,1), win2(:,2), '.') % plot the discrete winning set 1
bdwin2= extract_boundary2D(win2, filterSize);
hfill2=fill(bdwin2(:,1),bdwin2(:,2),...
    CM(4,:),'LineStyle','--','LineWidth',LThin,'EdgeColor',CM(4,:));
alpha(hfill2,0.5)

% % merged winning set for one step
% lowerN= size(bdwin1,1)/2;
% bdwin= [bdwin1(1:lowerN,:); bdwin2(lowerN+1:end,:)];
% hfill=fill(bdwin(:,1),bdwin(:,2),CM(4,:),...
%     'LineStyle','--','LineWidth',LThin,'EdgeColor',CM(4,:),...
%     'Parent', axBase);
% alpha(hfill,0.5)

% % plot guard set
fill(bdguard(:,1),bdguard(:,2),CP(2,:),...
    'LineStyle','-','LineWidth',LThin,'EdgeColor',CP(2,:),...
    'Parent', axBase);
plot(axBase, guard(:,1), guard(:,2), '.k','MarkerSize',3)

% % plot embedded figure
% hfill=fill(bdwin(:,1),bdwin(:,2),CM(4,:),...
%     'LineStyle','--','LineWidth',LThin,'EdgeColor',CM(4,:),...
%     'Parent', axSub);
% alpha(hfill,0.5)
% fill(bdguard(:,1),bdguard(:,2),CP(2,:),'LineStyle','-',...
%     'LineWidth',LThin,'EdgeColor',CP(2,:), 'Parent', axSub);


%% simulation
simNum= 5;
scaleD= 1;

% time limit of each semi step
tp1= 0.25;
tp2= 0.8;

% Q, R: the weighting factors for selecting control values
Q= 0.5;
R= 0.5;

% intersection of winning set and initial set
xid= 1:size(xgrid,1);
wxid= xid(any(leastctlr1,2));
xinit= xgrid(intersect(idx1,wxid),:);
% plot(xinit(:,1), xinit(:,2), '.')
% bdinitwin= extract_boundary2D(xinit, filterSize);
% fill(bdinitwin(:,1),bdinitwin(:,2),CP(6,:),'LineStyle','none');


rng('shuffle');
for iter=1:simNum
    
    dsim= [];
    xsim= [];
    usim= [];
    tsim= [];
    msim= [];
    
    % % initial condition
    k= randi(size(xinit,1));
%     k= ceil(size(xinit,1)/2);
    x0= xinit(k,:);
    x= x0;
    w= omega1;
    t= 0;
    
    % % first semi-step: a do while loop
    ind= get_stateid(x, xgrid);
    while (~any(ismember(ind, idx3)))
        % % compute control input
        uall= find(leastctlr1(ind,2:end));
        [uv,uid]= min(Q*(ugrid(uall)-omega1).^2 + R*(ugrid(uall)-w).^2);
        w= ugrid(uall(uid));
        if(isempty(w))
            plot(xsim(:,1), xsim(:,2),'-k','LineWidth', LThin, 'Parent', axBase)
            error('Not within the winning set.')
        end
        %     u= optctlr1(ind,:);
        %     if(u(3)==inf)
        %         error('Not within the winning set.')
        %     end
        %     w= ugrid(u(2),:);
        
        % % generate a random disturbance
        d= scaleD*D.*([-1;-1] + 2*rand(2,1));
        % % record
        dsim= [dsim; d'];
        xsim= [xsim; x];
        usim= [usim; w];
        tsim= [tsim; t];
        msim= [msim; 1];
        % % update
        x= vectorfield('pipm', x1, inc_t, x, w, d);
        t= t+inc_t;
        % % check update state
        ind= get_stateid(x, xgrid);
    end
    % % continue mode 1 until t>tp1 or get out of guardset
    while (any(ismember(ind, idx3)) && t<tp1)
        % % generate a random disturbance
        d= scaleD*D.*([-1;-1] + 2*rand(2,1));
        dsim= [dsim; d'];
        xsim= [xsim; x];
        usim= [usim; w];
        tsim= [tsim; t];
        msim= [msim; 1];
        % % update
        x= vectorfield('pipm', x1, inc_t, x, w, d);
        t= t+inc_t;
        ind= get_stateid(x, xgrid);
    end
    if (~any(ismember(ind, idx3)))
        x= xsim(end,:);
        t= t-inc_t;
    end
    
    % % second semi-step
    [sigma, zeta]= manifold('ppm', x, omega2, x2, x2_0, zeta_0);
    while (sigma > q_ppm(1,2) || sigma < q_ppm(1,1) ||...
            zeta > q_ppm(2,2) || zeta < q_ppm(2,1))
        ind= get_stateid(x, xgrid);
        % % compute control input
        uall= find(leastctlr2(ind,2:end));
        [uv,uid]= min(Q*(ugrid(uall)-omega2).^2 + R*(ugrid(uall)-w).^2);
        w= ugrid(uall(uid));
        if(isempty(w))
            plot(xsim(:,1), xsim(:,2),'-k','LineWidth', LThin, 'Parent', axBase)
            error('Not within the winning set.')
        end
        %     u= optctlr2(ind,:);
        %     if(u(3)==inf)
        %         error('Not within the winning set.')
        %     end
        %     w= ugrid(u(2),:);
        
        % % generate a random disturbance
        d= scaleD*D.*([-1;-1] + 2*rand(2,1));
        dsim= [dsim; d'];
        xsim= [xsim; x];
        usim= [usim; w];
        tsim= [tsim; t];
        msim= [msim; 2];
        
        % % update
        x= vectorfield('ppm', x2, inc_t, x, w, d);
        t= t+inc_t;
        ind= get_stateid(x, xgrid);
        [sigma, zeta]= manifold('ppm', x, omega2, x2, x2_0, zeta_0);
    end
    % % continue the last control config for ~0.05s
    while (sigma <= q_ppm(1,2) && sigma >= q_ppm(1,1) &&...
            zeta <= q_ppm(2,2) && zeta >= q_ppm(2,1) &&...
            t<tp2)
        % % generate a random disturbance
        d= scaleD*D.*([-1;-1] + 2*rand(2,1));
        dsim= [dsim; d'];
        xsim= [xsim; x];
        usim= [usim; w];
        tsim= [tsim; t];
        msim= [msim; 2];
        % % update
        x= vectorfield('ppm', x2, inc_t, x, w, d);
        t= t+inc_t;
        [sigma, zeta]= manifold('ppm', x, omega2, x2, x2_0, zeta_0);
        ind= get_stateid(x, xgrid);
    end
    if (sigma > q_ppm(1,2) || sigma < q_ppm(1,1) ||...
            zeta > q_ppm(2,2) || zeta < q_ppm(2,1))
        x= xsim(end,:);
        t= t-inc_t;
    end
    
    plot(xsim(:,1), xsim(:,2),'-k','LineWidth', LThin,...
        'Parent', axBase)
%     plot(xsim(:,1), xsim(:,2),'-k','LineWidth', LThin,...
%         'Parent', axSub)
end

% % format display
axis(axBase, 'square')
% rectangle(axBase, 'Position', [X(1,1) X(2,1) X(1,2)-X(1,1) X(2,2)-X(2,1)])
axis([-0.2 0.9 0.1 1.9])
% axis(axBase, [X(1,:) X(2,:)])
% axis(axSub, [0.1 0.25 0.7 0.9])
% axis(axBase,[-0.1 0.7 0.4 1.8])
xlabel(axBase, {'$x$[m]'}, 'interpreter','latex',...
    'FontUnits','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')

ylabel(axBase, {'$\dot{x}$[m/s]'}, 'interpreter','latex',...
    'FontUnits','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')
