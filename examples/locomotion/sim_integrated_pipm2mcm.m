%%%
% Simulation of local control strategy
% PIPM-->MCM
%%%


%% display setting
global CP CM FS LThick LThin filterSize;

FS= 16; % fontsize
LThick= 2; % linewidth
LThin= 1;

CP = get(0, 'DefaultAxesColorOrder'); % b,o,y,p,g,c,r
CM = parula(5); % p->b->g->o->y
filterSize= 1;

Cgoal= 5;
Cinit= 6;


%% problem setup
% % bound of disturbance
% D= [0.05; 0.1];
D= [0; 0];
% % discretization parameters
inc_t= 0.02; % [s]
% % control space
mu= 0.02;
U= -3:mu:-1;
% % state space
X= [-0.2, 0.8; 0.2, 1.2];
eta= [0.003; 0.003];

% % keyframe states and nominal controls
x1= [0; 0.55];
x2= [0.6; 0.3];
omega1= 2.6592;
omega2= -2;

% % constants for mapping from Euclidean to manifold
Dx= 0.0002;
zeta_0 = 1e-5; 
x1_0= [Dx; sqrt(x1(2)^2+omega1^2*Dx^2)];

% % robust margins [sigma, zeta]
ds1= 0.002;
dz1= 0.002;
ds2= 0.15;
dz2= 0.9*zeta_0;
zmin= 0.6*zeta_0;


%% plot problem setting
% figure
% axBase= axes;
% hold(axBase, 'on')
% plot(x2(1),x2(2),'*k')
% plot(x1(1),x1(2),'*k')
% 
% N1=5; 
% N2=5;
% xpos= (X(1,1):0.001:X(1,2))';
% 
% plot_robustset_cells_pipm(N1,N2,xpos,ds1,dz1,omega1,x1,Dx,x1_0,zeta_0);
% plot_robustset_cells_mcm(N1,N2,xpos,ds2,dz2,omega2,x2,zeta_0,zmin,'left');
% 
% 
% % % format display
% axis(axBase, 'square')
% % rectangle(axBase, 'Position', [X(1,1) X(2,1) X(1,2)-X(1,1) X(2,2)-X(2,1)])
% axis([-0.2 0.8 0.1 1.3])
% xlabel(axBase, {'$x$[m]'}, 'interpreter','latex',...
%     'FontSize',FS,...
%     'FontName','Times', 'FontWeight','bold')
% 
% ylabel(axBase, {'$\dot{x}$[m/s]'}, 'interpreter','latex',...
%     'FontSize',FS,...
%     'FontName','Times', 'FontWeight','bold')


%% simulating cells
load('integrated/pipm2mcm/data_grids_pipm2mcm.mat');

prefix= 'integrated/pipm2mcm/data_integrate_pipm2mcm_q';
midfix= '_p';
suffix= '.mat';

rng('shuffle');
simNum= 3;

% time limit of each semi step
dt1= 2*inc_t;
dt2= 5*inc_t;

% % Q, R: the weighting factors for selecting control values
Q= 0.5;
R= 0.5;

N1=5; N2=5;
xpos= (X(1,1):0.001:X(1,2))';

figure
axBase= axes;
hold(axBase, 'on')

for i=1:N1*N2
    if (i==2)
    % % plot problem setting
    plot_robustset_cells_pipm(N1,N2,xpos,ds1,dz1,omega1,x1,Dx,x1_0,zeta_0);
    plot_robustset_cells_mcm(N1,N2,xpos,ds2,dz2,omega2,x2,zeta_0,zmin,'left');
    
    qid= [];
    for j=1:N1*N2
        
        if(j==13)
        % % load data
        filename=strcat(prefix, num2str(i), midfix, num2str(j), suffix);
        load(filename)
        
        % % plot robust sets and R_inters
        % plot initial set
        idx1= initset+ones(size(initset));
        initial= xgrid(idx1,:);
        bdinit= extract_boundary2D(initial, filterSize);
        hinit= fill(bdinit(:,1),bdinit(:,2),CP(6,:),...
            'LineStyle','none', 'Parent', axBase);
        alpha(hinit, 0.5)
        % plot goal set
        idx2= goalset+ones(size(goalset));
        goal= xgrid(idx2,:);
        bdgoal= extract_boundary2D(goal, filterSize);
        hgoal= fill(bdgoal(:,1),bdgoal(:,2),CP(5,:),...
            'LineStyle','none', 'Parent', axBase);
        alpha(hgoal, 0.5)
        
        % plot win1
        win1= xgrid(any(leastctlr1,2),:);
        if(~isempty(win1))
            bdwin1= extract_boundary2D(win1, filterSize);
            hfill1=fill(bdwin1(:,1),bdwin1(:,2),...
                CM(4,:),'LineStyle','--','LineWidth',LThin,'EdgeColor',CM(4,:));
            alpha(hfill1,0.5)
        end
        % plot win2
        win2= xgrid(any(leastctlr2,2),:);
        if(~isempty(win1))
            bdwin2= extract_boundary2D(win2, filterSize);
            hfill2=fill(bdwin2(:,1),bdwin2(:,2),...
                CM(4,:),'LineStyle','--','LineWidth',LThin,'EdgeColor',CM(4,:));
            alpha(hfill2,0.5)
        end
        % plot guard set
        idx3= guardset+ones(size(guardset));
        guard= xgrid(idx3,:);
        if(~isempty(guard))
            bdguard= extract_boundary2D(guard, filterSize);
            hguard= fill(bdguard(:,1),bdguard(:,2),CP(2,:),...
                'LineStyle','-','LineWidth',LThin,'EdgeColor',CP(2,:),...
                'Parent', axBase);
            alpha(hguard, 0.5)
        end
        
        % % prepare controller
        ctlr.m1= 'pipm';
        ctlr.q1= x1;
        ctlr.u1= omega1;
        ctlr.id_goal1= idx3;
        ctlr.id_init1= idx1;
        ctlr.least1= leastctlr1;
        ctlr.dt1= dt1;
        ctlr.ugrid1= ugrid1;
        ctlr.m2= 'mcm';
        ctlr.q2= x2;
        ctlr.u2= omega2;
        ctlr.id_goal2= idx2;
        ctlr.id_init2= idx3;
        ctlr.least2= leastctlr2;
        ctlr.dt2= dt2;
        ctlr.ugrid2= ugrid2;
        
        % simulations
        xid= 1:size(xgrid,1);
        wxid= xid(any(leastctlr1,2));
        xinit= xgrid(intersect(idx1,wxid),:);
        k= randi(size(xinit,1));
        x0= xinit(k,:);
        [tsim,xsim,usim,dsim]=simulate_trajectory(inc_t, x0, xgrid, ctlr, Q, R, D);
        plot(xsim(:,1), xsim(:,2),'-k','LineWidth', LThin);
        
        qid= union(qid, wxid);
        end  % if(j==)
    end
    
    if(size(qid,1)<size(idx1,1))
        plot(xgrid(setdiff(qid,idx1),1), xgrid(setdiff(qid,idx1),2), '.k')
        disp(strcat('q', num2str(i), ' is not reachable'))
    end
    
    % % format display
    axis(axBase, 'square')
    % rectangle(axBase, 'Position', [X(1,1) X(2,1) X(1,2)-X(1,1) X(2,2)-X(2,1)])
    axis([-0.2 0.8 0.1 1.2])
    xlabel(axBase, {'$x$[m]'}, 'interpreter','latex',...
        'FontSize',FS,...
        'FontName','Times', 'FontWeight','bold')
    
    ylabel(axBase, {'$\dot{x}$[m/s]'}, 'interpreter','latex',...
        'FontSize',FS,...
        'FontName','Times', 'FontWeight','bold')
    end % if(i==)
end

% % % backward computation of nominal trajectory
% dec_t= -0.02;
% % xsim1= xsim(end,:);
% xsim1= x2';
% tsim1= 0;
% x= xsim1;
% t=tsim1;
% while (t>-0.4)
%     % % update
%     x= vectorfield('mcm', ctlr.q2, dec_t, x, -2, [0;0]);
%     t= t+dec_t;
%     
%     xsim1= [xsim1; x];
%     tsim1= [tsim1; t];
% end
% plot(xsim1(:,1), xsim1(:,2),'-b','LineWidth', LThick)
