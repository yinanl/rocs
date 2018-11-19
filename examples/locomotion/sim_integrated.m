%%%
% Integrated simulation of the following sequence of walking steps:
% PIPM-->PIPM-->PPM-->PIPM-->MCM-->HM-->MCM-->PIPM
%
%%%


%% display setting
global CP CM FS LThick LThin filterSize;

FS= 16; % fontsize
LThick= 2; % linewidth
LThin= 1;

CP = get(0, 'DefaultAxesColorOrder'); % b,o,y,p,g,c,r
CM = parula(5); % p->b->g->o->y
filterSize= 1;


%% Problem data
% % bound of disturbance
D= [0.05; 0.1];
% % discretization parameters
% inc_t= 0.02; % [s]
% % control space
mu= 0.02;
U= 2:mu:4;
Umcml= -3:mu:-1;
Umcmr= 1:mu:3;

% % keyframe states
% % PIPM
q_pipm= [0; 0.55];
w_pipm= 2.6592;
% % PPM
q_ppm= [0; 1.7];
w_ppm= 3.2300;
% % MCM
q_mcm= [0; 0.3];
w_mcml= -2;
w_mcmr= 2;

% % the distance between two keyframe states
% % PIPM-->PIPM
dx_pipms= 0.5;
% % PIPM-->PPM/PPM-->PIPM
dx_pipm_ppm= 0.7;
% % PIPM-->MCM/MCM-->PIPM
dx_pipm_mcm= 0.6;
% % MCM-->HM-->MCM
dx_mcm_hm= 0.8;

% % constants for mapping from Euclidean to manifold
Dx= 0.0002;
x0_pipm= [Dx; sqrt(q_pipm(2)^2+w_pipm^2*Dx^2)];
x0_ppm= [Dx; sqrt(q_ppm(2)^2-w_ppm^2*Dx^2)];

zeta_0= 1e-5;
zmin_mcml= 0.6*zeta_0;
zmin_mcmr= -zeta_0;

% % robust margins [sigma, zeta]
ds_pipm= 0.002; dz_pipm= 0.002;
ds_ppm= 0.04; dz_ppm= 0.003;
ds_mcm= 0.15; dz_mcml= 0.9*zeta_0; dz_mcmr= 0.9*zeta_0;
ds_hm= 0.2; % 

% % zeta boundaries for R_inter
zeta1_Rinter= [0.5, 10];
zeta2_Rinter= [-1000000, -0.1];


%% define a sequence
% % PIPM-->PIPM-->PPM-->PIPM
% % the state space for the walking sequence
% X= [-0.2 2.1; 0.2 1.9]; 
X= [-0.2 3.8; 0.2 1.9];
eta= [0.003; 0.003];
xpts= (X(1,1):0.001:X(1,2))';

figure
ax1= axes;
hold(ax1, 'on')

N1=5; N2=5;

% used for prioritizing cells
i_s= [3,2,4,1,5];
i_z= [3,2,4,1,5];
% i_s= randperm(N1);
% i_z= randperm(N2);

% % plot problem setting
x1= [0; q_pipm(2)];
x2= [dx_pipms; q_pipm(2)];
x3= [x2(1)+dx_pipm_ppm; q_ppm(2)];
x4= [x3(1)+dx_pipm_ppm; q_pipm(2)];
x5= [x4(1)+dx_pipm_mcm; q_mcm(2)];
x6= [x5(1)+dx_pipm_mcm; q_pipm(2)];
x7= [x6(1)+dx_pipms; q_pipm(2)];
plot_robustset_cells_pipm(N1,N2,xpts,ds_pipm,dz_pipm,w_pipm,...
    x1,Dx,x0_pipm,zeta_0);
plot_robustset_cells_pipm(N1,N2,xpts,ds_pipm,dz_pipm,w_pipm,...
    x2,Dx,x0_pipm,zeta_0);
plot_robustset_cells_ppm(N1,N2,xpts,ds_ppm,dz_ppm,w_ppm,...
    x3,Dx,x0_ppm,zeta_0);
plot_robustset_cells_pipm(N1,N2,xpts,ds_pipm,dz_pipm,w_pipm,...
    x4,Dx,x0_pipm,zeta_0);
plot_robustset_cells_mcm(N1,N2,xpts,ds_mcm,dz_mcml,w_mcml,...
    x5,zeta_0,zmin_mcml,'left');
plot_robustset_cells_mcm(N1,N2,xpts,ds_mcm,dz_mcml,w_mcmr,x5,...
    zeta_0,zmin_mcmr,'right');
plot_robustset_cells_pipm(N1,N2,xpts,ds_pipm,dz_pipm,w_pipm,...
    x6,Dx,x0_pipm,zeta_0);
plot_robustset_cells_pipm(N1,N2,xpts,ds_pipm,dz_pipm,w_pipm,...
    x7,Dx,x0_pipm,zeta_0);

% % format display
% axis(axBase, 'square')
axis([X(1,:) 0.1 2.0])
xlabel(ax1, {'$x$[m]'}, 'interpreter','latex',...
    'FontUnits','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')

ylabel(ax1, {'$\dot{x}$[m/s]'}, 'interpreter','latex',...
    'FontUnits','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')


%% simulating walking trajectories
rng('shuffle');

% % simulation time step 1ms
% inc_t= 0.001;
inc_t= 0.02;
% % time limit of each semi step
dt1= 2*inc_t;
dt2= 5*inc_t;
% % Q, R: the weighting factors for selecting control values
Q= 0.5;
R= 0.5;

% % strings
str1= 'integrated/';
str2= '/data_integrate_';
str3= '/data_grids_';
midfix1= '_q';
midfix2= '_p';
suffix= '.mat';

% % define the mode sequence
m_seq= {'pipms', 'pipm2ppm', 'ppm2pipm',...
    'pipm2mcm', 'mcm2pipm', 'pipms'}; % the mode sequence
q1= q_pipm;
q2= [dx_pipms+q1(1);q_pipm(2)];
q3= [dx_pipm_ppm+q2(1); q_ppm(2)];
q4= [dx_pipm_ppm+q3(1); q_pipm(2)];
q5= [dx_pipm_mcm+q4(1); q_mcm(2)];
q6= [dx_pipm_mcm+q5(1); q_pipm(2)];
q7= [dx_pipms+q6(1);q_pipm(2)];
q_seq= [q1, q2, q3, q4, q5, q6, q7]; % the keyframe state sequence

% % loop initial state
x_inits= [-0.06 0.47;...
    0.07 0.38;...
    -0.03 0.63;...
    -0.06 0.57;...
    0 0.52;...
    0.09 0.57];
result1.tsim=[];
result1.xsim=[];
result1.usim=[];
result1.tswitch=[];
result= repmat(result1,size(x_inits,1),1);
for k=1:size(x_inits, 1)
    tic
    % % initilization
    x0= x_inits(k,:);
    x= x0;

    t= 0;
    q0= q_seq(:,1);
    dsim= []; xsim= []; usim= []; tsim= []; msim= {};
    tswitch= [];
    % % loop the modes
    for m=1:size(m_seq,2)
        
        % % load the local grid for m_th walking step
        % xgrid, ugrid1, 2
        prefix= strcat(str1, m_seq{m}, str2, m_seq{m});
        filegrid= strcat(str1, m_seq{m}, str3, m_seq{m}, suffix);
        load(filegrid);
        id_xall= 1:size(xgrid,1);
        
        % % parameters setup
        switch m_seq{m}
            case 'pipms'
                m1= 'pipm';
                u1= w_pipm;
                x1_0= x0_pipm;
                m2= 'pipm';
                u2= w_pipm;
                x2_0= x0_pipm;
                
                ds1= ds_pipm;
                dz1= dz_pipm;
                zmin= 0;
                i_s= [randperm(N1-2)+1,1,N1];
                i_z= [randperm(N2-2)+1,1,N2];
                
            case 'pipm2ppm'
                m1= 'pipm';
                u1= w_pipm;
                x1_0= x0_pipm;
                m2= 'ppm';
                u2= w_ppm;
                x2_0= x0_ppm;
                
                ds1= ds_pipm;
                dz1= dz_pipm;
                zmin= 0;
                i_s= [randperm(N1-2)+1,1,N1];
                i_z= [randperm(N2-2)+1,1,N2];
                
            case 'ppm2pipm'
                m1= 'ppm';
                u1= w_ppm;
                x1_0= x0_ppm;
                m2= 'pipm';
                u2= w_pipm;
                x2_0= x0_pipm;
                
                ds1= ds_ppm;
                dz1= dz_ppm;
                zmin= 0;
                i_s= [randperm(N1-2)+1,1,N1];
                i_z= [randperm(N2-2)+1,1,N2];
                
            case 'pipm2mcm'
                m1= 'pipm';
                u1= w_pipm;
                x1_0= x0_pipm;
                m2= 'mcm';
                u2= w_mcml;
                
                ds1= ds_pipm;
                dz1= dz_pipm;
                zmin= zmin_mcml;
                i_s= [randperm(N1-2),4,5];
                i_z= [randperm(N2-2),4,5];
%                 i_s= [1,2,3,4,5];
%                 i_z= [2,1,3,4,5];
                
            case 'mcm2pipm'
                m1= 'mcm';
                u1= w_mcmr;
                m2= 'pipm';
                u2= w_pipm;
                x2_0= x0_pipm;
                
                ds1= ds_mcm;
                dz1= dz_mcmr;
                zmin= zmin_mcmr;
                i_s= [randperm(N1-2)+1,1,N1];
                i_z= [randperm(N2-2)+1,1,N2];
        end
        
        % % calculate the local initial and final keyframe states
        q1_rel= [q_seq(1,m)-q0(1); q_seq(2,m)];
        q2_rel= [q_seq(1,m+1)-q0(1); q_seq(2,m+1)];
        
        % % get the cell number of the initial state
        x_rel= [x(1)-q0(1), x(2)];
        [sigma, zeta]= manifold(m1, x_rel, u1, q1_rel, x1_0, zeta_0);
        is_q= floor((floor(sigma/ds1)+N1)/2)+1;
%         iz_q= floor((floor(zeta/dz1)+N2)/2)+1;
        iz_q= floor((floor((zeta-zmin)/dz1)+N2)/2)+1;
        id_q= (is_q-1)*N1+iz_q;
        
        % % N1*N2 sets of controller data for one cell in initial robust set
        % % choose the one controls the state closest to the nominal state
        % % data:
        % % goalset, initset, guardset,
        % % leastctlr1, leastctlr2, optctlr1, optctlr2
        ctlr.id_goalset=[];
        ctlr.id_initset=[];
        ctlr.id_guardset=[];
        ctlr.least1=[];
        ctlr.least2=[];
        ctlr.opt1=[];
        ctlr.opt2=[];
        b= 0;
        for j=1:size(i_z,2)
            for i=1:size(i_s,2)
                id_p= (i_s(i)-1)*N1+i_z(j);
                filectlr=strcat(prefix,midfix1,num2str(id_q),...
                    midfix2,num2str(id_p),suffix);
                load(filectlr)
%                 id_initset= initset+ones(size(initset));
%                 id_win1= intersect(id_initset,id_xall(any(leastctlr1,2)));
                id_win1= id_xall(any(leastctlr1,2));
                id_xgrid= get_stateid(x_rel, xgrid);
                if(ismember(id_xgrid, id_win1))
                    ctlr.id_goalset= goalset+ones(size(goalset));
                    ctlr.id_initset= initset+ones(size(initset));
                    ctlr.id_guardset= guardset+ones(size(guardset));
                    ctlr.least1= leastctlr1;
                    ctlr.least2= leastctlr2;
                    ctlr.opt1= optctlr1;
                    ctlr.opt2= optctlr2;
                    b= 1;
                    break;
                end
            end
            if(b)
                break;
            end
        end
        
        % % plot the winsets
        % plot initial set
        initial= xgrid(ctlr.id_initset,:)+repmat([q0(1) 0], size(ctlr.id_initset,1),1);
        if(~isempty(initial))
            bdinit= extract_boundary2D(initial, filterSize);
            hinit= fill(ax1, bdinit(:,1),bdinit(:,2),CP(6,:),'LineStyle','none');
            alpha(hinit, 0.5)
        end
        % plot goal set
        goal= xgrid(ctlr.id_goalset,:)+repmat([q0(1) 0], size(ctlr.id_goalset,1),1);
        if(~isempty(goal))
            bdgoal= extract_boundary2D(goal, filterSize);
            hgoal= fill(ax1, bdgoal(:,1),bdgoal(:,2),CP(6,:),'LineStyle','none');
            alpha(hgoal, 0.5)
        end
        
        % plot win1
        win1= xgrid(any(ctlr.least1,2),:);
%         if(~isempty(win1))
%             win1= win1+repmat([q0(1) 0], size(win1,1),1);
%             bdwin1= extract_boundary2D(win1, filterSize);
%             hfill1=fill(ax1, bdwin1(:,1),bdwin1(:,2),...
%                 CM(4,:),'LineStyle','--','LineWidth',LThin,'EdgeColor',CM(4,:));
%             alpha(hfill1,0.5)
%         end
        % plot win2
        win2= xgrid(any(ctlr.least2,2),:);
%         if(~isempty(win2))
%             win2= win2+repmat([q0(1) 0], size(win2,1),1);
%             bdwin2= extract_boundary2D(win2, filterSize);
%             hfill2=fill(ax1, bdwin2(:,1),bdwin2(:,2),...
%                 CM(4,:),'LineStyle','--','LineWidth',LThin,'EdgeColor',CM(4,:));
%             alpha(hfill2,0.5)
%         end
        % plot guard set
        guard= xgrid(ctlr.id_guardset,:)+repmat([q0(1) 0], size(ctlr.id_guardset,1),1);
%         if(~isempty(guard))
%             bdguard= extract_boundary2D(guard, filterSize);
%             hguard= fill(ax1, bdguard(:,1),bdguard(:,2),CP(2,:),...
%                 'LineStyle','-','LineWidth',LThin,'EdgeColor',CP(2,:));
%             alpha(hguard, 0.5)
%         end
        
        % % execute control simulation
        x= x_rel; w= u1;
        % % first semi-step: m1
        id_xgrid= get_stateid(x_rel, xgrid);
        while (~any(ismember(id_xgrid, ctlr.id_guardset)))
            % % compute control input
            uall= find(ctlr.least1(id_xgrid,2:end));
            [~,uid]= min(Q*(ugrid1(uall)-u1).^2 + R*(ugrid1(uall)-w).^2);
            w= ugrid1(uall(uid));
            if(isempty(w))
                if (~isempty(xsim))
                    plot(xsim(:,1), xsim(:,2),'-b','LineWidth', LThick)
                end
                keyboard
                error('Not within the winning set 1.')
            end
            
            % % generate a random disturbance
            d= D.*([-1;-1] + 2*rand(2,1));
            % % record
            dsim= [dsim; d'];
            xsim= [xsim; [x(1)+q0(1), x(2)]]; % shift in position
            usim= [usim; w];
            tsim= [tsim; t];
            msim{end+1}= m1;
            % % update
            x= vectorfield(m1, q1_rel, inc_t, x, w, d);
            t= t+inc_t;
            % % check update state
            id_xgrid= get_stateid(x, xgrid);
        end
        t1= t;
        tswitch= [tswitch;t1];
        
        % % continue mode 1 for additional dt1 time
        while (any(ismember(id_xgrid, ctlr.id_guardset)) && t<t1+dt1)
            % % generate a random disturbance
            d= D.*([-1;-1] + 2*rand(2,1));
            dsim= [dsim; d'];
            xsim= [xsim; [x(1)+q0(1), x(2)]]; % shift in position
            usim= [usim; w];
            tsim= [tsim; t];
            msim{end+1}= m1;
            % % update
            x= vectorfield(m1, q1_rel, inc_t, x, w, d);
            t= t+inc_t;
            id_xgrid= get_stateid(x, xgrid);
        end
        if (~any(ismember(id_xgrid, ctlr.id_guardset)))
            x= [xsim(end,1)-q0(1), xsim(end,2)];
            t= t-inc_t;
        end
        
        % % second semi-step
        w= u2;
        id_xgrid= get_stateid(x, xgrid);
        while (~any(ismember(id_xgrid, ctlr.id_goalset)))
            % % compute control input
            id_xgrid= get_stateid(x, xgrid);
            uall= find(ctlr.least2(id_xgrid,2:end));
            [~,uid]= min(Q*(ugrid2(uall)-u2).^2 + R*(ugrid2(uall)-w).^2);
            w= ugrid2(uall(uid));
            if(isempty(w))
                if (~isempty(xsim))
                    plot(xsim(:,1), xsim(:,2),'-b','LineWidth', LThick)
                end
                keyboard
                error('Not within the winning set 2.')
            end
            
            % % generate a random disturbance
            d= D.*([-1;-1] + 2*rand(2,1));
            dsim= [dsim; d'];
            xsim= [xsim; [x(1)+q0(1), x(2)]]; % shift in position
            usim= [usim; w];
            tsim= [tsim; t];
            msim{end+1}= m2;
            % % update
            x= vectorfield(m2, q2_rel, inc_t, x, w, d);
            t= t+inc_t;
            id_xgrid= get_stateid(x, xgrid);
        end
        t2= t;
        tswitch= [tswitch;t2];
        
        % % continue the last control config for ~0.05s
        while (any(ismember(id_xgrid, ctlr.id_goalset)) && t<t2+dt2)
            % % generate a random disturbance
            d= D.*([-1;-1] + 2*rand(2,1));
            dsim= [dsim; d'];
            xsim= [xsim; [x(1)+q0(1), x(2)]]; % shift in position
            usim= [usim; w];
            tsim= [tsim; t];
            msim{end+1}= m2;
            % % update
            x= vectorfield(m2, q2_rel, inc_t, x, w, d);
            t= t+inc_t;
            id_xgrid= get_stateid(x, xgrid);
        end
        if (~any(ismember(id_xgrid, ctlr.id_goalset)))
            %         x= [xsim(end,1)-q0(1), xsim(end,2)];
            t= t-inc_t;
        end
        
        plot(ax1, xsim(:,1), xsim(:,2),'-k','LineWidth', LThin)
        
        % % update the current keyframe state
        x= [xsim(end,1), xsim(end,2)]; % change back to absolute position
        q0= q_seq(:,m+1);
    end
    toc
    
    % % save simulation results
    result(k).tsim= tsim;
    result(k).xsim= xsim;
    result(k).usim= usim;
    result(k).tswitch= tswitch;
    
    % % plot the position and velocity-time trajectory
    figure
%     result= [tsim xsim usim];
    axsub1= subplot(2,1,1);
    plot(axsub1, tsim, xsim, 'LineWidth', LThick);
    legend(axsub1,'position','velocity')
    axsub2= subplot(2,1,2);
    plot(axsub2, tsim, usim, 'LineWidth', LThick);
    
    title(axsub1, ['x_0=[',num2str(x0(1)),', ',num2str(x0(2)),']']);
    xlabel(axsub1, {'$t$[s]'}, 'interpreter','latex',...
        'FontUnits','points', 'FontSize',FS,...
        'FontName','Times', 'FontWeight','bold')
    ylabel(axsub1, {'$x[m],\dot{x}$[m/s]'}, 'interpreter','latex',...
        'FontUnits','points', 'FontSize',FS,...
        'FontName','Times', 'FontWeight','bold')
    xlabel(axsub2, {'$t$[s]'}, 'interpreter','latex',...
        'FontUnits','points', 'FontSize',FS,...
        'FontName','Times', 'FontWeight','bold')
    ylabel(axsub2, {'$u$'}, 'interpreter','latex',...
        'FontUnits','points', 'FontSize',FS,...
        'FontName','Times', 'FontWeight','bold')
end


%% trajectory under large disturbance
