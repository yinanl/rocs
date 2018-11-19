%%%
% Compare the size of winning set under different levels of uncertainty
%%%

clear

%% display setting
FS= 16; % fontsize
LW= 1.5; % lineweight
CP = get(0, 'DefaultAxesColorOrder');
CNum= 6;
CM= parula(CNum);
filterSize= 1;


%% read problem data from file
setting_pipms

% % load specifications
% xgrid: a grid of the state space
% ugrid: a set of control values
% goalset: target set indices
% initset: initial set indices
load('pipms/data_walk_pipms_D0_004.mat');

% % inital set
idx1= initset+ones(size(initset));
initial= xgrid(idx1,:);
% % goal set
idx2= goalset+ones(size(goalset));
goal= xgrid(idx2,:);
% % guard set
idx3= guardset+ones(size(guardset));
guard= xgrid(idx3,:);

% % permissible states: bounded by asymptotes
xasym= X(1,1):0.01:X(1,2);
yasym11= k1*(xasym-x1(1));
yasym12= -k1*(xasym-x1(1));
yasym21= k2*(xasym-x2(1));
yasym22= -k2*(xasym-x2(1));


%% plot specifications
hfig= figure;
hold on

% % plot the inital, goal and guard sets
% plot(initial(:,1), initial(:,2), '.')
bdinit= extract_boundary2D(initial,filterSize);
fill(bdinit(:,1),bdinit(:,2),CP(6,:),'LineStyle','none');
% plot(goal(:,1), goal(:,2), '.')
bdgoal= extract_boundary2D(goal,filterSize);
fill(bdgoal(:,1),bdgoal(:,2),CP(5,:),'LineStyle','none');
% plot(guard(:,1), guard(:,2), '.k')

% % plot the contact locations
plot(x2(1), x2(2), 'o', 'MarkerFaceColor', CP(7,:))
plot(x1(1), x1(2), 'o', 'MarkerFaceColor', CP(7,:))

% % plot the asymptotes
% plot(xasym,yasym11)
% plot(xasym,yasym12)
% plot(xasym,yasym21)
% plot(xasym,yasym22)


%% plot winning sets
for N=0:CNum-1
    % % identify file names
    filename= strcat('pipms/data_walk_pipms_D', num2str(N), '_004.mat');
    
    % % load control synthesis results
    % leastctlr: N x (1+M); least restrictive controller; 1st col denotes
    % winning set.
    % optctlr: N x 3; (time) optimal controller; winning set, control index,
    % optimal value.
    load(filename);
    
    % % plot winning set of semistep 1
    win1= xgrid(any(leastctlr1,2),:);
    % plot(win1(:,1), win1(:,2), '.') % plot the discrete winning set 1
    bdwin1= extract_boundary2D(win1,filterSize);
%     hfill1= fill(bdwin1(:,1),bdwin1(:,2),CM(CNum-N,:),...
%         'LineStyle','--','LineWidth',LW,'EdgeColor',CM(CNum-N,:));
%     alpha(hfill1,0.5)
    
    % % plot winning set of semistep 2
    prewin2= xgrid(any(leastctlr2,2),:);
    [s, z]= manifold('pipm', prewin2, omega1, x1, x1_0, zeta_0); % trim winning set 2
    win2= prewin2(s<=rbset1(1,2),:);
    % plot(win2(:,1), win2(:,2), '.') % plot the discrete winning set 1
    bdwin2= extract_boundary2D(win2,filterSize);
%     hfill2= fill(bdwin2(:,1),bdwin2(:,2),CM(CNum-N,:),...
%         'LineStyle','--','LineWidth',LW,'EdgeColor',CM(CNum-N,:));
%     alpha(hfill2,0.5)
    
    % % merge two winning sets
    N1= size(bdwin1,1)/2;
    win1lower= bdwin1(1:N1,:);
    win1upper= bdwin1((N1+1):end, :);
    N2= size(bdwin2,1)/2;
    win2lower= bdwin2(1:N2,:);
    win2upper= bdwin2((N2+1):end, :);
    [~, ial1, ibl1]= intersect(win1lower(:,1),win2lower(:,1));
    [~, ial2]= setdiff(win1lower(:,1),win2lower(:,1));
    [~, ibl2]= setdiff(win2lower(:,1),win1lower(:,1));
    winlower= [win1lower(ial2,:); min(win1lower(ial1,:),win2lower(ibl1,:));...
        win2lower(ibl2,:)];
    [~, iau1, ibu1]= intersect(win1upper(:,1),win2upper(:,1));
    [~, iau2]= setdiff(win1upper(:,1),win2upper(:,1));
    [~, ibu2]= setdiff(win2upper(:,1),win1upper(:,1));
    winupper= [win1upper(iau2,:); max(win1upper(iau1,:),win2upper(ibu1,:));...
        win2upper(ibu2,:)];
    bdwin= [winlower;flip(winupper)];
    hfill2= fill(bdwin(:,1),bdwin(:,2),CM(CNum-N,:),...
        'LineStyle','--','LineWidth',LW,'EdgeColor',CM(CNum-N,:));
    alpha(hfill2,0.5)
end


% % Format display
axis square
% axis([X(1,1) X(1,2) X(2,1) X(2,2)])
rectangle('Position', [X(1,1) X(2,1) X(1,2)-X(1,1) X(2,2)-X(2,1)])
axis([-0.2 0.8 0.1 1.2])
xlabel({'$x$[m]'}, 'interpreter','latex',...
    'FontUnits','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')

ylabel({'$\dot{x}$[m/s]'}, 'interpreter','latex',...
    'FontUnits','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')

cbar= colorbar;
cbar.Direction= 'reverse';
cbar.TicksMode= 'manual';
cbar.Ticks= (0:CNum-1)/(CNum-1);
cbar.TickLabelsMode= 'manual';
cbar.TickLabels= {'D_r=[0.1; 0.2]','D_r=[0.09; 0.16]','D_r=[0.07; 0.12]',...
    'D_r=[0.05; 0.1]','D_r=[0.01; 0.02]','D_r=[0; 0]'};

