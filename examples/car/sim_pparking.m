clear
clc
addpath('../../matlab/')
%% load spec & controller
%%% Define ODEs or DEs for car kinematics %%%
vf= @car; % ODEs
fm= @carflow; % DEs

example= 'parking';
%%% Controller data %%%
% - ts: Sampling time.
% - U : All input values.
% - X : Workarea.
% - A(can be empty): obstacles.
% - G(can be empty): Target area.
% - pavings: Tree-structrued controller.
% - tag: indicating if a cell is inside the winning set.
% - ctlr: all valid control inputs for each cell in pavings.

%%% Load from .mat file %%%
% load('data_carParking.mat')

%%% Load from .h5 file %%%
ctlrfile= 'controller_carParking.h5';
ts= h5read(ctlrfile, '/ts');
X= h5read(ctlrfile, '/X')';
U= h5read(ctlrfile, '/U')';
A= permute(h5read(ctlrfile, '/xobs'), [3,2,1]);
G= permute(h5read(ctlrfile, '/G'), [3,2,1]);
pavings= h5read(ctlrfile, '/pavings')';
tag= h5read(ctlrfile, '/tag');
ctlr= h5read(ctlrfile, '/ctlr')';
H= h5read(ctlrfile, '/H');
L= h5read(ctlrfile, '/L');
D= h5read(ctlrfile, '/D');
d= h5read(ctlrfile, '/d');


%% numerical simulation
% x0= [2.0; 2.5; 0];
x0= [5.0; 2.5; 0];
Tsim= 10;

x= x0;
t= 0;
isim = 1;

%%% for movie playback %%%
% tmov= [];
% xmov= [];

tsim= [];
xsim= [];
usim= [];
while(t < Tsim || ...
        x(1)>G(1,2) || x(1)<G(1,1) ||...
        x(2)>G(2,2) || x(2)<G(2,1) || ...
        x(3)>G(3,2) || x(3)<G(3,1))
    isim= isim+1;
    %%% compute control input %%%
    xid= find(x(1)>=pavings(:,1) & x(1)<=pavings(:,2) & ...
        x(2)>=pavings(:,3) & x(2)<=pavings(:,4) & ...
        x(3)>=pavings(:,5) & x(3)<=pavings(:,6));
    uid= find(ctlr(xid(1),:));
    if (isempty(uid))
        error("Invalid controller.");
    else
        %%% select a random u from all valid control values
        pick=randperm(numel(uid));
        u= U(uid(pick(1)),:);
%         %%% select the minimum value
%         uall= U(uid,:);
%         [val, ind]= min(abs(uall(:,1)));
%         %[val, ind]= min(norm(uall-repmat(usim(isim-1,:),size(uall,1),1), 2));
%         u= uall(ind,:);
%         %%% select the first/last value
%         u= U(uid(end),:);
    end
    
    %%% Store simulation data %%%
    xsim= [xsim; x'];
    tsim= [tsim; t];
    usim= [usim; u];
    
    %%% Update the state of the dyanmical system %%%
%     %%% use difference equations
%     xx= fm(ts,x',u);
%     t= t+ts;
%     x= xx';
    %%% use ode
    [tt, xx]= ode45(@(t,x) vf(t,x,u), [0, ts], x);
    
    %%% for movie playback %%%
%     tmov= [tmov; repmat(t, size(tt,1)-1, 1) + tt(1:end-1)];
%     xmov= [xmov; xx(1:end-1,:)];
    
    %%% update x, t %%%
    x= xx(end,:)';
    t= t + tt(end,:);
end


%% display workspace
% define color
pink = [255,182,193]/255;
gold = [1,0.84,0];
lightgold = [255,236,139]/255;
gray = [0.6,0.6,0.6];
lightblue = [176, 226, 255]/255;
orange= [0.8500, 0.3250, 0.0980];
green= [0.4660, 0.6740, 0.1880]*0.7;

FS= 16; % fontsize
LW= 1.5; % lineweight

hf1= figure;
hold on
% % curb
curbsize= 0.5;
hext= L/2 - H/2 + curbsize;
wext= 0.5;
rectangle('Position', [X(1,1)-wext, X(2,1)-hext, X(1,2)-X(1,1)+2*wext, curbsize],...
    'EdgeColor','k', 'LineWidth',2, 'FaceColor', 'k')

% % rear and front car positions
dfr= 2*L+d; % end-to-end distance between front and rear cars

% if (~isempty(xobs))
%     for i= 1:size(xobs, 3)
%         rectangle('Position', [xobs(1,1,i), xobs(2,1,i), ...
%             xobs(1,2,i)-xobs(1,1,i), xobs(2,2,i)-xobs(2,1,i)],...
%             'EdgeColor',gray, 'FaceColor',gray)
%     end
% end
prear= [0 0;L 0;L H;0 H];
pfront= [dfr 0;dfr+L 0;dfr+L H;dfr H];
xyrear= polygon(prear,0.25);
fill(xyrear(:,1), xyrear(:,2), 'k')
xyfront= polygon(pfront,0.25);
fill(xyfront(:,1), xyfront(:,2), 'k')

% % initial position
[xrec, yrec]= calc_rect_angle(x0(1),x0(2),L,H,x0(3));
xy= polygon([xrec' yrec'],0.25);
fill(xy(:,1), xy(:,2), lightblue)

% % % goal area
% rectangle('Position', [G(1,1), G(2,1), ...
%     G(1,2)-G(1,1), G(2,2)-G(2,1)],...
%     'EdgeColor',gold,'FaceColor',gold)

axis equal
rectangle('Position', [X(1,1)-wext, X(2,1)-hext, ...
    X(1,2)-X(1,1)+2*wext, X(2,2)-X(2,1)+hext],...
    'EdgeColor','k', 'LineWidth',2)
axis([X(1,1)-wext X(1,2)+wext X(2,1)-hext X(2,2)])


xlabel({'$x$ position'}, 'interpreter','latex',...
    'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')
ylabel({'$y$ position'}, 'interpreter','latex',...
    'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')


%% animated path (movie playback)
% movfile = strcat('movie_',example,'.avi');
% v = VideoWriter(movfile);
% open(v);
% 
% frames= 1:5:size(xmov,1);
% 
% for i= 1:size(frames,2)
%     [xrec, yrec]= calc_rect_angle(xmov(frames(i),1),xmov(frames(i),2),L,H,xmov(frames(i),3));
%     xy= polygon([xrec' yrec'],0.25);
%     ppark= fill(xy(:,1), xy(:,2), lightblue);
% %     ppark= plot_rectangle_angle(xmov(frames(i),1),xmov(frames(i),2),L,H,xmov(frames(i),3));
% 
%     M(i)= getframe;
%     writeVideo(v,M(i));
%     
%     delete(ppark);
% end
% 
% for i= 1:size(xsim,1)
%     [xrec, yrec]= calc_rect_angle(xsim(i,1),xsim(i,2),L,H,xsim(i,3));
%     xy= polygon([xrec' yrec'],0.25);
%     fill(xy(:,1), xy(:,2), lightblue);
% %     plot_rectangle_angle(xsim(i,1),xsim(i,2),L,H,xsim(i,3));
% end
% plot(xsim(:,1),xsim(:,2),...
%     'LineWidth', 1, ...
%     'Marker','o','MarkerEdgeColor','r','MarkerFaceColor','r')
% 
% writeVideo(v,M(i+1));
% close(v)


%% static path

% % plot path
for i= 1:size(xsim,1)
    [xrec, yrec]= calc_rect_angle(xsim(i,1),xsim(i,2),L,H,xsim(i,3));
    xy= polygon([xrec' yrec'],0.25);
    fill(xy(:,1), xy(:,2), lightblue);
end

% plot(xsim(:,1),xsim(:,2),...
%     'LineWidth', 1, ...
%     'Marker','o','MarkerEdgeColor','r','MarkerFaceColor','r')
% plot(xsim(:,1),xsim(:,2),'Color', orange, 'LineWidth', LW)
xfr= xsim(1:end-1,:);
xto= xsim(2:end,:);
quiver(xfr(:,1),xfr(:,2), xto(:,1)-xfr(:,1), xto(:,2)-xfr(:,2), 0,...
    'LineWidth', LW, 'Color', green)


%% time-control curves
usim= [usim(2:end, :); usim(end,:)];
tq= [0:0.01:tsim(end)]';
uq= interp1(tsim,usim,tq,'previous');

hf2= figure;
plot(tq, uq, 'LineWidth', 2);

axis([0, tsim(end), -2, 2])
xlabel({'$t(s)$'}, 'interpreter','latex',...
    'FontSize',FS, 'FontName','Times', 'FontWeight','bold')
ylabel({'$v,\;\phi$'}, 'interpreter','latex',...
    'FontSize',FS, 'FontName','Times', 'FontWeight','bold')
hl= legend({'$v(t)$','$\phi(t)$'}, 'Interpreter', 'latex',...
    'FontSize',FS, 'FontName','Times', 'FontWeight','bold');


%% save control curves
%%% save workspace
% sfile= strcat(example, '_space.eps');
% print(hf1, sfile, '-depsc2')

% pfile= strcat(example, '_path.eps');
% print(hf1, pfile, '-depsc2')

% cfile= strcat(example, '_controls.eps');
% print(hf2, cfile, '-depsc2')
