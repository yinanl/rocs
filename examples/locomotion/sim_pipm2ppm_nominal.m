%%%
% Random simulations for the nominal case
%%%


%% display setting
FS= 16; % fontsize
LThick= 1.5; % linewidth
LThin= 1;
CP = get(0, 'DefaultAxesColorOrder'); % b,o,y,p,g,c,r
CM = parula(5); % p->b->g->o->y
filterSize= 1;

Cgoal= 5;
Cinit= 6;


%% case setting
setting_pipm2ppm

% % bound of disturbance
% D= [0.15; 0.3];
D= [0.4; 0.8];

simNum= 100;


%% plot specification
figure
axBase= axes;
hold(axBase, 'on')

% % plot robust sets
xpos= (X(1,1):0.001:X(1,2))';
v1l= real(sqrt(omega1^2*(xpos-x1(1)).^2 + x1(2)^2 + omega1^2/x1(2)^2*q_pipm(1,1)));
v1u= real(sqrt(omega1^2*(xpos-x1(1)).^2 + x1(2)^2 + omega1^2/x1(2)^2*q_pipm(1,2)));
z1l= real(x1_0(2)*((q_pipm(2,1)/zeta_0-1)*Dx./(xpos-x1(1))).^(1/omega1^2));
z1u= real(x1_0(2)*((q_pipm(2,2)/zeta_0-1)*Dx./(xpos-x1(1))).^(1/omega1^2));

v2l= real(sqrt(-omega2^2*(xpos-x2(1)).^2 + x2(2)^2 - omega2^2/x2(2)^2*q_ppm(1,2)));
v2u= real(sqrt(-omega2^2*(xpos-x2(1)).^2 + x2(2)^2 - omega2^2/x2(2)^2*q_ppm(1,1)));
z2l= real(x2_0(2)*((q_ppm(2,1)/zeta_0-1)*Dx./(xpos-x2(1))).^(-1/omega2^2));
z2u= real(x2_0(2)*((q_ppm(2,2)/zeta_0-1)*Dx./(xpos-x2(1))).^(-1/omega2^2));

v1bdu= [xpos v1u];
v1bdl= [xpos v1l];
z1bdu= [xpos z1u];
z1bdl= [xpos z1l];
% % left and right boundary of the initial robust set
rb1left= z1bdu(v1l<=z1u & z1u<=v1u & z1bdu(:,1)<=x1(1),:);
rb1right= z1bdl(v1l<=z1l & z1l<=v1u & z1bdl(:,1)>=x1(1),:);
% % upper and lower boundary of the initial robust set
rb1upper= v1bdu(rb1left(end,1)<=v1bdu(:,1) & v1bdu(:,1)<=rb1right(1,1), :);
rb1lower= v1bdl(rb1left(1,1)<=v1bdl(:,1) & v1bdl(:,1)<=rb1right(end,1), :);

v2bdu= [xpos v2u];
v2bdl= [xpos v2l];
z2bdu= [xpos z2u];
z2bdl= [xpos z2l];
% % left and right boundary of the final robust set
rb2left= z2bdl(v2l<=z2l & z2l<=v2u & z2bdl(:,1)<=x2(1),:);
rb2right= z2bdu(v2l<=z2u & z2u<=v2u & z2bdu(:,1)>=x2(1),:);
% % upper and lower boundary of the final robust set
rb2upper= v2bdu(rb2left(1,1)<=v2bdu(:,1) & v2bdu(:,1)<=rb2right(end,1), :);
rb2lower= v2bdl(rb2left(end,1)<=v2bdl(:,1) & v2bdl(:,1)<=rb2right(1,1), :);

% % % plot robust tubes
% plot(xpos, v1u, '-', 'Color', CP(5,:), 'LineWidth', LThin)
% plot(xpos, v1l, '-', 'Color', CP(5,:), 'LineWidth', LThin)
% plot(xpos, v2u, '-', 'Color', CP(5,:), 'LineWidth', LThin)
% plot(xpos, v2l, '-', 'Color', CP(5,:), 'LineWidth', LThin)
% plot(xpos, z1u, '-', 'Color', CP(5,:), 'LineWidth', LThin)
% plot(xpos, z1l, '-', 'Color', CP(5,:), 'LineWidth', LThin)
% plot(xpos, z2u, '-', 'Color', CP(5,:), 'LineWidth', LThin)
% plot(xpos, z2l, '-', 'Color', CP(5,:), 'LineWidth', LThin)

% % plot the boundary of robust sets
plot(rb1left(:,1), rb1left(:,2), '-', 'Color', CP(Cinit,:), 'LineWidth', LThick)
plot(rb1right(:,1), rb1right(:,2), '-', 'Color', CP(Cinit,:), 'LineWidth', LThick)
plot(rb1upper(:,1), rb1upper(:,2), '-', 'Color', CP(Cinit,:), 'LineWidth', LThick)
plot(rb1lower(:,1), rb1lower(:,2), '-', 'Color', CP(Cinit,:), 'LineWidth', LThick)
plot(rb2left(:,1), rb2left(:,2), '-', 'Color', CP(Cgoal,:), 'LineWidth', LThick)
plot(rb2right(:,1), rb2right(:,2), '-', 'Color', CP(Cgoal,:), 'LineWidth', LThick)
plot(rb2upper(:,1), rb2upper(:,2), '-', 'Color', CP(Cgoal,:), 'LineWidth', LThick)
plot(rb2lower(:,1), rb2lower(:,2), '-', 'Color', CP(Cgoal,:), 'LineWidth', LThick)

bdinit= [rb1left; rb1upper; rb1right; flipud(rb1lower)];
fill(bdinit(:,1),bdinit(:,2),CP(Cinit,:),...
    'LineStyle','none', 'Parent', axBase);
bdgoal= [rb2left; rb2lower; rb2right; flipud(rb2upper)];
fill(bdgoal(:,1),bdgoal(:,2),CP(Cgoal,:),...
    'LineStyle','none', 'Parent', axBase);

% % plot contact locations
plot(x2(1), x2(2), 'o', 'MarkerFaceColor', CP(7,:), 'Parent', axBase)
plot(x1(1), x1(2), 'o', 'MarkerFaceColor', CP(7,:), 'Parent', axBase)


%% simulations: nominal case
t1= 0.3342;
t2= 0.2891;

outliers= 0;

% % generate the points within the initial robust set
[xrv, xl, gtbl] = grid_uniform_generate([-0.15 0.15; 0.4 0.62], eta);
num_v= size(xrv,1);
[s1, z1]= manifold('pipm', xrv, omega1, x1, x1_0, zeta_0);
idx= find(s1>q_pipm(1,1)+d_rb1(1) & s1<q_pipm(1,2)-d_rb1(1) &...
    z1>q_pipm(2,1)+d_rb1(2) & z1<q_pipm(2,2)-d_rb1(2));
xinit= xrv(idx,:);
% for vid= 1:num_v
%     cx= xrv(vid,:)';
%     
% end

for iter= 1:simNum
    
    % % initial condition
    k= randi(size(xinit,1));
    x0= xinit(k,:);
    x= x0;
    w= omega1;
    t= 0;
    
    dsim= [];
    xsim= [];
    usim= [];
    tsim= [];
    
    % % first semi-step: time duration t1
    while (t<t1)
        % % generate a random disturbance
        d= D.*([-1;-1] + 2*rand(2,1));
        % % record
        dsim= [dsim; d'];
        xsim= [xsim; x];
        usim= [usim; w];
        tsim= [tsim; t];
        % % update
        x= vectorfield('pipm', x1, inc_t, x, w, d);
        t= t+inc_t;
    end
    
    % % second semi-step
    w= omega2;
    while (t<t1+t2)
        % % generate a random disturbance
        d= D.*([-1;-1] + 2*rand(2,1));
        dsim= [dsim; d'];
        xsim= [xsim; x];
        usim= [usim; w];
        tsim= [tsim; t];
        
        % % update
        x= vectorfield('ppm', x2, inc_t, x, w, d);
        t= t+inc_t;
    end    
    plot(xsim(:,1), xsim(:,2),'-k','LineWidth', LThin,...
        'Parent', axBase)
    
    % % determine if it ends in the final robust set
    [s2, z2]= manifold('ppm', x, omega2, x2, x2_0, zeta_0);
    if (s2<=q_ppm(1,1)||s2>=q_ppm(1,2)||z2<=q_ppm(2,1)||z2>=q_pipm(2,2))
        outliers = outliers + 1;
    end
end




%% format display
axis(axBase, 'square')
% rectangle(axBase, 'Position', [X(1,1) X(2,1) X(1,2)-X(1,1) X(2,2)-X(2,1)])
axis(axBase, [-0.2 0.9 0.2 1.95])
xlabel(axBase, {'$x$[m]'}, 'interpreter','latex',...
    'FontUnits','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')

ylabel(axBase, {'$\dot{x}$[m/s]'}, 'interpreter','latex',...
    'FontUnits','points', 'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')
