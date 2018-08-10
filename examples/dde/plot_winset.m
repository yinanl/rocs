
addpath('../../matlab/')
%% load specification
load('data_ddeInv.mat')


%% diplay control synthesis results

% % plot the winset and/or its subset
win= pavings(tag==1,:);
uwin= ctlr(tag==1,:);

% % additional condition
cbox= win;
% cbox= win(win(:,5)>0,:);
% uid= uwin(win(:,5)>0,:);
% uval= zeros(size(U,1),size(U,2),size(cbox,1));

cb= cbox(1:500:end,:);
ccenter= [(cb(:,1)+cb(:,2))/2, (cb(:,3)+cb(:,4))/2];
plot(ccenter(:,1), ccenter(:,2), '.', 'MarkerSize', 6)
hold on

% % % plot the center of goal area
% wbox= pavings(tag==1,:);
% wcenter= [(wb(:,1)+wb(:,2))/2, (wb(:,3)+wb(:,4))/2, (wb(:,5)+wb(:,6))/2];
% plot3(wcenter(:,1), wcenter(:,2), wcenter(:,3), '.', 'markersize', 6)

axis equal
axis([X(1,:) X(2,:)])
% view(90,0)


%% check points
% x= [4, 2, 0];
% xid= find(x(1)>=win(:,1) & x(1)<=win(:,2) & ...
%         x(2)>=win(:,3) & x(2)<=win(:,4) & ...
%         x(3)>=win(:,5) & x(3)<=win(:,6));
% uid= find(uwin(xid(1),:));
% uall= U(uid,:);