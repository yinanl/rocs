
addpath('../../matlab')
%% load spec & controller
% data saved in .mat:
% - Tree-structrued controller: ctree, cindex, cvalue.
% - All input values: U.
% - Workarea: X.
% - Sampling time: ts.
% - Target area: G.
load('data_ipdlCobuchi.mat')


%% diplay winning set
% w= ctree(cindex(cindex(:,2)>0,1),2:5);
% % ac= [(pavings(:,1)+pavings(:,2))/2, (pavings(:,3)+pavings(:,4))/2];
% ct= [(w(:,1)+w(:,2))/2, (w(:,3)+w(:,4))/2];
% dia= [w(:,2)-w(:,1), w(:,4)-w(:,3)];

wbox= pavings(tag==1,:);
w= wbox(1:5:end,:);
ct= [(w(:,1)+w(:,2))/2, (w(:,3)+w(:,4))/2];

figure
hold on

for i= 1:size(ct,1)
    plot(ct(i,1), ct(i,2), '.', 'markersize', 12) 
end
% plot2_boxes(w)

% axis equal
axis([X(1,:) X(2,:)])