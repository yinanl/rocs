
addpath('../../matlab/')
%% load specification
load('data_car_spec1.mat')
load('data_car_ctree_spec1.mat')


%% diplay winning set
w= ctree(cindex(cindex(:,2)>0,1),2:7);

ct= [(w(:,1)+w(:,2))/2, (w(:,3)+w(:,4))/2, (w(:,5)+w(:,6))/2];
dia= [w(:,2)-w(:,1), w(:,4)-w(:,3), w(:,6)-w(:,5)];

figure
hold on

for i= 1:size(ct,1)
    
    plot3(ct(i,1), ct(i,2), ct(i,3), 'o')
%     plot3_box([sp(i,1); sp(i,3); sp(i,5)], ...
%         spdia(i,1), spdia(i,2), spdia(i,3), 'b', 0.3)
    
end

axis equal
axis([X(1,:) X(2,:)])