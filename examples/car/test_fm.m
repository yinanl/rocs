% test inclusion function computation
addpath('../../matlab/')

spec1


%% computation
x= [8 8.2; 
    0 0.15; 
    0.3927 0.5890];
xbox= interval(x,[],2)';

% u= [0.6000   -0.9000;
%     0.9000   -0.9000;
%     0.9000   -0.6000;
%     0.9000   -0.3000];
u= urv;
ybox= repmat(xbox,size(u,1),1);

for i=1:size(u,1)
    ybox(i,:)= fm(taus, xbox, u(i,:));
end


%% display
% % define color
cr= [0.6350 0.0780 0.1840];
cb= [0 0.4470 0.7410];
cy= [0.9290 0.6940 0.1250];

plot2_boxes(ybox, cr, cr, 0.5)