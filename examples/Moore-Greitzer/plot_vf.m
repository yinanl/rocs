addpath('../../matlab')

load('data_grid_U.mat')
s= 10;
X= [0.44, 0.54; 0.63*s, 0.67*s];
e= 0.003;
% G= [0.4519-e, 0.4519+e;
%     0.6513-e, 0.6513+e];
G= [0.5039-e, 0.5039+e;
    0.6605-e, 0.6605+e].*[1; s];

% % plot vector field
N= [10; 10];
w= (X(:,2)-X(:,1))./N;
u= U(17, :);
figure
plot2_vf(X, w, {@(t,x) MG2(t,x,u)}, 1)
rectangle('Position', [G(1,1), G(2,1), G(1,2)-G(1,1), G(2,2)-G(2,1)],...
            'EdgeColor','g')

axis([X(1,:) X(2,:)])