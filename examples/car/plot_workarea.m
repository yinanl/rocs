
addpath('../../matlab/')
%% load specification
load('data_car_spec1.mat') % X, U, ts, G, xobs
x0= [7.6; 0.4; pi/2];


%% define color
pink = [255,182,193]/255;
gold = [1,0.84,0];
lightgold = [255,236,139]/255;
gray = [0.6,0.6,0.6];
lightblue = [176 226 255]/255;

orange= [0.8500 0.3250 0.0980];


%% plot
figure
hold on

rectangle('Position', [X(1,1), X(2,1), X(1,2)-X(1,1), X(2,2)-X(2,1)],...
    'EdgeColor','k', 'LineWidth',2)

% % avoid area
if (~isempty(xobs))
    % real
    for i= 1:size(xobs, 3)
        rectangle('Position', [xobs(1,1,i), xobs(2,1,i), ...
            xobs(1,2,i)-xobs(1,1,i), xobs(2,2,i)-xobs(2,1,i)],...
            'EdgeColor',gray, 'FaceColor',gray)
    end
end

% % goal area
% real
rectangle('Position', [G(1,1),G(2,1),G(1,2)-G(1,1),G(2,2)-G(2,1)],...
    'EdgeColor',gold,'FaceColor',gold)

% initial area
w=0.16; % car width
h=0.08; % car length
plot_rectangle_angle(x0(1),x0(2),w,h,x0(3))
% plot(x0(1),x0(2),'Marker','p','MarkerEdgeColor','r','MarkerFaceColor','r')

% text
% text(7.7,0.4,'A_I','FontSize',16)
% text(8.25,0.8,'A_{o}','FontSize',16)
% text(9.2,0.2,'A_{g}','FontSize',16)

axis equal
axis([X(1,:) X(2,:)])

xlabel('x position')
ylabel('y position')