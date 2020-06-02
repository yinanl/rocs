
addpath('../../matlab/')
clear all
clc


%% load data and problem settings
load('data_carParking.mat')


%% plot setting
blue= [0 0.4470 0.7410];

% % state space
X= [0 8; 0 4; -pi/2.5 pi/2.5];

% % front and rear cars
olb= [0 4.3; 0 0];
oub= [2 6.3; 1 1];

% % % car size
% L=2;
% H=1;
% d=2;
% D=0.5;

% % 4 corners of the car
dx1= [L/2; H/2];
dx2= [L/2; -H/2];
dx3= [-L/2; -H/2];
dx4= [-L/2; H/2];
g= [H/2; L/2];

% % front and rear car corners
o1= [0;H];
o2= [L; H];
o3= [2*L+d; H];
o4= [3*L+d; H];

% % discretize theta
% N= 5;
theta= X(3,1):0.05:X(3,2);
% z= repmat(theta, N, 1);


%% plot analytical constraints
x11= zeros(1, numel(theta));
y11= x11;
x21= x11; y21= x11;
x31= x11; y31= x11;
x41= x11; y41= x11;

x12= x11; y12= x11;
x22= x11; y22= x11;
x32= x11; y32= x11;
x42= x11; y42= x11;

x13= x11; y13= x11;
x23= x11; y23= x11;
x33= x11; y33= x11;
x43= x11; y43= x11;

A1= [1 0 0; 0 1 0];
A2= [1 0 0; 0 1 0; -1 0 0];
A3= [0 1 0;1 0 0;-1 0 0];
P= [];
for i=1:numel(theta) %\theta\in[-72, 72]deg
    z= theta(i);
    
    % % car body area
    Q= [sin(z) -cos(z); cos(z) sin(z)];
    A= [[Q,[0;0]]; [-Q,[0;0]]];
    b1= [Q*o1+g; -Q*o1+g];
    b2= [Q*o2+g; -Q*o2+g];
    b3= [Q*o3+g; -Q*o3+g];
    b4= [Q*o4+g; -Q*o4+g];
    P= [P; ...
        Polyhedron('A',A, 'b', b1, 'Ae', [0 0 1], 'be', z);...
        Polyhedron('A',A, 'b', b2, 'Ae', [0 0 1], 'be', z);...
        Polyhedron('A',A, 'b', b3, 'Ae', [0 0 1], 'be', z);...
        Polyhedron('A',A, 'b', b4, 'Ae', [0 0 1], 'be', z)];
%     P= [P; ...
%         Polyhedron('A',A, 'b', b3, 'Ae', [0 0 1], 'be', z);...
%         Polyhedron('A',A, 'b', b4, 'Ae', [0 0 1], 'be', z)];
    
    
    R= [cos(z) -sin(z); sin(z) cos(z)];
    x1= R*dx1;
    x2= R*dx2;
    x3= R*dx3;
    x4= R*dx4;
    % % real car collision area
    br1= [[L;H]-R*dx1; x1(1)];
    br2= [[L;H]-R*dx2; x2(1)];
    br3= [[L;H]-R*dx3; x3(1)];
    br4= [[L;H]-R*dx4; x4(1)];
    P= [P;...
        Polyhedron('A',A2, 'b', br1, 'Ae', [0 0 1], 'be', z);...
        Polyhedron('A',A2, 'b', br2, 'Ae', [0 0 1], 'be', z);...
        Polyhedron('A',A2, 'b', br3, 'Ae', [0 0 1], 'be', z);...
        Polyhedron('A',A2, 'b', br4, 'Ae', [0 0 1], 'be', z)];
    
    % % front car collision area
    bf1= [[3*L+d;H]-x1;-(2*L+d)+x1(1)];
    bf2= [[3*L+d;H]-x2;-(2*L+d)+x2(1)];
    bf3= [[3*L+d;H]-x3;-(2*L+d)+x3(1)];
    bf4= [[3*L+d;H]-x4;-(2*L+d)+x4(1)];
    P= [P; Polyhedron('A',A2, 'b', bf1, 'Ae', [0 0 1], 'be', z);...
        Polyhedron('A',A2, 'b', bf2, 'Ae', [0 0 1], 'be', z);...
        Polyhedron('A',A2, 'b', bf3, 'Ae', [0 0 1], 'be', z);...
        Polyhedron('A',A2, 'b', bf4, 'Ae', [0 0 1], 'be', z)];
    
    % % collision with curb
    bc1= [-D-x1(2); X(1,2); X(1,1)];
    bc2= [-D-x2(2); X(1,2); X(1,1)];
    bc3= [-D-x3(2); X(1,2); X(1,1)];
    bc4= [-D-x4(2); X(1,2); X(1,1)];
    P= [P;...
        Polyhedron('A',A3, 'b', bc1, 'Ae', [0 0 1], 'be', z);...
        Polyhedron('A',A3, 'b', bc2, 'Ae', [0 0 1], 'be', z);...
        Polyhedron('A',A3, 'b', bc3, 'Ae', [0 0 1], 'be', z);...
        Polyhedron('A',A3, 'b', bc4, 'Ae', [0 0 1], 'be', z)];
end
P.plot('LineStyle', 'none', 'Color', [0.7 0.7 0.7])
% P.plot('Color', [0.7 0.7 0.7])
hold on

% % % rear and front car areas
% prear= [0 0;L 0;L H;0 H];
% pfront= [2*L+d 0; 3*L+d 0; 3*L+d H; 2*L+d H];
% xyrear= polygon(prear,0.25);
% fill3(xyfront(:,1),xyfront(:,2),repmat(theta(end), size(xyfront,1)), 'k')
% fill3(xyrear(:,1),xyrear(:,2),repmat(theta(end), size(xyrear,1)), 'k')


%% plot constraints from data

% % plot the center of constraints
cbox= pavings(tag==-1,:);
cb= cbox(1:100:end,:);
ccenter= [(cb(:,1)+cb(:,2))/2, (cb(:,3)+cb(:,4))/2, (cb(:,5)+cb(:,6))/2];
plot3(ccenter(:,1), ccenter(:,2), ccenter(:,3), '.', 'markersize', 6)
hold on

% % plot the center of goal area
wbox= pavings(tag==1,:);
wb= wbox(1:100:end,:);
wcenter= [(wb(:,1)+wb(:,2))/2, (wb(:,3)+wb(:,4))/2, (wb(:,5)+wb(:,6))/2];
plot3(wcenter(:,1), wcenter(:,2), wcenter(:,3), '.', 'markersize', 6)


%% Properties of the graph
FS= 16; % fontsize
LW= 1.5; % lineweight
xlabel({'$x$ position'}, 'interpreter','latex',...
    'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')
ylabel({'$y$ position'}, 'interpreter','latex',...
    'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')
zlabel({'$\theta$'}, 'interpreter','latex',...
    'FontSize',FS,...
    'FontName','Times', 'FontWeight','bold')
axis equal
axis([X(1,:) X(2,:) X(3,:)])
% view(0,0) % x-z
% view(0, 90) % x-y / view(2)
