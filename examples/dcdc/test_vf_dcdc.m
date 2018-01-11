% parameters
xc= 70;
xl= 3;
rc= 0.005;
rl= 0.05;
r0= 1;
vs= 1;
I= eye(2);
h= 0.5;

a11= -rl/xl;
a12= -1/(xc*(rc+r0));
A1= [a11 0; 0 a12];
b1= [vs/xl; 0];
F1= expm(A1*h);
G1= A1\I*(expm(A1*h)-I)*b1;

a21= (-1/xl)*(rl + r0*rc/(r0+rc));
a22= (-1/xl)*(r0/(r0+rc));
a23= (1/xc)*(r0/(r0+rc));
a24= (-1/xc)*(1/(r0+rc));
A2= [a21 a22; a23 a24];
b2= b1;
F2= expm(A2*h);
G2= A2\I*(expm(A2*h)-I)*b2;


% test intervals
% x1= interval([1.2, 1.25]);
% x2= interval(1.11, 1.12);
% x1= interval([1.4, 1.5]);
% x2= interval(1.15, 1.156);
x1= interval([1.15, 1.55]);
x2= interval(1.09, 1.17);

% interval to zonotope
xc = [mid(x1), mid(x2)];
gw = [width(x1), width(x2)]/2;

nbox= size(x1,1); % number of boxes
n= 2; % system dimension

yt1 = interval(zeros(nbox,n),[],0); % dim=2
yt2 = interval(zeros(nbox,n),[],0);
for i= 1:nbox
    xzono1= [xc(i,:)', diag(gw(i,:)')];
    yc1= F1*xzono1;
    % zonotope to interval
    yintr1= sum(abs(yc1(:, 2:end)), 2);
    yint1= [yc1(:,1)+G1-yintr1, yc1(:,1)+G1+yintr1];
    
    yt1(i, :)= (interval(yint1,[],2))';
    
    xzono2= [xc(i,:)', diag(gw(i,:)')];
    yc2= F2*xzono2;
    % zonotope to interval
    yintr2= sum(abs(yc2(:, 2:end)), 2);
    yint2= [yc2(:,1)+G2-yintr2, yc2(:,1)+G2+yintr2];
    
    yt2(i, :)= (interval(yint2,[],2))';
end

addpath('../../matlab')
figure
hold on

plot2_boxes([x1, x2], 'b', 'b', 0.7)
plot2_boxes(yt1,'r','r',0.5)
plot2_boxes(yt2,'y','y',0.5)