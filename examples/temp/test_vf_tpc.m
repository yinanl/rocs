% parameters
r1= 0.002;
r2= 0.1;
c= 16;

ts= 10;

% test intervals
x1= interval([18, 20]);
x2= interval(20, 22);


%% simulate system trajectory
vf1= @(t, x)([-r1*(x(1)-c); 0]);      % off
vf2= @(t, x)([-r1*(x(1)-x(2)); r2]);  % heating
vf3= @(t, x)([-r1*(x(1)-x(2)); -r2]); % cooling
vf4= @(t, x)([-r1*(x(1)-x(2)); 0]);   % on

[t1, ysim1]= ode45(@(t, x) vf1(t,x), [0, ts], [19, 21]);
[t2, ysim2]= ode45(@(t, x) vf2(t,x), [0, ts], [19, 21]);
[t3, ysim3]= ode45(@(t, x) vf3(t,x), [0, ts], [19, 21]);
[t4, ysim4]= ode45(@(t, x) vf4(t,x), [0, ts], [19, 21]);

% interval propogation
y1= [exp(-r1*ts)*x1+(1-exp(-r1*ts))*c, x2];
y2= [exp(-r1*ts)*x1+(1-exp(-r1*ts))*(x2-r2/r1)+r2*ts, x2+r2*ts];
y3= [exp(-r1*ts)*x1+(1-exp(-r1*ts))*(x2+r2/r1)-r2*ts, x2-r2*ts];
y4= [exp(-r1*ts)*x1+(1-exp(-r1*ts))*x2, x2];


%% plot and compare simulations and intervals
addpath('../../matlab')
cb= [0    0.4470    0.7410];
co= [0.8500    0.3250    0.0980];
cy= [0.9290    0.6940    0.1250];
cp= [0.4940    0.1840    0.5560];
cg= [0.4660    0.6740    0.1880];

figure
hold on
plot2_boxes([x1, x2], cb, cb, 0.7)
plot2_boxes(y1,co,co,0.5)
plot2_boxes(y2,cy,cy,0.5)
plot2_boxes(y3,cp,cp,0.5)
plot2_boxes(y4,cg,cg,0.5)

plot(ysim1(:,1), ysim1(:,2))
plot(ysim2(:,1), ysim2(:,2))
plot(ysim3(:,1), ysim3(:,2))
plot(ysim4(:,1), ysim4(:,2))