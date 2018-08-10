function df=vdp_disturbed(t,x,delta)
% reversed Van de Pol equations with additive disturbance

s= rng;
w= delta*(-1 + 2*rand(2,1));
% w= delta*rand(2,1);
df= [-x(2); x(1) + (x(1)^2-1)*x(2)]+w;