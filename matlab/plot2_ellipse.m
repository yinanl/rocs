function plot2_ellipse(P, c)
% plot ellipse defined by {x: x'*P*x=c}

N= 100;
h= 2*pi/N;
a= (1:N+1)*h;

% Method 1
[U,S,~]= svd(P);

z= diag(sqrt(S\[c;c]))*[cos(a);sin(a)];
x= U\z;

plot(x(1,:),x(2,:))


% % Method 2
% M= sqrtm(inv(P/c));
% y= M*[cos(a);sin(a)];
% plot(y(1,:), y(2,:))