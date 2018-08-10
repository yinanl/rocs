function h= plot2_ellipse(P, c, varargin)
% PLOT2_ELLIPSE(P,c,VARARGIN)
% ------------------------- 
% plot ellipse defined by {x: x'*P*x=c}
%
% Inputs:
%
if(isempty(varargin) || nargin~= 5)
    color= 'k';
    linestyle= '-';
    linewidth= 1.5;
else
    color= varargin{1};
    linestyle= varargin{2};
    linewidth= varargin{3};
end

N= 100;
h= 2*pi/N;
a= (1:N+1)*h;

% Method 1
[U,S,~]= svd(P);

z= diag(sqrt(S\[c;c]))*[cos(a);sin(a)];
x= U\z;

h= plot(x(1,:),x(2,:),...
    'Color',color,'LineStyle',linestyle, 'LineWidth', linewidth);


% % Method 2
% M= sqrtm(inv(P/c));
% y= M*[cos(a);sin(a)];
% plot(y(1,:), y(2,:))