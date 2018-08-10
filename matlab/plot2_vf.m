function plot2_vf(X, h, F, scale, varargin)
% PLOT2_VF(X,h,F)
% display vector field of 2-d ordinary differential equations
% 
% @params
%       X       an interval of n-dim state space
%       h       meshgrid size (2 x 1 or 1 x 2)
%       F       the functions determining the vector field
%       scale   quiver scale
%       varargin parameters of the vector field


if(size(X,1)~=2 && size(X,2)~=2)
    error('X is not a 2-d interval.')
end

if(~iscell(F))
    error('Function handles should be in a cell.')
end

nf= numel(F);

[mx,my]= meshgrid(X(1,1):h(1):X(1,2), X(2,1):h(2):X(2,2));
u=zeros(size(mx,1), size(mx,2), size(F,1));
v=zeros(size(my,1), size(my,2), size(F,1));

% figure
% hold on
nin= nargin;
if(nin > 4)
    for i=1:size(mx,1)
        for j=1:size(my,2)
%             plot(mx(i,j), my(i,j),'.')

            x=[mx(i,j); my(i,j)];
            for k= 1:nf
                f= F{k};
                dx= f(0,x,varargin{1});
                u(i,j,k)=dx(1);
                v(i,j,k)=dx(2);
            end

        end
    end
else
    for i=1:size(mx,1)
        for j=1:size(my,2)
%             plot(mx(i,j), my(i,j),'.')

            x=[mx(i,j); my(i,j)];
            for k= 1:nf
                f= F{k};
                dx= f(0,x);
                u(i,j,k)=dx(1);
                v(i,j,k)=dx(2);
            end

        end
    end
end

for k= 1:nf
    quiver(mx,my,u(:,:,k),v(:,:,k),scale);
end