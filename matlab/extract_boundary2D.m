function r = extract_boundary2D(x, windowSize)
% FUNCTION
% -----------
% Extract and smooth the boundary of a given array of 2D data points.
%
% INPUT
% -----------
% x: An (N x 2) array of data points.
%
% OUTPUT
% -----------
% r: An (M x 2) array of closed boundary points.
%

b= (1/windowSize)*ones(1, windowSize);
a= 1;

z= unique(x(:,1)); % extract boundaries
ymax= z;
ymin= z;
for i= 1:size(z,1)
    ymax(i)= max(x(x(:,1)==z(i),2));
    ymin(i)= min(x(x(:,1)==z(i),2));
end
yfmax= filter(b,a,[repmat(ymax(1), windowSize-1,1); ymax]);
yfmin= filter(b,a,[repmat(ymin(1), windowSize-1,1); ymin]);

r= [z,yfmin(windowSize:end); flipud(z),flipud(yfmax(windowSize:end))];