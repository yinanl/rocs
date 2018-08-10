function [x,y] = calc_rect_angle(xc,yc,w,h,a)
% plot an rectangle with an rotation.

r= sqrt((w/2)^2+(h/2)^2);

b_tr= atan(h/w);
b_br= atan(-h/w);
b_bl= b_tr-pi;
b_tl= b_br+pi;


% after applying angle a
tr= wrap_angle(b_tr+a);
br= wrap_angle(b_br+a);
bl= wrap_angle(b_bl+a);
tl= wrap_angle(b_tl+a);

% obtain the position of the 4 points
n= 4;
x= zeros(1,n);
y= zeros(1,n);
x(1)= xc+r*cos(tr); y(1)= yc+r*sin(tr);
x(2)= xc+r*cos(br); y(2)= yc+r*sin(br);
x(3)= xc+r*cos(bl); y(3)= yc+r*sin(bl);
x(4)= xc+r*cos(tl); y(4)= yc+r*sin(tl);


function beta= wrap_angle(alpha)
% translate angle into range [-pi,pi]
b= mod(alpha,2*pi);
if(b>pi)
    beta= b-2*pi;
else
    beta= b;
end