function dx = MGaug(t, x, u)
a = 1/3.5;
B = 2.0;
H = 0.18;
W = 0.25;
lc = 8.0;
cx = 1.0/lc;
cy = 1.0/(4*lc*B*B);
aH = a+H;
H2 = H/(2.0*W*W*W);
W2 = 3*W*W;

dx= [cx * (aH+H2*(x(1)-W).*(W2-(x(1)-W).^2) - x(2)) + u(1);
    cy * (x(1) - x(3).*sqrt(x(2)));
    u(2)];
end