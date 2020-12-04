function dx = engineMG(t, x, u)
global aH H2 W W2 cx cy

dx= [cx * (aH+H2*(x(1)-W).*(W2-(x(1)-W).^2) - x(2)) + u(1);
    cy * (x(1) - u(2).*sqrt(x(2)))];
end