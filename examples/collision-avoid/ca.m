function dx = ca(t, x, u, d)
dx = [-u(1)+d(1)*cos(x(3))+u(2)*x(2);
    d(1)*sin(x(3))-u(2)*x(1);
    d(2)-u(2)];
end