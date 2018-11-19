function y= vectorfield(mode, xc, inc_t, x, u, d)
% The vectorfield for different modes.
% INPUTS
% ---------
% mode: the walking mode (pipm, ppm, slm, mcm, hm, sm).
% xc: the contact location
% inc_t: incremental time horizon
% x: initial state (2D)
% u: control input (the angular velocity omega)
% d: disturbance (2D)

switch mode
    case 'pipm'
        A= [0 1; u^2 0];
        b= [0; -u^2*xc(1)];
        F= expm(A*inc_t);
        G= A\(F-eye(2));
        y= (F*x'+G*b + G*d)';
        
    case 'ppm'
        A= [0 1; -u^2 0];
        b= [0; u^2*xc(1)];
        F= expm(A*inc_t);
        G= A\(F-eye(2));
        y= (F*x'+G*b + G*d)';
        
    case 'mcm'
        y= [x(1)+(x(2)+d(1))*inc_t+0.5*inc_t^2*(u+d(2)), ...
            x(2)+(d(2)+u)*inc_t];
        
    case 'hm'
        y= [x(1)+(x(2)+d(1))*inc_t+0.5*inc_t^2*d(2), ...
            x(2)+d(2)*inc_t];
        
    otherwise
        y= x;
end

