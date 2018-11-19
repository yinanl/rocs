function [sigma, zeta] = manifold(mode, x, w, xnom, x0, zeta_0)
% FUNCTION
% -----------
% Convert the Euclidean space [x,vx] to the manifold space [sigma, zeta].
%
% INPUT
% -----------
% mode: a string, the locomotion mode.
% x: An (N x 2) array of data points.
% w: the nominal control value.
% xnom: the keyframe state.
% x0, zeta_0: the starting point for zeta computation for PIPM & PPM.
%
% OUTPUT
% -----------
% sigma/zeta: An (N x 2) array of states in manifold space.
%

switch mode
    case 'pipm'
        sigma= (xnom(2)/w)^2 * (x(:,2).^2-xnom(2)^2-w^2*(x(:,1)-xnom(1)).^2);
        zeta= zeta_0 + zeta_0*(x(:,2)./x0(2)).^(w^2).*(x(:,1)-xnom(1))./x0(1);
        % zeta= zeta_0 + zeta_0*exp(w^2*log(x(:,2)./x0(2))).*(x(:,1)-xnom(1))./(x0(1)-xnom(1));
    case 'ppm'
        sigma= -(xnom(2)/w)^2 * (x(:,2).^2-xnom(2)^2+w^2*(x(:,1)-xnom(1)).^2);
        zeta= zeta_0 + zeta_0*(x(:,2)./x0(2)).^(-w^2).*(x(:,1)-xnom(1))./x0(1);
    case 'mcm'
        sigma= 2*w*(x(:,1)-xnom(1))- (x(:,2).^2-xnom(2)^2);
        zeta= zeta_0 + zeta_0*((-w)*log(x(:,2)./xnom(2))-(x(:,1)-xnom(1)));
    otherwise
        sigma= xnorm(1);
        zeta= xnorm(2);
end