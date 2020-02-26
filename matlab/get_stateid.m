function ind= get_stateid(x, xgrid)
% FUNCTION
% -----------
% Obtain the index in a grid of a state.
%
% INPUT
% -----------
% x:        An (1 x n) state 
% xgrid:    A (N x n) state grid. N is the number of states, n is dimension
% OUTPUT
% -----------
% ind: the index of x in xgrid.
%

dx= xgrid - repmat(x,size(xgrid,1),1);
[~, ind]= min(sum(dx.^2,2));
end