function [xv, xl, table] = grid_uniform_generate(area, gw)
% GRID_UNIFROM_GENERATE - generate a n-dimensinal uniform grid space
% Two steps:
%   1) converts the value to the label for each dimension:
%       label [xl1, xl2, ..., xln] ([N1, N2,..., Nn])
%   2) relates each combination of labels to a unique index by
%       ind=x1 + N1(x2-1) + N1*N2(x3-1)... + N1*...*Nn-1(xn-1)
%
% SYNTEX
% ---------
% [xv, xl, table] = grid_uniform_generate(area, gw)
% 
% INPUTS
% ---------
% area:         (n x 2) [min, max]
% gw:           (n x 1) grid width
%
% RETURNS
% ---------
% xv:           indiced array with real value
% xl:           indiced array with index in each dimension
% table
%   .xmin:      (n x 1) the minimum value of the table
%   .sizeDim:   (n x 1)
%   .gw:        (n x 1) grid width
%   .clv:       (n x 1) column of factors converting xl to index


xv = [];
xl = [];
table = [];
% check input area
if( any(area(:,1) > area(:,2)) )
    error('grid_uniform_subset:InputIllegal','area=[min, max].'); 
end

n = size(area, 1);

% check input grid width
if(size(gw, 1) == 1 && n > 1)
    gw = gw * ones(n, 1);
elseif(size(gw, 1) ~= n)
    error('dimensions of two inputs should match.');
end



% compute the number of grids in each dimension (sizeDim)
% and the overall number of grids in the space (sizeAll)
sizeDim = floor( (area(:,2)-area(:,1))./gw ) + 1;
sizeAll = prod(sizeDim);



% xmin = xv[1] = a + mod(b-a, gw)/2, area = [a, b]
xmin = mod(area(:,2)-area(:,1), gw) / 2 + area(:,1);



% relate each index to corresponding n dimensional labels(start from 1)
xl = (f(ones(n,1), sizeDim))';

% xv = xl * clv
clv = ones(n,1);
for i = 2 : n
    clv(i) = clv(i - 1) * sizeDim(i - 1);
end

% xl --> xv  (labels to real values)
% x[n] = x[1] + (n-1) * gw, n=1,2,...num_v
xv = repmat(xmin',sizeAll,1) + (xl - 1) .* repmat(gw',sizeAll,1);

% table
table.sizeDim = sizeDim;
table.xmin = xmin;
table.xb = area;
table.gw = gw;
table.clv = clv;


function q = f(xlb, xub)
% traverse all the combination between xlb and xub
% recursive method
%
% @params:  xlb     lower bound (column)
%           xub     upper bound (column)
% @returns: q       a row of all the combinations

n = size(xlb, 1);

q = [];
if( n > 1 )
    
    xlb_sub = xlb(1: end-1);
    xub_sub = xub(1: end-1);
    
    q_sub = f(xlb_sub, xub_sub);
    nq = size(q_sub, 2);
    
    for i = xlb(end): 1: xub(end)
        q_a = [q_sub; i * ones(1,nq)];
        q = [q, q_a];
    end
    
else
    
    q = xlb: 1: xub;
    
end
