function u= ctlr_feasible(x, ctree, cindex, cvalue)
% Read current state, return admissible control ids

% search control tree ctree
i=1;
while(ctree(i,1)>=0)
    
    d= ctree(i,1)+1;
    if(x(d)>ctree(2*i,2*d) && x(d)<ctree(2*i,2*d+1))
        i= 2*i;
    else
        i= 2*i+1;
    end
end

% get the control
us= cvalue(cindex(:,1)==i,:);
u= find(us);