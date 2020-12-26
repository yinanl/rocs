function p = labeling(x, G, A)

if(x(1)>=G(1,1) && x(1)<=G(1,2) &&...
        x(2)>=G(2,1) && x(2)<=G(2,2))
    p= 1;
elseif(x(1)>=A(1,1) && x(1)<=A(1,2) &&...
        x(2)>=A(2,1) && x(2)<=A(2,2))
    p= -1;
else
    p= 0;
end

end