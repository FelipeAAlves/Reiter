function J = getjac(x)
    if (isa(x.d,'ssparse'))
        if (prod(size(x.d))<=1000)
            J = full(x.d);
        else
            J = sparse(x.d);
        end
    else
        J = x.d;
    end
