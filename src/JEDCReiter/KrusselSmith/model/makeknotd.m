function [knots,logshift] = makeknotd(kmin,kmax,n,logshift)
    if(nargin>3)
        knots = logspaceshift(kmin,kmax,n,logshift);
    else
        [knots,logshift] = logspaceshift(kmin,kmax,n,1,n/4);
    end
    
    knots = knots';
