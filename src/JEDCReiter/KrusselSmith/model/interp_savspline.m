
%%% INPUT
%       (1) S  :        spline with savin fnc
%       (2) Xq :        nodes for the function to be eval
%%% OUTPUT
%   Yq : policy at nodes

function Yq = interp_savspline(S,Xq)

    %== INITIALIZE - allow to hold deriv1 ==%
    Yq = initsize(Xq, S.x, S.y);

    n = length(S.y);
    if (any(Xq>=S.x(n)))
        error('too big CoH');
    end

    Yq(Xq<S.x(1)) = 0;  % for lower CoH, no saving

    iMid = Xq>S.x(1);
    indXq = lookup(S.x,Xq(iMid),0);

    Yq(iMid) = S.y(indXq) + ( Xq(iMid) - S.x(indXq) ) .* S.Slope(indXq);
