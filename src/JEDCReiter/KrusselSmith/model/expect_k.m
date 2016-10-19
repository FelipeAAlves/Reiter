

%%% Description
%       Compute j-th central moment of the cross-sectional distribution of k,
%
%%% INPUTS:
%   pvec:  vector of probabilities
%          with 2 columns: first column vector of k-nodes, second column vector of probabilities
%   j:     which_ moment
function EV = expect_k(pvec,j,k)
    global MP;
    if (nargin==1)
        j = 1;
    end

    if (nargin<3)
        k = MP.AssetsGridFine;
    end
    doCentral = 0;

    if (j<-1)
        doCentral = 1;
    end
    j = abs(j);


    subtr = 0;
    if (doCentral) %compute central moment
        subtr = expect_k(pvec,1,k);
    end

    % subtract mean if j>1:
    k = k-subtr;
    dy = k.^j;
    dy = repmat(dy,MP.neps,1);
    if (isempty(pvec))  %return integrating vector:
        EV = dy;
    else
        EV = dot(dy,pvec);
    end
