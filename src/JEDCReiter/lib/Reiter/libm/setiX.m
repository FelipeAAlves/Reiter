%%% Description
%       Set indexes for computing derivatives in Reiter case:
%
%%% INPUT
%       (1) Number of variables
%       (2) Number of expectational shocks
%       (3) Number of shocks
%%% OUTPUT:
%    X2  : vector with ALL variables (including shocks)
%    indx: index structure (where in the vector the variables are)
%
function [X2,indx] = setiX(X, nEta, nEps)

    nx     = length(X);

    % variables indexes
    indx.x    = 1:nx;
    indx.xlag = nx + indx.x;

    nn = max(indx.xlag);

    indx.eta = nn+1:nn+nEta;
    indx.eps = nn+nEta+1:nn+nEta+nEps;

    % nX2 = max([max(indx.x), max(indx.xlag), max(indx.eps), max(indx.eta)]);
    X2 = zeros(nn+nEta+nEps,1);
    X2([indx.x indx.xlag]) = repmat(X,2,1);
