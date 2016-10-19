% Matlab file for_ paper 
% "Approximate and Almost-Exact Aggregation in Dynamic Stochastic Heterogeneous-Agent Models"
% Michael Reiter, IHS, December 2009
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% computes "Gramian", solution of discrete Lyapunov equation
%   G = A*G*A' + B*B'
% Inputs:
%    A, B as above
%    n: horizon for gramianspace
%    eps: threshold to determine rank
% Outputs:
%   gstr:  Gramian structure with SVD of observability matrix
%   L:     square root of Gramian
function [G,GramSpace] = gramian(A,B,n,eps)
  % transform structures into gramian class:
  if(nargin==1 && isstruct(A))
    G = class(A,'gramian');
    return;
  end
  [rb,cb] = size(B);
  % transform U and s into gramian class:
  if(rb==1 && cb==size(A,2))
    G.U = A;
    G.s = B;
    G = class(G,'gramian');
    return;
  end
  
  % from here on do real computations:
  if(nargin<4)
    eps = 1e-14;
  end
  assert(size(A,1)==size(A,2) & size(A,1)==size(B,1));
  if(n==0 || n>=size(A,1))
    if(1)
      L = doubling_chol(A,B,30);
      [G.U,s] = svd(L,'econ');
      G.s = diag(s)';
      G = class(G,'gramian');
    else
      S = doubling(A,B*B',30);
      S = (S+S')/2;
      [G.U,s] = svd(S);
      G.s = sqrt(diag(s)');
      G = class(G,'gramian');      
    end
    if(nargout>1)
      GramSpace = gramianspace_ex(A,B);    
    end
  else
    GramSpace = gramianspace(A,B,n);
    G = gramian_lowrank(A,B,GramSpace,eps);
  end

