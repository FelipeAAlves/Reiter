% Matlab file for_ paper 
% "Approximate and Almost-Exact Aggregation in Dynamic Stochastic Heterogeneous-Agent Models"
% Michael Reiter, IHS, December 2009
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% SuperSparse matrix class
function S = ssparse(ir,ic,v,nr,nc)
  if(nargin==1)
    if(isa(ir,'ssparse'))
      S = ir;
      return;
    end
    n = size(ir);
    if(length(n)~=2)
      error('ssparse must be matrix');
    end
    S.i = find(ir);  %also handles case_ of sparse (not ssparse) ir!
    S.v = full(ir(S.i));
    S.n = n;
  elseif(nargin==2)
    S.i = [];
    S.v = [];
    S.n = [ir ic];
  else
    if(isempty(ic))
      ii = ir(:);
      if(any(ii<1) | any(ii>nr*nc))
	error('wrong size in ssparse');
      end
    else
      if(any(ir<1) | any(ir>nr) | any(ic<1) | any(ic>nc))
	error('wrong size in ssparse');
      end
      ii = sub2ind2([nr nc],ir(:),ic(:));
    end
    [ii,indx] = sort(ii);
    v = v(indx);
    n = length(ii);
    same = find(diff(ii)==0);
    if(isempty(same))
      S.i = ii;
      S.v = v;
    else
      for i=1:length(same)
	j = same(i);
	ii(j)=0;
	v(j+1) = v(j+1)+v(j);
      end
      keep = (ii~=0);
      S.i = ii(keep);
      S.v = v(keep);
    end
    keep = (S.v~=0);
    S.i = S.i(keep);
    S.v = S.v(keep);
    S.n = [nr nc];
  end
  S = class(S,'ssparse');
