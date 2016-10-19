% Matlab file to implement class "deriv1",
% forward mode of automatic differentation, first derivatives, possibly ssparse
% Michael Reiter, Universitat Pompeu Fabra, April 2007
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
function xOut=vertcat(varargin)
  nv = length(varargin);
  for i=1:nv
    if(isa(varargin{i},'deriv1'))
      np = nindep(varargin{i});
      break;
    end
  end
  istart = 0;
  iSparse = 0;
  for i=1:nv
    if(~isa(varargin{i},'deriv1'))
      varargin{i} = deriv1(varargin{i},0,np);
    end
    if(~isempty(varargin{i}.v) & istart==0)
      istart = i;
      nc = size(varargin{i}.v,2);
    end
    if(issparse(varargin{i}.d))
      iSparse = max(1,iSparse);
    end
    if(isssparse(varargin{i}.d))
      iSparse = max(2,iSparse);
    end
  end

  iR=[];iC=[];iJ=[];V=[];newv=[];newd=[];
  if(istart==0)
    xOut = [];
    return;
  else
    x=varargin{istart};
    nr = 0;
    for i=istart:nargin
      s2=varargin{i};
      if(isempty(s2.v))
	continue;
      end
      [nr2,nc2] = size(s2);
      if(nc2~=nc & ~isempty(newv))
	error('unequal column number in vertical concatenation');
      end
      newv=[newv;s2.v];
      newd2 = s2.d;

      [ii,ij,v] = find(ssparse(newd2));
      irc = index([nr2 nc],ii(:));
      iR = [iR;irc(:,1)+nr];
      iC = [iC;irc(:,2)];
      iJ = [iJ;ij];
      V = [V;v];

      nr = nr+nr2;
    end
  end
  %ii = index([nr nc],[iR iC]);
  ii = sub2ind([nr nc],iR, iC);
  if(iSparse==2)
    newd = ssparse(ii,iJ,V,nr*nc,np);
  else
    newd = sparse(ii,iJ,V,nr*nc,np);
    if(iSparse==0)
      newd = full(newd);
    end
  end

  xOut.v = newv;
  xOut.d = newd;
  xOut=deriv1(xOut);
