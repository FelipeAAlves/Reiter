% Matlab file to implement class "deriv1",
% forward mode of automatic differentation, first derivatives, possibly ssparse
% Michael Reiter, Universitat Pompeu Fabra, April 2007
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
% 
function xOut=horzcat(varargin)
  nv = length(varargin);
  for i=1:nv
    if(isa(varargin{i},'deriv1'))
      np = nindep(varargin{i});
      break;
    end
  end
  istart = 0;
  for i=1:nv
    if(~isa(varargin{i},'deriv1'))
      varargin{i} = deriv1(varargin{i},0,np);
    end
    if(~isempty(varargin{i}.v) & istart==0)
      istart = i;
    end
  end
  if(istart==0)
    xOut = [];
    return;
  else
    x=varargin{istart};
    newv = x.v;
    newd = x.d;
    [nr,nc,np] = size(x);
    for i=istart+1:nargin
      s2=varargin{i};
      if(isempty(s2.v))
	continue;
      end
      [nr2,nc2] = size(s2);
      nc = nc+nc2;
      if(nr2~=nr & ~isempty(newv))
	error('unequal row number in horizontal concatenation');
      end
      newv=[newv s2.v];
      newd=[newd;s2.d];
    end
  end
  xOut.v = newv;
  xOut.d = newd;
  xOut=deriv1(xOut);


