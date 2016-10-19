function [n,varargout] = size(S,i)
  if(nargin==1)
    if(nargout<=1)
      n = S.n;
    else
      n = size(S,1);
      for j=2:nargout
	varargout{j-1} = size(S,j);
      end
    end
  else
    if(i>2)
      n = 1;
    else
      n = S.n(i);
    end
  end
