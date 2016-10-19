function [n,varargout] = handle_size(nn,i)
  if(nargin>1)
    assert(nargout<=1);
    if(i<=0)
      error('wrong index in size');
    elseif(i>2)
      n = 1;
    else
      n = nn(i);
    end
  else
    if(nargout<=1)
      n = nn;
    else
      n = nn(1);
      varargout{1} = nn(2);
      for j=3:nargout
	varargout{j} = 1;
      end
    end
  end
