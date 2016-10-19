function varargout = size(m,i)
  nn = size(m.v);
  if(nargout<=1)
    if(nargin>1)
      varargout = {nn(i)};
    else
      varargout = {nn};
    end
  else
    d = length(nn);
    varargout = cell(1,d);
    for i=1:d
      varargout{i} = nn(i);
    end
    if(nargout>d)
      varargout{d+1} = nindep(m);
    end
  end

    
