function M = blockdiag(x,y)
  if(nargin>1)
    assert(nargin==2)
    M.data = {x,y};
  else
    assert(iscell(x));
    M.data = x;
  end

  M = class(M,'blockdiag');