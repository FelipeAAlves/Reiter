function ind=lookup(tabvals,x,varargin)
  ind = lookup(getval(tabvals),getval(x),varargin{:});
