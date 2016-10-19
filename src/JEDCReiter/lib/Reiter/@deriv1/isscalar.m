function b = isscalar(x)
  n = size(x);
  b = prod(n)==1;