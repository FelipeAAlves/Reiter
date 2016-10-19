function i = isvector(x)
  n = size(x);
  i = length(n)==2 & (n(1)==1 | n(2)==1);