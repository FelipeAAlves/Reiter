function M = spdiag(d)
  n = length(d);
  n2 = n*n;
  M.v  = sparse(1:n,1:n,getval(d),n,n);
  [i,j,v] = find(getjacob(d));
  if(isempty(i))
    M = M.v;
  else
    i2 = index([n n],[i i]);
    M.d  = ssparse(i2,j,v,n2,nindep(d));
    M = deriv1(M);
  end
