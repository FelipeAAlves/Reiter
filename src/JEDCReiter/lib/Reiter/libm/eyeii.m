function x = eyeii(n,i);
  assert(n>0);
  assert(all(i>0 & i<=n));
  m = length(i);
  x = sparse(i,(1:m)',ones(size(i)),n,m);
