function m = min(M)
  n = size(M);
  if(n(2)==1 | n(1)==1)
    m = min(M.v);
    if(nnz(M)<n(1))
      m = min(m,0);
    end
  else
    error('not yet implemented');
  end