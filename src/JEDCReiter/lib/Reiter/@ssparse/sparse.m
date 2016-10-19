function M = sparse(S)
  [nr,nc] = size(S);
  [ir,ic,v] = find(S);
  M = sparse(ir,ic,v,nr,nc);
