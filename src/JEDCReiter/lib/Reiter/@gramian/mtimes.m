function q = mtimes(S,b)
  US = S.U.*repmat(S.s,size(S.U,1),1);
  q = US*(US'*b);

