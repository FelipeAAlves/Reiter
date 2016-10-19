function q = full(S)
  US = S.U.*repmat(S.s,size(S.U,1),1);
  q = US*US';
