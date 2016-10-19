% L*L' = S
function L = sqrt(S)
  L = S.U.*repmat(S.s,size(S.U,1),1);