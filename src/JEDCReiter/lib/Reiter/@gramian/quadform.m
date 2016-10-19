function q = quadform(H,S,H2)
  US = S.U.*repmat(S.s,size(S.U,1),1);
  if(nargin==3)
    q = H*US*(US'*H2);
  else
    q = H*US*(US'*H');
  end
