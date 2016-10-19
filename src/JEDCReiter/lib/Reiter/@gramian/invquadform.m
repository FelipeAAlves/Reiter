function q = invquadform(H,S,H2,eps)
  s = S.s;
  if(nargin>3)
    s = max(s,s(1)*eps);
  end
  US = S.U.*repmat(1./s,size(S.U,1),1);
  if(nargin>2 & ~isempty(H2))
    q = H*US*(US'*H2);
  else
    q = H*US*(US'*H');
  end
