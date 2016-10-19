function q = quadform(H,S,H2)
  if(nargin==3)
    q = H*S*H2;
  else
    q = H*S*H';
  end
