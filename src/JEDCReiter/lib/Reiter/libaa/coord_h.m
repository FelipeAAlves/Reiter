function C = coord_h(H,Hbasic,eps)
if(nargin<3)
  eps = 1e-10;
end
C = H\Hbasic;
assert(distance(Hbasic,H*C,1)<eps);

