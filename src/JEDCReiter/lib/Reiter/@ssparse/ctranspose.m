function z = ctranspose(x)
  [irx,icx] = find(x);
  z = ssparse(icx,irx,x.v,x.n(2),x.n(1));

