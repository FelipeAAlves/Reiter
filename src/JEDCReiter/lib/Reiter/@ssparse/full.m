function z = full(x);
  z = zeros(x.n);
  z(x.i) = x.v;
