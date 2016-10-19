function z = mltimes(S,v)

  if(isreal(v) | isa(v,'deriv1'))
    z = mltimes(S.a,v) + mltimes(S.b,v);
  else
    assert(0);
  end
