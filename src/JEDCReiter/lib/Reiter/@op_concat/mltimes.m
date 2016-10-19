function z = mltimes(S,v)

  if(isreal(v) | isa(v,'deriv1'))
    z = mltimes(S.b,mltimes(S.a,v));
  else
    assert(0);
  end
