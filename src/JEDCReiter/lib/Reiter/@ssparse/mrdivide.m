function z = mrdivide(x1,x2)
  if(isscalar(x1) | isscalar(x2))
    z = rdivide(x1,x2);
    return;
  else
    assert(0);
  end