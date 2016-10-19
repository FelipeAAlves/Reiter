function z = mtimes(x,y)
  if(isa(x,'op_concat'))
    z = x.a*(x.b*y);
  else
    z = (x*y.a)*y.b;
  end
