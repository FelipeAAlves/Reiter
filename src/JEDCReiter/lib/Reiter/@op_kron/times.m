function z = times(x,y)
  if(isa(x,'op_kron'))
    z = x;
    z.a = z.a*y;
  end
  if(isa(y,'op_kron'))
    z = y;
    z.a = z.a*x;
  end

