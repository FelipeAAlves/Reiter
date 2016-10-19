function z = mtimes(x,v)
  if(isreal(v) || isa(v,'deriv1'))
    assert(size(v,1)==size(x,2));
    z = initsize(sparse(size(x,1),size(v,2)),x.data{:},v);
    i1 = 0;
    i2 = 0;
    for i=1:length(x.data)
      [nr,nc] = size(x.data{i});
      z(i1+1:i1+nr,:) = x.data{i}*v(i2+1:i2+nc,:);
      i1 = i1+nr;
      i2 = i2+nc;
    end
  else
    error('wrong second argument type in mtimes');
  end