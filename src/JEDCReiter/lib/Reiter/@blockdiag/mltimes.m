function z = mltimes(x,v)
  if(isreal(v) | isa(v,'deriv1'))
    assert(size(v,1)==size(x,2));
    z = initsize(zeros(size(x,2),size(v,2)),x.data{:},v);
    i1 = 0;
    i2 = 0;
    for i=1:length(x.data)
      [nc,nr] = size(x.data{i});  % interchange rows and columns
      z(i1+1:i1+nr,:) = mltimes(x.data{i},v(i2+1:i2+nc,:));
      i1 = i1+nr;
      i2 = i2+nc;
    end
  else
    error('wrong second argument type in mtimes');
  end