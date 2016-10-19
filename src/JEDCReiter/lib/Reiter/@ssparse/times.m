function z = times(x,y)
  if(isscalar(x))
    x = full(x);
    if(isa(y,'ssparse'))
      z = y;
      z.v = y.v*x;
    else
      z = y*x;
    end
  elseif(isscalar(y))
    y = full(y);
    if(isa(x,'ssparse'))
      z = x;
      z.v = x.v*y;
    else
      z = y*x;
    end
  else
    nx = size(x);
    ny = size(y);
    if(any(nx~=ny))
      error('incompatible size in times');
    end
    if(~isa(x,'ssparse'))
      x = ssparse(x);
    end
    if(~isa(y,'ssparse'))
      y = ssparse(y);
    end
    [irx,icx,vx] = find(x);
    [iry,icy,vy] = find(y);
    [ii,ix,iy] = intersect(x.i,y.i);
    z = ssparse(ii,[],x.v(ix).*y.v(iy),nx(1),nx(2));
  end