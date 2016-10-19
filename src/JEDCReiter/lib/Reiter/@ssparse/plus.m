function z = plus(x,y)
  nx = size(x);
  ny = size(y);
  if(prod(nx)==1)
    x = full(x);
    if(x==0)
      z = y;
    else
      y = full(y);
      z = x+y;
    end
    return;
  end
  if(prod(ny)==1)
    y = full(y);
    if(y==0)
      z = x;
    else
      x = full(x);
      z = x+y;
    end
    return;
  end

  [irx,icx,vx] = find(x);
  [iry,icy,vy] = find(y);
  z = ssparse([irx;iry],[icx;icy],[vx;vy],nx(1),nx(2));
  