function [R,nr] = index2ranges(ii)
  nr = 0;
  if(isempty(ii))
    R = [1 0 0];
    return;
  end
  if(all(diff(ii)==1))
    R = [ii(1) ii(end) ii(1)-1];
    nr = length(ii);
    return;
  end
  R = [];
  l = 0;
  n = length(ii);
  r = zeros(1,3);
  while(l<n)
    l = l+1;
    r(1) = ii(l);
    r(3) = ii(l)-l;
    while(l<n & (ii(l+1)-ii(l)==1) )
      l = l+1;
    end;
    r(2) = ii(l);
    nr = nr + r(2)-r(1)+1;
    R = [R;r];
  end
