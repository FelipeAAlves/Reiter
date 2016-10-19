function [x2,p2] = ptd_merge(x,p)
if(iscell(p))
  error('not yet impl.');
else
  if(size(p,2)>1 && size(x,2)==1)
    x2 = x;
    p2 = sum(p,2);
  else
    x = x(:);
    p = p(:);
    [x,I] = sort(x);
    p = p(I);
    x2 = x;
    p2 = p;
    iKeep = true(size(x2));
    j = 1;  %currently valid i;
    for i=2:length(x2)
      if(x2(i)==x2(i-1))
        p2(j) = p2(j) + p2(i);  % add mass of i to j, and drop i;
        iKeep(i) = false;
      else
        j = i;
      end
    end
    x2 = x2(iKeep);
    p2 = p2(iKeep);
  end
end