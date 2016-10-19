function xOut = sum(x,isum)
  n = size(x);
  if(nargin<2)
    if(n(1)==1 & n(2)>1)
      isum = 2;
    else
      isum = 1;
    end
  end
  [i,j,v] = find(x);
  if(isum==2)
    j = ones(size(i));
    xOut = ssparse(i,j,v,n(1),1);
  else
    i = ones(size(j));
    xOut = ssparse(i,j,v,1,n(2));
  end
