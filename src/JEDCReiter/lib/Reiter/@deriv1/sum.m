function xOut=sum(xIn,is)
  n = size(xIn);
  if(nargin==1)
    is = 1;
    if(n(1)==1 & n(2)>1)
      is = 2;
    end
  end
  if(is>2 | is<1)
    error('sum not implemented for_ dimension>2');
  end
  xOut.v=sum(xIn.v,is);
  if(n(3-is)==1)  %is is only non-singleton dimension
    xOut.d = sum(xIn.d);
  else
    [nr nc np] = size(xIn);
    [i,j,v] = find(xIn.d);
    [ir,ic] = ind2sub2([nr nc],i);
    if(is==1)
      xOut.d = ssparse(ic,j,v,nc,np);
    else  %if(is==2)
      xOut.d = ssparse(ir,j,v,nr,np);
    end
    if(~dissparse(xIn))
      xOut.d = full(xOut.d);
    end
  end
  xOut=deriv1(xOut);
