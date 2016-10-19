function [xOut,I]=min(x1,x2)
  if(nargin==1)
    [xOut.v,I] = min(getval(x1));
    [nr,nc] = size(x1.v);
    indx = index([nr nc],[I(:) (1:nc)']);
    xOut.d = x1.d(indx,:);
  else
    x1v = getval(x1);
    x2v = getval(x2);
    xOut.v = min(x1v,x2v);
    i1 = x1v<=x2v;
    i2 = ~i1;
    if(isa(x1,'deriv1'))
      xOut.d = x1.d;
      if(isa(x2,'deriv1'))
	xOut.d(i2,:) = x2.d(i2,:);
      else
	xOut.d(i2,:) = 0;
      end
    else
      xOut.d = x2.d;
      xOut.d(i1,:) = 0;
    end
  end
  xOut=deriv1(xOut);
