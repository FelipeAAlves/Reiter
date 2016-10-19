
function xOut=plus(x1,x2)
  n1 = size(x1);
  n2 = size(x2);
  if(prod(n1)==1 & prod(n2)>1)
      x1 = repmat(x1,n2(1),n2(2));
  end
   if(prod(n2)==1 & prod(n1)>1)
      x2 = repmat(x2,n1(1),n1(2));
  end
  x = getval(x1);
  y = getval(x2);
  dfdx1 = 1;
  dfdx2 = 1;
  dfdx1 = dfdx1(:);
  dfdx2 = dfdx2(:);
  if(notisa(x1,'deriv1'))
    np = nindep(x2);
  else
    np = nindep(x1);
  end
  if(notisa(x1,'deriv1'))
    xOut.v = plus(x1,x2.v);
    if isscalar(x2)
        d2 = repmat(x2.d, n1(1),n1(2));
    else
        d2 = x2.d;
    end
    xOut.d = elelmult_eachcol(dfdx2,d2);        %falves CHANGED
  elseif(notisa(x2,'deriv1'))
    xOut.v = plus(x1.v,x2);
    if isscalar(x1)
        d1 = repmat(x1.d, n2(1), n2(2));
    else
        d1 = x1.d;
    end
    xOut.d = elelmult_eachcol(dfdx1,d1);
  else  % both deriv1
    xOut.v = plus(x1.v,x2.v);
    if isscalar(x1) && ~isscalar(x2)
        d1 = repmat(x1.d, n2(1), n2(2));
        d2 = x2.d;
    elseif ~isscalar(x1) && isscalar(x2)
        d1 = x1.d;
        d2 = repmat(x1.d, n1(1), n1(2));
    else
        d1 = x1.d;
        d2 = x2.d;
    end
    xOut.d = elelmult_eachcol(dfdx1,d1) + elelmult_eachcol(dfdx2,d2);
  end
  xOut=deriv1(xOut);
