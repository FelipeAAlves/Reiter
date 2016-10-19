
function xOut=rdivide(x1,x2)
  n1 = size(x1);
  n2 = size(x2);
  x = getval(x1);
  y = getval(x2);
  dfdx1 = 1./y;
  dfdx2 = -x./y.^2;
  dfdx1 = dfdx1(:);
  dfdx2 = dfdx2(:);
  if(notisa(x1,'deriv1'))
    np = nindep(x2);
  else
    np = nindep(x1);
  end
  if(notisa(x1,'deriv1'))
    xOut.v = rdivide(x1,x2.v);
    xOut.d = elelmult_eachcol(dfdx2,x2.d);
  elseif(notisa(x2,'deriv1'))
    xOut.v = rdivide(x1.v,x2);
    xOut.d = elelmult_eachcol(dfdx1,x1.d);
  else  % both deriv1
    n1d = size(x1.d);
    np = n1d(end);
    xOut.v = rdivide(x1.v,x2.v);
    xOut.d = elelmult_eachcol(dfdx1,x1.d) + elelmult_eachcol(dfdx2,x2.d);
  end
  xOut=deriv1(xOut);
