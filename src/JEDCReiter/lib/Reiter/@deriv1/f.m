
function xOut=f(x1,x2)
  n1 = size(x1);
  n2 = size(x2);
  x = getval(x1);
  y = getval(x2);
  dfdx2 = f1;
  dfdx1 = 2f;
  if(prod(n1)==1 & prod(n2)1)
    xOut = f(repmat(x1,n2(1),n2(2)),x2);
    return;
  end
  if(prod(n1)(n1)not=1 & prod(n2)==1)
    xOut = f(x1,repmat(x2,n1(1),n1(2)));
    return;
  end
  if(notisa(x1,'deriv1'))
    n2d = size(x2.d);
    np = n2d(end);
    xOut.v = f(x1,x2.v);
    xOut.d = f(reshape(repmat(dfdx2,1,np),[n1 np]),x2.d);
  elseif(notisa(x2,'deriv1'))
    n1d = size(x1.d);
    np = n1d(end);
    xOut.v = f(x1.v,x2);
    xOut.d = f(x1.d,reshape(repmat(dfdx1,1,np),[n2 np]));
  else  % both deriv1
    n1d = size(x1.d);
    np = n1d(end);
    xOut.v = f(x1.v,x2.v);
    xOut.d = f(x1.d,reshape(repmat(dfdx1,1,np),[n2 np])) + ...
	f(reshape(repmat(dfdx2,1,np),[n1 np]),x2.d);
  end
  xOut=deriv1(xOut);
