
function xOut=power2(xIn)
  x = xIn.v;
  xOut.v = x.^2;
  d = 2.*x;
  xOut.d = elelmult_eachcol(d(:),xIn.d);

  xOut=deriv1(xOut);
