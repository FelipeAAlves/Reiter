
function xOut=sinh(xIn)
  x = xIn.v;
  xOut.v=sinh(x);
  d = cosh(x);
  xOut.d = elelmult_eachcol(d(:),xIn.d);

  xOut=deriv1(xOut);
