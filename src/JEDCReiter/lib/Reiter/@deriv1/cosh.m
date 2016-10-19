
function xOut=cosh(xIn)
  x = xIn.v;
  xOut.v=cosh(x);
  d = sinh(x);
  xOut.d = elelmult_eachcol(d(:),xIn.d);

  xOut=deriv1(xOut);
