
function xOut=asinh(xIn)
  x = xIn.v;
  xOut.v=asinh(x);
  d = 1./sqrt(x.^2+1);
  xOut.d = elelmult_eachcol(d(:),xIn.d);

  xOut=deriv1(xOut);
