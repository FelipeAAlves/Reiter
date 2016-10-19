
function xOut=log10(xIn)
  x = xIn.v;
  xOut.v=log10(x);
  d = 1./(log(10).*x);
  xOut.d = elelmult_eachcol(d(:),xIn.d);

  xOut=deriv1(xOut);
