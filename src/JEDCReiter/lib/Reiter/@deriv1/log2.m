
function xOut=log2(xIn)
  x = xIn.v;
  xOut.v=log2(x);
  d = 1./(log(2).*x);
  xOut.d = elelmult_eachcol(d(:),xIn.d);

  xOut=deriv1(xOut);
