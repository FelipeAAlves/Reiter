
function xOut=exp(xIn)
  x = xIn.v;
  xOut.v=exp(x);
  d = exp(x);
  xOut.d = elelmult_eachcol(d(:),xIn.d);

  xOut=deriv1(xOut);
