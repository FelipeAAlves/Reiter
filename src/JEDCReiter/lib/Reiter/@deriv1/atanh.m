
function xOut=atanh(xIn)
  x = xIn.v;
  xOut.v=atanh(x);
  d = 1./(1-x.^2);
  xOut.d = elelmult_eachcol(d(:),xIn.d);

  xOut=deriv1(xOut);
