
function xOut=sin(xIn)
  x = xIn.v;
  xOut.v=sin(x);
  d = cos(x);
  xOut.d = elelmult_eachcol(d(:),xIn.d);

  xOut=deriv1(xOut);
