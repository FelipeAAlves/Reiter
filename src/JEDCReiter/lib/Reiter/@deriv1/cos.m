
function xOut=cos(xIn)
  x = xIn.v;
  xOut.v=cos(x);
  d = -sin(x);
  xOut.d = elelmult_eachcol(d(:),xIn.d);

  xOut=deriv1(xOut);
