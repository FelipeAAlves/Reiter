
function xOut=atan(xIn)
  x = xIn.v;
  xOut.v=atan(x);
  d = 1./(x.^2+1);
  xOut.d = elelmult_eachcol(d(:),xIn.d);

  xOut=deriv1(xOut);
