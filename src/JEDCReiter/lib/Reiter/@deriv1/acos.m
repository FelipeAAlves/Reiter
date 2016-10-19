
function xOut=acos(xIn)
  x = xIn.v;
  xOut.v=acos(x);
  d = -1./sqrt(1-x.^2);
  xOut.d = elelmult_eachcol(d(:),xIn.d);

  xOut=deriv1(xOut);
