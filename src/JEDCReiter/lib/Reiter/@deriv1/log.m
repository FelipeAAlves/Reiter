
function xOut=log(xIn)
  x = xIn.v;
  xOut.v=log(x);
  d = 1./x;
  xOut.d = elelmult_eachcol(d(:),xIn.d);

  xOut=deriv1(xOut);
