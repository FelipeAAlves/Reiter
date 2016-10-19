
function xOut=tan(xIn)
  x = xIn.v;
  xOut.v=tan(x);
  d = sec(x).^2;
  xOut.d = elelmult_eachcol(d(:),xIn.d);

  xOut=deriv1(xOut);
