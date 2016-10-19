function xOut=erf(xIn)
  x = xIn.v;
  xOut.v=erf(x);
  d = (2/sqrt(pi))*exp(-x.^2);
  xOut.d = elelmult_eachcol(d(:),xIn.d);

  xOut=deriv1(xOut);
