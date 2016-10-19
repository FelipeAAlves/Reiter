
function xOut=tanh(xIn)
  x = xIn.v;
  xOut.v=tanh(x);
  d = sech(x).^2;
  xOut.d = elelmult_eachcol(d(:),xIn.d);

  xOut=deriv1(xOut);
