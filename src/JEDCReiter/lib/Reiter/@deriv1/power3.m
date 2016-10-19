
function xOut=power3(xIn)
  x = xIn.v;
  xOut.v=x.^3;
  d = 3.*x.^2;
  xOut.d = elelmult_eachcol(d(:),xIn.d);

  xOut=deriv1(xOut);
