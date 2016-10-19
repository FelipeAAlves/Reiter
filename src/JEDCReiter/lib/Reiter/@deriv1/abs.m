
function xOut=abs(xIn)
  x = xIn.v;
  xOut.v = abs(x);

  d = ones(size(x));
  d(x<0) = -1;
  d(x==0) = NaN;
  xOut.d = elelmult_eachcol(d(:),xIn.d);
  xOut = deriv1(xOut);
