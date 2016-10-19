
function xOut=polyval(p,xIn)
%added Feb 3, 2012 by Alisdair McKay

  x = xIn.v;
  xOut.v=polyval(p,x);
  k = polyder(p);
  d = polyval(k,x);
  xOut.d = elelmult_eachcol(d(:),xIn.d);

  xOut=deriv1(xOut);
