function sout=paren(s,indx)
  sout.v = s.v(indx);
  sout.d = s.d(indx,:);
  sout=deriv1(sout);

