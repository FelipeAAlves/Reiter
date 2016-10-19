function sout=uminus(s1)
sout.v=-s1.v;
sout.d=-s1.d;
sout=deriv1(sout);
