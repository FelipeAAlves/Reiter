function str = svdstruct(A)
[str.u,s,str.v] = svd(A,'econ');
str.s = diag(s)';
