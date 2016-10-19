function X = doubling(A,Sig,n)
  A = A';

  for i=1:n-1
    Sig = Sig+A'*(Sig*A);
    A = A*A;
  end
  X = Sig+A'*(Sig*A);