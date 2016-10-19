% m is logical number of columns of QR
function x = update_house(QR,m,x,eps)

  n = size(QR,1);
  for k = 1:m
    v = [1;QR(k+1:n,k)];
    x(k:n) = householder_mult(x(k:n),v);
  end
  v = householder(x(m+1:n));
  if(isempty(v))
      x=[];
      return;
  end
  x(m+1:n) = householder_mult(x(m+1:n),v);

  rmax = max(diag(QR(1:m,1:m)));
  rmin = rmax*eps;
  if(abs(x(m+1))<rmin)  % x is almost linearly dependent with QR
    x = [];
  else
    x(m+2:n) = v(2:end);
  end



  