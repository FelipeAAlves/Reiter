% compute low-rank Gramian in factorized form
%   str is output from gramianspace
function gstr = gramian_lowrank(A,B,str,eps)
  if(nargin<4)
    eps = 1e-15;
  end
  n = max(find(str.s>str.s(1)*eps));

  if(n<size(str.u,2))
    u = str.u(:,1:n);
    AU = A*u;
    UAU = (AU'*u)';
    lam = (B'*u)';
  else
    AU = A*str.u;
    UAU = (AU'*str.u)';
    lam = (B'*str.u)';
  end

  if(any(abs(eig(UAU))>=1))
    % warning('matrix not stable');
    error('matrix not stable');
    [vv,ee] = eig(UAU);
    ed = abs(diag(ee));
    fac = 0.999./max(ed,1);
    UAU = real(vv*diag(fac.*diag(ee))*inv(vv));
  end
  L = doubling_chol(UAU,lam,30);

  if(n<size(str.u,2))
    [gstr.U,s] = svd(u*L,'econ');  
  else
    [gstr.U,s] = svd(str.u*L,'econ');  
  end
  gstr.s = diag(s)';
  gstr = gramian(gstr);
  