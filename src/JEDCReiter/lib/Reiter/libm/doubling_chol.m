% solve for square root of discrete Lyapunof equation
% S = L*L' solves
% S = A*S*A' + B0*B0'
function L = doubling_chol(A,B0,n)

  nh = size(B0,2);
  nA = size(A,2);

  tol = 1e-14;
  % B'*Pi = Q*R;
  B = B0;
  for i=1:n
    if(max(abs(A*ones(nA,1)))<1e-50)
      break;
    end
    if(1)
      [Q,R,p,r] = rrqr(B',tol);
      % B(p,:)' = Q*R
      % disp(distance(B(p,:)',Q*R));
      % disp(distance(B',Q*R(:,pinv)))
      pinv = iperm(p);
      L = R(1:r,pinv)';
      %disp(distance(L*L',B*B'));
    else
      L = B;
    end
    if(i<n)
      B = [L A*L];
      A = A*A;
    end
  end


% inverts a permutation:
function ip = iperm(p)
  ip = zeros(size(p));
  ip(p) = 1:length(p);
