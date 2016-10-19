% solves A*X + B*X*T = C12;
% works best if dimension of T small!
function S = sylvschur(A,B,T,C12)

  [Q,R] = schur(T,'real');
  C = C12*Q;

  n = size(Q,1);
  m = size(B,1);
  i=1;
  S = zeros(size(C12));
  id2 = kron(eye(2),A);
  while(i<=n)
    if(i<n & R(i+1,i)~=0)
      ii = i:i+1;
      Binv = id2 + kron(R(ii,ii)',B);
    else
      ii = i;
      Binv = A + R(i,i)*B;
    end
    if(i==1)
      Ci = C(:,ii);
    else
      Ci = C(:,ii) - B*S(:,1:i-1)*R(1:i-1,ii);
    end
    if(length(ii)==1)
      S(:,ii) = Binv\Ci;
    else
      S(:,ii) = reshape(Binv\Ci(:),m,2);
    end
    i = i + length(ii);
  end

  S = S*Q';
