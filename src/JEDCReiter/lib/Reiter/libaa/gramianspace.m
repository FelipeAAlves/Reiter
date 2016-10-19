% compute subspace in which the Gramian Q of (A,h) lies
% Q satisfies
%    Q = A*Q*A' + h*h';
% subspace str=(U,S) satisfies
%    Q = U S S' U'  with U'U = I
% diagonal s os S on output gives singular values of observability matrix;
%    to get singular values of Gramian; use gramien_lowrank
function str = gramianspace(A,h0,n)

  %if(n==0 || n*length(A)<1e7)
  %  str = gramianspace_ex(A,h0,n);
  %  return;
  %end
  maxColumn = round(6e7/size(A,1)); % maximal number of columns of basis; restrict for memory reasons;
                                    % changed from 2e7 to 6e7 on 8/7/12
  
  h = h0;
  if(1)
    n1 = ceil(n/4);
    n2 = n - 3*n1;
    nn = [n1 n1 n1 n2];
    steps = [1 2 4 8];
    steps2 = [repmat(steps(1),1,nn(1)) repmat(steps(2),1,nn(2)) ...
              repmat(steps(3),1,nn(3)) repmat(steps(4),1,nn(4))];
  else
    n1 = ceil(n/2);
    n2 = ceil(n/4);
    n3 = n - n1 - n2;
    nn = [n1 n2 n3];
    steps = [1 1 1];
    steps2 = [repmat(steps(1),1,nn(1)) repmat(steps(2),1,nn(2)) repmat(steps(3),1,nn(3))];
  end
  nh = size(h,2);
  QR = zeros(size(A,1),min(nh*n,50));
  QR(:,1:nh) = qr_house(h);
  offs = nh;
  norm0 = norm(h);
  iDone = 0;
  for iter=1:n
    if(iDone)
      break;
    end
    for l=1:steps2(iter)
      h = A*h;
    end
    if(max(abs(h))>max(abs(h0))*1e6)
      fprintf(1,'WARNING: stop in gramianspace because of likely instability of A');
      break;
    end
    for ih=1:nh
      x = update2_house(QR,offs,full(h(:,ih)),1e-15);
      if(~isempty(x))
	if(any(isfinite(x)==0))
	  error('NaN in x');
	end
	offs = offs+1;
	if(size(QR,2)<offs)
	  QR = morecols(QR,50);
	end
	QR(:,offs) = x;
        if(offs==size(QR,1))
          iDone = 1;
          break;
        end
        if(offs>=maxColumn)
          disp('Warning in gramianspace: computation finished for lack of space');
          iDone = 1;
          break;
        end
      end
    end
  end
  QR = QR(:,1:offs);
  [Q,R] = full_house(QR);
  str = svdstruct(R);
  str.u = Q*str.u;


function m2 = morecols(m,n)
  m2 = [m zeros(size(m,1),n)];
