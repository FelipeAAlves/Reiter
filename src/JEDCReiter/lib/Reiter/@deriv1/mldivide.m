function xOut=mldivide(x1,x2)
  if(all(all(tril(getval(x1),-1)==0)))  % upper triangular matrix
    xOut = x1*x2*0;  % INITIALIZATION; ACHTUNG: INEFFICIENT !!!???
    n = length(x1);
    m = size(x2,2);
    for i=n:-1:1
      xii = mysubsref(x1,'()',i,i);
      % rest = x2(i,1:m);
      rest = mysubsref(x2,'()',i,1:m);
      if(i<n)
	xin = mysubsref(x1,'()',i,i+1:n);
	rest = rest - xin*mysubsref(xOut,'()',i+1:n,1:m);
      end
      xOut = mysubsasgn(xOut,rest ./ repmat(xii,1,m),'()',i,1:m);
    end
  elseif(all(all(triu(getval(x1),1)==0)))  % lower triangular matrix
    xOut = x1*x2*0;  % INITIALIZATION; ACHTUNG: INEFFICIENT !!!???
    n = length(x1);
    m = size(x2,2);
    for i=1:n
      xii = mysubsref(x1,'()',i,i);
      % rest = x2(i,1:m);
      rest = mysubsref(x2,'()',i,1:m);
      if(i>1)
	xin = mysubsref(x1,'()',i,1:i-1);
	rest = rest - xin*mysubsref(xOut,'()',1:i-1,1:m);
      end
      xOut = mysubsasgn(xOut,rest ./ repmat(xii,1,m),'()',i,1:m);
    end
  else
    assert(0);
  end
