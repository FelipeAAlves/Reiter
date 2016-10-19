function X = setrowrange(X,ica,icb,ReduceRowIndex,Y)
  if(any(ica<1 | icb>X.n(1)))
    error('index out of range');
  end
  [ir,ic] = ind2sub2(X.n,X.i);
  v = X.v;

  iDrop = []; nr=0;
  nI = length(ica);
  Keep = true(size(ir));
  for j=1:nI
    % eliminate old values:
    if(length(ir)>0)
      Keep = Keep & (ir<ica(j) | ir>icb(j));
    end
    nr2 = icb(j) - ica(j) + 1;
    iay(j) = nr+1;
    iby(j) = nr+nr2;
    nr = nr + nr2;
  end

  if(nI>1)
    Rsort = sortrows([ica icb],1);
    if(any(icb(1:nI-1)>=ica(2:nI)))
      error('multiple assignment of rows');
    end
  end


  if(isscalar(Y) & ~(nr==1 & size(X,2)==1))
    if(full(Y)==0)
      iR = [ir(Keep)];
      iC = [ic(Keep)];
      V =  [v(Keep)];
    else
      error('assigns scalar to the rows of sparse matrix');
    end
  else
    if(nr~=size(Y,1) | size(X,2)~=size(Y,2))
      error('subscripted assignment dimension mismatch');
    end
    [irY,icY,vY] = getrowrange(Y,iay,iby,-ReduceRowIndex);
    iR = [irY;ir(Keep)];
    iC = [icY;ic(Keep)];
    V =  [vY ;v(Keep)];
  end
  X = ssparse(iR,iC,V,X.n(1),X.n(2));

