function X = setrange(X,ica,icb,ReduceRowIndex,Y)
  if(any(ica<1 | icb>X.n(1)*icb>X.n(2)))
    error('index out of range');
  end

  iDrop = []; nr=0;
  for j=1:length(ica)
    % eliminate old values:
    if(length(X.i)>0)
      ia = firstindexge(X.i,ica(j));
      ib = lastindexle(X.i,icb(j));
      if(ib>=ia)
	iDrop = [iDrop ia:ib];
      end
    end
    nr2 = icb(j) - ica(j) + 1;
    iay(j) = nr+1;
    iby(j) = nr+nr2;
    nr = nr + nr2;
  end
  iDrop = sort(iDrop);
  if(any(diff(iDrop(:))==0))
    error('multiple assignment of rows');
  end
  Keep = true(size(X.i));
  Keep(iDrop)=0;
  if(isscalar(Y) & full(Y)==0)
    iRC = [X.i(Keep)];
    V =  [X.v(Keep)];
  else
    nn = size(Y,1)*size(Y,2);
    if(isscalar(Y))
      Y = ssparse(1:nn,[],Y*ones(nn,1),nn,1);
    end
    if(~isa(Y,'ssparse'))
      Y = ssparse(Y);
    end
    if(nr~=nn)
      error('subscripted assignment dimension mismatch');
    end
    [ircY,vY] = getrange(Y,iay,iby,-ReduceRowIndex);
    iRC = [ircY;X.i(Keep)];
    V =  [vY ;X.v(Keep)];
  end
  X = ssparse(iRC,[],V,X.n(1),X.n(2));

