function S = getcolrange(X,ica,icb)
  if(nargin==2)
    if(any(ica<1) | any(ica>X.n(2)))
      error('index out of range');
    end
    [ir,ic] = ind2sub2(X.n,X.i);
    assert(0);
    ia = firstindexge(ic,ica);
    ib = lastindexle(ic,icb);
    if(ib<ia)
      S = ssparse(0,0);
    else
      ii = ia:ib;
      S = ssparse(ir(ii),ic(ii)-(ica-1),X.v(ii),X.n(1),icb-ica+1);
      %S.v = X.v(ii);
      %S.n = [X.n(1) icb-ica+1];
      %S.i = sub2ind2(S.n,ir(ii),ic(ii)-(ica+1));
      %S = class(S,'ssparse');
    end 
    return;
  end

  if(ica<1 | icb>X.n(2))
    error('index out of range');
  end
  [ir,ic] = ind2sub2(X.n,X.i);
  ia = firstindexge(ic,ica);
  ib = lastindexle(ic,icb);
  if(ib<ia)
    S = ssparse(0,0);
  else
    ii = ia:ib;
    S = ssparse(ir(ii),ic(ii)-(ica-1),X.v(ii),X.n(1),icb-ica+1);
  end



    