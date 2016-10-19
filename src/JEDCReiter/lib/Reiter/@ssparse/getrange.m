function [iRC,V,nr] = getrowrange(X,ica,icb,ReduceRowIndex)
  if(any(ica<1 | icb>X.n(1)*X.n(2)))
    error('index out of range');
  end

  iRC=[]; V=[];nr=0;
  for j=1:length(ica)
    ia = firstindexge(X.i,ica(j));
    ib = lastindexle(X.i,icb(j));
    ii = ia:ib;
    if(ib>=ia)
      iRC = [iRC;X.i(ii)-ReduceRowIndex(j)];
      V = [V;X.v(ii)];
    end
    nr = nr + icb(j) - ica(j) + 1;
  end

