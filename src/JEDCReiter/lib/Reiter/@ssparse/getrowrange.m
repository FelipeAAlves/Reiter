function [iR,iC,V,nr] = getrowrange(X,ica,icb,ReduceRowIndex)
  if(any(ica<1 | icb>X.n(1)))
    error('index out of range');
  end
  if(nnz(X)==0)
    iR=[]; iC=[]; V=[];nr=0;
    return;
  end
  [ir,ic] = ind2sub2(X.n,X.i);

  nca = length(ica);
  iR=cell(nca,1); iC=cell(nca,1); V=cell(nca,1);
  nr = 0;
  for j=1:nca
    ii = ir>=ica(j) & ir<=icb(j);
    ii=ii(:);
    if(~isempty(ii))
      iR{j} = ir(ii)-ReduceRowIndex(j);
      iC{j} = ic(ii);
      V{j}  = X.v(ii);
    end
    nr = nr + icb(j) - ica(j) + 1;
  end
  iR=vcatcell(iR);
  iC=vcatcell(iC);
  V=vcatcell(V);

function X = vcatcell(x)
  n = length(x);
  nn = 0;
  for i=1:n
    nn = nn+length(x{i});
  end
  X = zeros(nn,1);
  lauf = 0;
  for i=1:n
    m = length(x{i});
    X(lauf+1:lauf+m) = x{i};
    lauf = lauf+m;
  end
  assert(lauf==nn);



