function S = getcolrange(X,iRange,offset)
  if(nargin<3)
    offset = 0;
  end
  [ir,ic] = ind2sub2(X.n,X.i);
  ia = firstindex(ic,iRange(1));
  ib = lastindex(ic,iRange(2));
  ii = ia:ib;
  S.v = X.v(ii);
  S.n = [X.n(1) iRange(2)-iRange(1)+1];
  S.i = sub2ind2(S.n,ir(ii),ic(ii)-(iRange(1)+1+offset));
  S = class(S,'ssparse');
