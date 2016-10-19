function S = getallcols(X)
  nc = size(X,2);
  [ir,ic] = ind2sub2(X.n,X.i);
  lastc = 0;
  for j=1:length(ic)
    if(ic(j)>lastc)
      for l=lastc+1:ic(j)
	iBeg(l)=j;
      end
      lastc = ic(j);
    end
  end
  for l=lastc+1:nc+1
    iBeg(l)=length(ic)+1;
  end
  for j=1:nc
    ii = iBeg(j):iBeg(j+1)-1;
    S{j} = struct('ir',ir(ii),'v',X.v(ii));
  end
