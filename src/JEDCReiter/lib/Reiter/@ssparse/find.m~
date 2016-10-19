% sortrows:
%    0 (or empty):  leave sorted columnwise
%    1 : sort rows, columns undefined
%    2 : sort rows first, columns second
function [i,j,v] = find(x,iSortrows)
  [i,j] = ind2sub2(x.n,x.i);
  i = i(:);
  j = j(:);
  if(nargout>2)
    v = x.v;
    v = v(:);
  end

  if(nargin>1 & iSortrows & ~isempty(i))
    if(iSortrows==1)
      [i,indx] = sort(i);
      j = j(indx);
      v = v(indx);
    elseif(iSortrows==2)
      nm = max(j)+1;
      % implementation of sort here: minimize memory requirements, don't generate too many temporaries:
      [i,indx] = sort(i+j/nm);
      i = floor(i);  
      j = j(indx);
      v = v(indx);
    else
      error('wrong iSortrows in find');
    end
  end
