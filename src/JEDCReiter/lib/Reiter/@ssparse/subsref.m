function xout=subsref(x,s)
  nx = size(x);
  switch s.type
    case '()'
      if(length(s.subs)==1)
	s1 = s.subs{1};
	if(ischar(s1))
	  assert(s1(1)==':');
	  xout = x;
	  xout.n = [prod(x.n) 1];
	else
	  if(islogical(s1))
	    s1 = find(s1);
	  end
	  Ranges = index2ranges(s1);
	  [iRC,V,nr] = getrange(x,Ranges(:,1),Ranges(:,2),Ranges(:,3));
	  xout = ssparse(iRC,[],V,size(s1,1),size(s1,2));
	end
      elseif(length(s.subs)==2)
	s1 = s.subs{1};
	if(islogical(s1))
	  s1 = find(s1);
	end
	s2 = s.subs{2};
	if(islogical(s2))
	  s2 = find(s2);
	end
	if(ischar(s2))
	  assert(s2(1)==':');
	  assert(isvector(s1));
	  Ranges = index2ranges(s1);
	  [iR,iC,V,nr] = getrowrange(x,Ranges(:,1),Ranges(:,2),Ranges(:,3));
	  xout = ssparse(iR,iC,V,nr,x.n(2));
	elseif(ischar(s1))
	  error('not yet implemented');
	else
	  assert(isvector(s1));
	  assert(isvector(s2));
	  ii = gridmake(s1(:),s2(:));
	  Ranges = index2ranges(sub2ind2(nx,ii(:,1),ii(:,2)));
	  [iRC,V,nr] = getrange(x,Ranges(:,1),Ranges(:,2),Ranges(:,3));
	  xout = ssparse(iRC,[],V,length(s1),length(s2));	  
	end
      else
	error('subscript to sparse matrix must be 2-dim');
    end
  otherwise
    error('wrong index type_ for_ deriv1');
  end

