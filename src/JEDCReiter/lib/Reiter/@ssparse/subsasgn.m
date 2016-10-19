function x=subsasgn(x,s,y)
  switch s.type
    case '()'
      if(length(s.subs)==1)
	s1 = s.subs{1};
	if(islogical(s1))
	  s1 = find(s1);
	end
	Ranges = index2ranges(s1);
	x = setrange(x,Ranges(:,1),Ranges(:,2),Ranges(:,3),y);
      elseif(length(s.subs)==2)
	s1 = s.subs{1};
	if(islogical(s1))
	  s1 = find(s1);
    end
    if(isempty(s1))
        return;
    end
	s2 = s.subs{2};
	if(ischar(s2))
	  assert(s2(1)==':');
	  assert(isvector(s1));
	  Ranges = index2ranges(s1);
	  x = setrowrange(x,Ranges(:,1),Ranges(:,2),Ranges(:,3),ssparse(y));
	else
	  error('not yet implemented');
	end
      else
	error('subscript to sparse matrix must be 2-dim');
    end
  otherwise
    error('wrong index type_ for_ deriv1');
  end

