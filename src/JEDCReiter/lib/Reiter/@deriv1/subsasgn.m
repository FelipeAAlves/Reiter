function a=subsasgn(a,s,x)

  %if (isa(I1,'deriv1') | isa(I2,'deriv1'))
    %error('index must be integer, not deriv1');
    %end

  na = size(a);
  if length(s.subs) == 1
    I = handle_index(s.subs{1},prod(na));;
    if ~isa(x,'deriv1')
      a.v(I)=x; 
      a.d(I,:) = 0;
    else
      a.v(I)=x.v; 
      a.d(I,:) = x.d;
    end
  elseif length(s.subs)==2
    I1 = handle_index(s.subs{1},na(1));;
    I2 = handle_index(s.subs{2},na(2));;
    ii = gridmake(I1,I2);
    I = index(na,ii);
    if ~isa(x,'deriv1')
      a.v(I1,I2)=x; 
      a.d(I,:) = 0;
    else
      a.v(I1,I2)=x.v; 
      a.d(I,:)=x.d; 
    end
  else
	  n = prod(size(a));
	  indx = reshape(1:n,size(a));
	  lini = subsref(indx,s);
	  s2.type = '()';
	  s2.subs = {lini};
	  %ci = cell(1);
	  %ci{1} = lini;
	  a = subsasgn(a,s2,x);
    % error('wrong index in assignment to deriv1');
  end
      
function ii = handle_index(subs,nmax)
  if(ischar(subs))
    assert(subs(1)==':');
    ii = (1:nmax)';
  else
    ii = subs(:);
  end