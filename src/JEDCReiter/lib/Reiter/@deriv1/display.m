function []=display(s)
  
  value = s.v;
  n = size(s.v);
  disp('Value:');
  disp( value( 1:min(5,n(1)) ) );
  if n(1)>5
      disp('...')
  end
  
  if(n(2)==1)
    deriv = s.d
  else
    disp('Deriv:');
    if( dissparse(s) )
      [ii,ider,vd] = find(s.d);
      [ir,ic] = ind2sub2(n,ii);
      for i=1:length(ii)
	    disp(sprintf('Row %d Col %d Var %d: %g',ir(i),ic(i),ider(i),vd(i)));
      end
    else
      for i=1:n(1)
	for j=1:n(2)
	  disp(sprintf('Row %d Col %d',i,j));
	  disp(s.d(index(n,[i j]),:));
	end
      end
    end
  end

