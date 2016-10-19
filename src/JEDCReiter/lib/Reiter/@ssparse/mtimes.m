function z = mtimes(x,y)
  nx = size(x);
  ny = size(y);
  if(prod(nx)==1 | prod(ny)==1)
    z = times(x,y);
    return;
  end
  if(all([nx ny]<=1000000))  % use Matlab sparse multiplication
    if(isa(x,'ssparse'))
      x = sparse(x);
    end
    if(isa(y,'ssparse'))
      y = sparse(y);
    end
    z = ssparse(x*y);
    return;
  end

  if(issparse(x))
    x = ssparse(x);
  end
  if(issparse(y))
    y = ssparse(y);
    return;
  end
  if(isa(x,'ssparse'))
    if(~isa(y,'ssparse'))
      [irx,icx,vx] = find(x,1);  % sort the rows of x
      iR=[];iC=[];V=[];
      xrow=0;
      nrx = length(irx);
      for i=1:nrx
	irow = irx(i);
	xrow = xrow + vx(i)*y(icx(i),:);
	if(i==nrx | irow~=irx(i+1))
	  ic = find(xrow); 
	  if(~isempty(ic))
	    iR=[iR;repmat(irow,length(ic),1)]; iC=[iC;ic']; V=[V;xrow(ic)'];
	  end
	  xrow = 0;
	end
      end
      z = ssparse(iR,iC,V,nx(1),ny(2));
    else
      [irx,icx,vx] = find(x);
      [iry,icy,vy] = find(y,2);  % sort rows first, then columns
      [ir,ic,v] = helpmult(irx,icx,vx,iry,icy,vy);
      z = ssparse(ir,ic,v,nx(1),ny(2));
    end
  else  %x not sparse, so y must be sparse
    [iry,icy] = ind2sub2(y.n,y.i);
    iR=[];iC=[];V=[];
    xcol=0;
    nry = length(iry);
    for i=1:nry
      icol = icy(i);
      xcol = xcol + y.v(i)*x(:,iry(i));
      if(i==nry | icol<icy(i+1))
	ir = find(xcol);
	if(~isempty(ir))
	  iR=[iR;ir]; iC=[iC;repmat(icol,length(ir),1)]; V=[V;xcol(ir)];
	end
	xcol = 0;
      end
    end
    z = ssparse(iR,iC,V,nx(1),ny(2));
  end
      