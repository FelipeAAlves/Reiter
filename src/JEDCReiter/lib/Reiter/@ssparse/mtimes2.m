function z = mtimes(x,y)
  nx = size(x);
  ny = size(y);
  if(prod(nx)==1 | prod(ny)==1)
    z = times(x,y);
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
      nry = length(iry);
      for i=1:nry
	irow = iry(i);
	xrow = xrow + vx(i)*x(icy(i),:);
	if(i==nry | irow<iry(i+1))
	  ir = find(xrow);
	  if(~isempty(ir))
	    iR=[iR;repmat(irow,length(ir),1)]; iC=[iC;ir]; V=[V;xrow(ir)];
	  end
	  xrow = 0;
	end
      end
      z = ssparse(iR,iC,V,nx(1),ny(2));
    else
      [irx,icx,vx] = find(x);
      nnx = nnz(x);
      nny = nnz(y);
      [iry,icy,vy] = find(y);
      M = [iry icy vy];
      M = sortrows(M,[1 2]);
      iry = M(:,1);
      icy = M(:,2);
      vy = M(:,3);

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
      