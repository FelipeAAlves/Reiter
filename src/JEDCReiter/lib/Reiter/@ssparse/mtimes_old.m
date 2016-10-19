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
    [irx,icx,vx] = find(x);
    if(~isa(y,'ssparse'))
      if(0)  % output dense:
	z = zeros(nx(1),ny(2));
	for i=1:length(irx)
	  for j=1:ny(2)
	    z(irx(i),j) = z(irx(i),j) + vx(i)*y(icx(i),j);
	  end
	end
      else  % output sparse:
	iR = zeros(length(irx)*ny(2),1);	iC = iR; V = iR;
	l = 1;
	for i=1:length(irx)
	  for j=1:ny(2)
	    iR(l) = irx(i);
	    iC(l) = j;
	    V(l) = vx(i)*y(icx(i),j);
	    l = l+1;
	  end
	end
	z = ssparse(iR,iC,V,nx(1),ny(2));
      end
    else
      nnx = nnz(x);
      nny = nnz(y);
      [iry,icy,vy] = find(y);
      M = [iry icy vy];
      M = sortrows(M,[1 2]);
      iry = M(:,1);
      icy = M(:,2);
      vy = M(:,3);

      if(0)
	II = zeros(1000,3);
	n = 0;
	ix1 = 1;
	iy1 = 1;
	while(ix1<=nnx & iy1<=nny)
	  jc = icx(ix1);
	  ix2 = nnx+1;
	  for l=ix1+1:nnx
	    if(icx(l)>jc)
	      ix2 = l;
	      break;
	    end
	  end
	  iix = ix1:(ix2-1);

	  if(iry(iy1)<jc)
	    for l=iy1+1:nny
	      if(iry(l)>=jc)
		iy1 = l;
		break;
	      end
	    end
	  end
	  iy2 = iy1;
	  if(iry(iy1)==jc)
	    iy2 = nny+1;
	    for l=iy1+1:nny
	      if(iry(l)>jc)
		iy2=l;
		break
	      end
	    end
	    iiy = iy1:(iy2-1);
	    n0 = length(iix)*length(iiy);
	    if(n+n0>size(II,1))
	      II = [II;zeros(n+n0,3)];
	    end
	    v = kron(vx(iix),vy(iiy));
	    ir = kron(irx(iix),ones(length(iiy),1));
	    ic = kron(ones(length(iix),1),icy(iiy));
	    II(n+1:n+n0,:) = [ir ic v];
	    n = n+n0;
	  end
	  ix1 = ix2;
	  iy1 = iy2;
	end
	z = ssparse(II(1:n,1),II(1:n,2),II(1:n,3),nx(1),ny(2));
      else
	[ir,ic,v] = helpmult(irx,icx,vx,iry,icy,vy);
	z = ssparse(ir,ic,v,nx(1),ny(2));
      end
    end
  else  %x not sparse, so y must be sparse
    if(0)  %output dense
      z = zeros(nx(1),ny(2));
      [iry,icy] = ind2sub2(y.n,y.i);
      for i=1:length(iry)
	for j=1:nx(1)
	  z(j,icy(i)) = z(j,icy(i)) + y.v(i)*x(j,iry(i));
	end
      end
    else
      [iry,icy] = ind2sub2(y.n,y.i);
      iR = zeros(length(iry)*nx(1),1);	iC = iR; V = iR;
      l = 1;
      for i=1:length(iry)
	for j=1:nx(1)
	  iR(l) = j;
	  iC(l) = icy(i);
	  V(l) = y.v(i)*x(j,iry(i));
	  l = l+1;
	end
      end
      z = ssparse(iR,iC,V,nx(1),ny(2));
    end
  end
      