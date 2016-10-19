function xOut=mtimes(x1,x2)
  nSparse = 1;
  n1 = size(x1);
  n2 = size(x2);
  if(prod(n1)==1 | prod(n2)==1)
    xOut = x1.*x2;
    return;
  end
  x1v = getval(x1);
  x2v = getval(x2);

  xOut.v = x1v*x2v;
  nout = size(xOut.v);
  if(~isa(x2,'deriv1'))
    D1 = sparse(0);
  else
    if(n2(2)==1)
      D1 = sparse(x1v*x2.d);
    else
      if(isssparse(x2.d))
        x2cols = getallcols(x2.d);
        D1r = []; D1c = []; D1v = [];
        np = size(x2.d,2);
        assert(length(x2cols)==np);
        for i=1:np
          [ir,ic] = ind2sub2(n2,x2cols{i}.ir);
          tmp = x1v*sparse(ir,ic,x2cols{i}.v,n2(1),n2(2));
          [ir,ic,v] = find(tmp);
          D1r = [D1r;colvec(sub2ind2(nout,ir,ic))]; D1c = [D1c;i*ones(length(ir),1)]; D1v = [D1v;colvec(v)];
        end
	    clear x2cols;
	    D1 = ssparse(D1r,D1c,D1v,prod(nout),np);
	    clear D1r D1c D1v
      else
	    assert(0);
        ir = []; ic = []; v = [];
        for j=1:n2(2)
          [i0,j0,v0] = find(x1v*reshape(x2.d(:,j),n2));
          ir = [ir;sub2ind2(nout,i0,j0)];
          ic = [ic;j*ones(length(i0),1)];
          v = [v;v0];
        end
        D1 = sparse(ir,ic,v,nout(1),nout(2));
      end
    end
  end
  
  if(~isa(x1,'deriv1'))
    D2 = sparse(0);
  else
    if(n1(1)==1)
      D2 = x2v'*sparse(x1.d);
    else
      if(isssparse(x1.d))
	x1cols = getallcols(x1.d);
	D2r = []; D2c = []; D2v = [];
	np = size(x1.d,2);
	assert(length(x1cols)==np);
	for i=1:np
	  [ir,ic] = ind2sub2(n1,x1cols{i}.ir);
	  tmp = sparse(ir,ic,x1cols{i}.v,n1(1),n1(2))*x2v;
	  [ir,ic,v] = find(tmp);
	  D2r = [D2r;sub2ind2(nout,ir,ic)]; D2c = [D2c;i*ones(length(ir),1)]; D2v = [D2v;v];
	end
	clear x1cols;
	D2 = ssparse(D2r,D2c,D2v,prod(nout),np);
	clear D2r D2c D2v;
      else
	ir = []; ic = []; v = [];
	for j=1:n1(2)
	  [i0,j0,v0] = find(reshape(x1.d(:,j),n1)*x2v);
	  ir = [ir;sub2ind2(nout,i0,j0)];
	  ic = [ic;j*ones(length(i0),1)];
	  v = [v;v0];
	end
	D2 = sparse(ir,ic,v,nout(1),nout(2));
      end
    end
  end

  xOut.d = D1+D2;
  if(size(xOut.d,2)>=nSparse)
    xOut.d = ssparse(xOut.d);
  end
  xOut=deriv1(xOut);


function Z = kronmult(M,T)
  nM = size(M);
  nT = size(T);
  k = nT(1)/nM(2);
  if(isa(T,'ssparse') | issparse(T))
    [i,j,v] = find(T);
    im = floor((i-1)/nM(2)) + 1;
    Ts = cell(k,1);
    for l=1:length(im)
      Ts{im(l)} = [Ts{im(l)} l];
    end
    I = [];    J = [];    V = [];
    for l=1:k
      indx = Ts{l};
      if(~isempty(indx))
	[i2,j2,v2] = find(M*ssparse(i(indx)-(l-1)*nM(2),j(indx),v(indx),nM(2),nT(2)));
	I = [I;i2+(l-1)*nM(1)];
	J = [J;j2];
	V = [V;v2];
      end
    end
    Z = ssparse(I,J,V,k*nM(1),nT(2));
  else
    Z = zeros(nM(1)*k,nT(2));
    for i=1:k
      Z((i-1)*nM(1)+1:i*nM(1),:) = M*T((i-1)*nM(2)+1:i*nM(2),:);
    end
  end