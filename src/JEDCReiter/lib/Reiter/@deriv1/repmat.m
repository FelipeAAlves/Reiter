function r = simplerepmat(m,nr,nc)

  n = size(m);
  r.v = repmat(m.v,nr,nc);
  [ii,ideriv,v] = find(m.d);
  [ir,ic] = ind2sub2(n,ii);
  
  % replicate nr times vertically:
  iR = repmat(ir,nr,1) + kron(n(1)*(0:nr-1)',ones(length(ir),1));
  iC = repmat(ic,nr,1);


  % replicate nc times vertically:
  iCC = repmat(iC,nc,1) + kron(n(2)*(0:nc-1)',ones(length(iC),1));
  iRR = repmat(iR,nc,1);

  iDeriv = repmat(ideriv,nr*nc,1);
  V = repmat(v,nr*nc,1);

  r.d = ssparse(sub2ind2([n(1)*nr n(2)*nc],iRR,iCC),iDeriv,V,n(1)*nr*n(2)*nc,nindep(m));
  if(~dissparse(m))
    r.d = full(r.d);
  end
  r = deriv1(r);

  %c = cell(nr,1);
  %for i=1:nr
  %  c{i} = m;
  %end
  %m2 = vertcat(c{:});
  %c = cell(nc,1);
  %for i=1:nc
  %  c{i} = m2;
  %end
  %r = horzcat(c{:});

  
  