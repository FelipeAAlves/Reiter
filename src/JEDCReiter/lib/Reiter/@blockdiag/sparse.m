function M = sparse(B)
  nn = size(B);
  nr = 0;
  nc = 0;
  iR=[]; iC=[];VV=[];
  
  for i=1:length(B.data)
    [nir,nic] = size(B.data{i});
    [ir,ic,vv] = find(sparse(B.data{i}));
    iR = [iR;ir+nr];
    iC = [iC;ic+nc];
    VV = [VV;vv];
    nr = nr + nir;
    nc = nc + nic;
  end
  assert(nr==nn(1));
  assert(nc==nn(2));
  M = sparse(iR,iC,VV,nr,nc);
    