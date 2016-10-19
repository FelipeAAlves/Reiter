function xOut = transpose(xIn)
  xOut.v = xIn.v';
  if(isvector(xIn.v) | nnz(xIn.d)==0)
    xOut.d = xIn.d;
  else
    if(dissparse(xIn))
      nn = size(xIn.d);
      [i,j,v] = find(xIn.d);
      n = size(xIn);
      indx = index(n,i);
      i2 = sub2ind2([n(2) n(1)],indx(:,2),indx(:,1));
      xOut.d = ssparse(i2,j,v,nn(1),nn(2));
    else
      nI = size(xIn.v);
      nO = [nI(2) nI(1)];
      nn = prod(nI);
      indxO = index(nO,(1:nn)');
      indxI = index(nI,[indxO(:,2) indxO(:,1)]);
      xOut.d = xIn.d(indxI,:);
    end
  end
  xOut=deriv1(xOut);
