function xOut = sparse(iR,iC,xIn,nr,nc)
  if(nargin<4)
    nr = max(iR);
  end
  if(nargin<5)
    nc = max(iC);
  end
  xOut.v = sparse(iR,iC,xIn.v,nr,nc);
  np = nindep(xIn);


  [ird,icd,vd] = find(xIn.d);
  iRd = iR(ird);
  iCd = iC(ird);
  indx = sub2ind2([nr nc],iRd,iCd);
  xOut.d = ssparse(indx,icd,vd,nr*nc,np);
  xOut=deriv1(xOut);
