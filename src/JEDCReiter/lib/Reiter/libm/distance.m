function d = distance(m1,m2,iRel)
  if(nargin<3)
    scal = 1;
  elseif(iRel==1)
    scal = max(abs(m1(:)));
  elseif(iRel==2)
    scal = max(abs(m2(:)));
  else
    error('wrong iRel in distance');
  end
  d = max(abs(m1(:)-m2(:))) / (scal+1e-6);
