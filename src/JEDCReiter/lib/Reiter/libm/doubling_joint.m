function [S1,S2,S12,cov1,cov12,maxerror] = doubling_joint(A1,B1,A2,B2,Sigma,nIter,h1,h2,H)

  n2 = length(A2);
  for i=1:n2
    S12{i} = B1*Sigma*B2{i}';
    S2{i} = B2{i}*Sigma*B2{i}';
  end
  S1 = B1*Sigma*B1';
  for iIter=1:nIter
    S1 = S1 + A1*(S1*A1');
    S2 = apply2cell(@plus,S2,apply2cell(@quadform,A2,S2));
    for i=1:n2
      S12{i} = S12{i} + A1*(S12{i}*A2{i}');
    end
    if(nargout>5)
      for i=1:n2
        maxerror(iIter,i) = norm(h1*A1 - h2{i}*A2{i}*H{i});
      end
    end
    if(iIter<nIter)
      A1 = A1*A1;
      A2 = apply2cell(@square,A2);
    end
  end
  
  if(nargout>3)
    cov1 = quadform(h1,S1);
  end
  
  if(nargout>4)
    for i=1:n2
      cov12{i} =  quadform(h1,S1) + quadform(h2{i},S2{i}) - 2*h1*S12{i}*h2{i}';
    end
  end

  
function x2 = square(x)
x2 = x*x;


