function [ObservGram,PCA,balred] = reduc3(A,B,SqrtCovmat,Sigma,Hbasic,nGramian)

  [U,S] = svd(SqrtCovmat,'econ');
  s = (diag(S)').^2;
  PCA = gramian(U,s);
  
  ObservGram = gramian(A',full(Hbasic),nGramian(1),1e-14);
  % CEA = gramianspace(A',full(Hbasic),nGramian(2));

  R = SqrtCovmat;
  Q = sqrt(ObservGram);
  [U,S,V] = svd(R'*Q,'econ');
  s = diag(S)';
  s2 = max(s,s(1)*1e-8);
  Ht = Q*(V.*repmat(1./sqrt(s2),size(V,1),1));
  Htinv = R*(U.*repmat(1./sqrt(s2),size(U,1),1));

  balred = gramian(Ht,s);


