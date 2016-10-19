function z = elelmult_eachcol(vec,mat)
  
  if(isscalar(vec))
    z = vec.*mat;
  elseif(size(mat,1)==1)
    z = vec*mat;
  else
    assert(isa(mat,'ssparse'));
    assert(size(vec,2)==1);
    n = length(vec);
    if(issparse(vec) | isa(vec,'ssparse'))
      [ir,ic,v] = find(vec);
      M = ssparse(ir,ir,v,n,n);
    else
      M  = ssparse(1:n,1:n,vec,n,n);
    end
    z = M*mat;
  end

  % DOES THE FOLLOWING:
  %z = zeros(size(mat));
  %for j=1:size(mat,2)
  %  z(:,j) = vec.*mat(:,j);
  %end