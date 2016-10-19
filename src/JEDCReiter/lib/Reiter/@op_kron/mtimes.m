function z = mtimes(M,v)

  [mr,mc] = size(M);
  [ar,ac] = size(M.a);
  [br,bc] = size(M.b);
  [vr,vc] = size(v);
  assert(mc==vr);

  zz = cell(ar,1);
  for j=1:ac
    offs = (j-1)*bc;
    mv = M.b*v(offs+1:offs+bc,:);
    for i=1:ar
      if(M.a(i,j)~=0)
	if(isempty(zz{i}))
	  zz{i} = M.a(i,j)*mv;
	else
	  zz{i} = zz{i} + M.a(i,j)*mv;
	end
      end
    end
  end
  z = vertcat(zz{:});




