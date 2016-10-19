function xOut=power(x1,x2)
  if(~isa(x2,'deriv1') & isscalar(x2))
    if(x2==2)
      xOut = power2(x1);
      return;
    elseif(x2==3)
      xOut = power3(x1);
      return;
    end
  end
  xOut = power_gen(x1,x2);
