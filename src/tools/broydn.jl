# computes Jacobian by forward differences
# Input:
# func: string, name of function
# x: point at which to take Jacobian
#    if function value at x is already known, then give x={starting point, function value}
# step: scalar, relative stepwidth;
# more arguments will be passed on to the function;
#
function jacob(func,x,step);
  f0 = func(x);
  n = length(x);
  m = length(f0);
  jac = zeros(m,n);
  x0 = copy(x);
  for i=1:n
    step2 = step*max(1,abs(x0[i]));
    y = copy(x0);
    y[i] = x0[i] + step2;
    jac[1:m,i] = (func(y) - f0)/step2;
  end;
  return jac;
end;
# compute Jacobian by central differences:
function jacobCentral(func,x,step)
  f0 = func(x);
  n = length(x);
  m = length(f0);
  jac = zeros(m,n);
  x0 = copy(x);
  for i=1:n
    step2 = step*max(1,abs(x0[i]));
    y = copy(x0);
    y[i] = x0[i] + step2;
    fUpp = func(y)
    y[i] = x0[i] - step2;
    fLow = func(y)
    jac[1:m,i] = (fUpp-fLow)/(2*step2);
  end;
  return jac;
end;



broydn_init_B = [];
broydn_out_B = [];
broydn_output_file = "C:\\Users\\falves\\Git\\Reiter\\src\\tools\\broydn_output.dat"; # []


# """
# Broydn's method to solve system of nonlinear equations
# ### INPUT:
# - `fname`:  name of vector function that should be set to zero
# - `xold` :   starting value for parameters
# - `opts` :   a vector of options; this can be empty, then just default values are taken;
#             otherwise:
#             opts[1]: tolerance level on function value (should be around 10^(-6))
#             opts[2]: 1 if Jacobian through automatic differentiation
#                (class deriv1) should be used, zero otherwise
#             opts[3]: 1 if output on iterations is desired
#             opts[4]: steplength for forward difference Jacobian, if used (should be around 10^(-5))
# - `dfname`
# ### OUTPUT:
# - `x`:      parameter vector that solves the system
# - `check`:  0 if ok, 1 if there is some problem (no solution found)
# The function fname should set the first element of its return vector to
# 1e100 in order to indicate overflow.
# """
function broydn(fname, xold, opts, dfname=x->x);
    global broydn_init_B;
    global broydn_output_file;
    # global Parameter broydn_init_B does not interfere with recursive calling of broydn,
    #  since only used in the beginning to initialize B, must be set anew before
    #  each call of broydn

    if (length(opts)>0)
        TOLF = opts[1];
    else  #default
        TOLF = 1e-6;
    end
    if (length(opts)>1)
        iAD = opts[2];
    else  #default
        iAD = 0;
    end
    if (length(opts)>2)
        iprint = opts[3];
    else  #default
        iprint = 0;
    end
    if (length(opts)>3)
        step_jacob = opts[4];
    else  #default
        step_jacob = 1e-6;
    end
    if (length(opts)>=5)
        MAXITS = opts[5];
    else  #default
        MAXITS = 200;
    end

    # initialize outputs:
    local B;
    B = broydn_init_B;
    broydn_init_B = [];  #must be set each time before calling broydn!

    xold = xold[:];
    x = xold;

    EPS  = 1.0e-7;
    TOLX = EPS;
    STPMX = 100.0;
    TOLMIN = 1.0e-6;
    n=size(x,1);

    (f,fvec) = fmin_br(x,fname);
    if (f > 1e4)
        FF = open("broydnfvec.dat","w");
        print(FF,fvec);
        close(FF);
    end
    fvcold = fvec
    fold   = f
    fx = fvec
    if abs( fvec[1] )>=1e100
        warn("overflow in function given to broydn at initial vector");
        return (x, -1)
    end

    test = maximum(abs(fvec));
    if (test<TOLF)  #changed from 0.01*TOLF to TOLF;
        return (x, 0)
    end
    #stpmax=STPMX*max(sqrt(dot(x,x)),n);

    stpmax=1000.;
    if (isempty(B))
        restrt=1;
    else
        restrt=0;
    end

    since_restrt = 0;
    alam = 1;
    for its=1:MAXITS
        since_restrt = since_restrt+1;
        txt = @sprintf("  its: %2d  f: %e    step size taken: %e\n",its, f, alam);
        if (isfile("broydn_stop"))
            return (x,4);
        end
        if (iprint>0)
            if (!isempty(broydn_output_file))
                FF = open(broydn_output_file,"a");
                (its==1) && print(FF,"ITERATIONS: \n");
                print(FF,txt);
                if (iprint>1)
                    print(txt)
                    println(FF,x);
                    println(FF,fvec);
                end
                close(FF)
            end
        end
        # determine whether restart (new computation of Jacobian):
        if (isfile("broydn_restrt.txt"))
            system("del broydn_restrt.txt");
            restrt=1;
        end

        if (restrt==1 || since_restrt>=3*n) #compute Jacobian at 3n iterations since last restart
            since_restrt = 0;
            if (iAD==1)
                (iprint>1) && println("now do Jacobian")

                B = dfname(x);
                # saveT("bf.bmt",B,fvec);stop()
                if (!all(isfinite(B)))
                    ii = find(!all(isfinite(sum(B,2)))); # columns not all finite
                    FF = open("bnotfinite.dat","w");
                    println(FF,ii);
                    println(FF,fvec[ii]);
                    println(FF,B[ii[1],:]);
                    close(FF);
                end
                # fvec = getval(fvec); # not needed??
            else
                B = jacob(fname,x,step_jacob);
                @assert(size(B,1)==size(B,2),@printf("size(B) = %d,%d\n",size(B,1),size(B,2)))

                #saveT("b.bmt",B);
                # println(cond(B)); #done()
            end

            global broydn_out_B;
            broydn_out_B = B;

        elseif (its>1)
            s = x - xold;
            skip = 1;
            w = (fvec-fvcold) - B*s;
            for i=1:n
                if abs(w[i]) >= EPS*(abs(fvec[i])+abs(fvcold[i]))
                    skip=0
                else
                    w[i]=0.0
                end
            end

            if skip==0
                B = B + w*s' / dot(s,s);
            end
        end
        if (length(B)<=20 && cond(B)>1e14)
            println(B)
            println(cond(B))
            error("condition B")
        end
        # maxr = maximum(abs(B),2);
        # p = - B \ fvec;
        # # p = - (B./(repmat(maxr,1,size(B,1)))) \ (fvec./maxr);
        Bscaled = B  + 1e-20*eye(size(B,1));
        maxr = 1.0 ./ maximum(abs(Bscaled),2);
        Bscaled = broadcast(*,Bscaled,maxr);
        maxc = 1.0 ./ maximum(abs(Bscaled),1);
        Bscaled = broadcast(*,Bscaled,maxc);
        if (iprint>1)
            condB = cond(Bscaled);
            txt = @sprintf("condition number of scaled Jacobian in broydn is %e",condB);
            println(txt);
            if (!isempty(broydn_output_file))
                FF = open(broydn_output_file,"a");
                print(FF,txt);
                close(FF)
            end
            condB>1e14 && error("condition B")
            saveT("bmat.bmt",B,Bscaled);
        end
        #error("done")
        p = - (Bscaled\(fvec.*maxr)) .* maxc';
        p = p[:]
        g = B'*fvec;
        xold = x;
        fvcold = fvec;
        fold=f;
        (x,f,fvec,check,alam) = lnsrch_br(xold,fold,g,p,stpmax,fname);
        fx = fvec;
        test = maximum(abs(fvec));
        if (test < TOLF)
            return (x, 0)
        end
        if check==1
            if (restrt>0)
                return (x, 1)
            else
                test=0.0;
                den=max(f,0.5*n);
                test = maximum(abs(g) .* maximum([abs(x');ones(1,n)])') / den;
                if (test < TOLMIN)
                    return (x, 0)
                else
                    restrt=1;
                end
            end
        else
            alam<1e-3 ? (restrt=1) : (restrt=0)
            test= maximum(abs(x-xold) ./ maximum([abs(x');ones(1,n)])');
            if (test < TOLX)
                return (x, 0)
            end
        end
    end

    (iprint>0) && println("MXITS exceeded in broydn")
    return (x, 3)
end

function broydnUnscaled(fname,xold,opts,dfname=x->x);
  global broydn_init_B;
  global broydn_output_file;
  # global Parameter broydn_init_B does not interfere with recursive calling of broydn,
  #  since only used in the beginning to initialize B, must be set anew before
  #  each call of broydn

  if(length(opts)>0)
    TOLF = opts[1];
  else  #default
    TOLF = 1e-6;
  end
  if(length(opts)>1)
    iAD = opts[2];
  else  #default
    iAD = 0;
  end
  if(length(opts)>2)
    iprint = opts[3];
  else  #default
    iprint = 0;
  end
  if(length(opts)>3)
    step_jacob = opts[4];
  else  #default
    step_jacob = 1e-6;
  end
  if(length(opts)>=5)
    MAXITS = opts[5];
  else  #default
    MAXITS = 200;
  end

  # initialize outputs:
  local B;
  B = broydn_init_B;
  broydn_init_B = [];  #must be set each time before calling broydn!

  xold = xold[:];
  x = xold;

  EPS = 1.0e-7;
  TOLX = EPS;
  STPMX = 100.0;
  TOLMIN = 1.0e-6;
  n=size(x,1);

  (f,fvec)=fmin_br(x,fname);
  if(f>1e4)
    FF = open("broydnfvec.dat","w");
    print(FF,fvec);
    close(FF);
  end
  fvcold = fvec;
  fold=f;
  fx = fvec;
  if abs(fvec[1])>=1e100
    warn("overflow in function given to broydn at initial vector");
    return (x, -1)
  end;
  test = maximum(abs(fvec));
  if (test<TOLF)  #changed from 0.01*TOLF to TOLF;
    return (x, 0)
  end;
  #stpmax=STPMX*max(sqrt(dot(x,x)),n);
  stpmax=1000.;
  if(isempty(B))
    restrt=1;
  else
    restrt=0;
  end
  since_restrt = 0;
  alam = 1;
  for its=1:MAXITS
    since_restrt = since_restrt+1;
    txt = @sprintf("its:  %d; f:  %e; step size taken: %e\n",its,getval(f),alam);
    if(isfile("broydn_stop"))
      return (x,4);
    end
    if(iprint>0)
      print(txt)
      if(!isempty(broydn_output_file))
        FF = open(broydn_output_file,"a");
        print(FF,txt);
        if(iprint>2)
          println(FF,x);
          println(FF,fvec);
        end
        close(FF)
      end;
    end
    # determine whether restart (new computation of Jacobian):
    if(isfile("broydn_restrt.txt"))
      system("del broydn_restrt.txt");
      restrt=1;
    end
    if(restrt==1 || since_restrt>=3*n) #compute Jacobian at 3n iterations since last restart
      since_restrt = 0;
      if(iAD==1)
        if(iprint>1)
          println("now do Jacobian")
        end
        @time B = dfname(x);
        # saveT("bf.bmt",B,fvec);stop()
        if(!all(isfinite(B)))
          ii = find(!all(isfinite(sum(B,2)))); # columns not all finite
          FF = open("bnotfinite.dat","w");
          println(FF,ii);
          println(FF,fvec[ii]);
          println(FF,B[ii[1],:]);
          close(FF);
        end
        # fvec = getval(fvec); # not needed??
      else
        B = jacob(fname,x,step_jacob);
        @assert(size(B,1)==size(B,2),@printf("size(B) = %d,%d\n",size(B,1),size(B,2)))

        #saveT("b.bmt",B);
        # println(cond(B)); #done()
      end;
      global broydn_out_B = B;
    elseif(its>1)
      s = x - xold;
      skip=1;
      w = (fvec-fvcold) - B*s;
      for i=1:n
        if abs(w[i]) >= EPS*(abs(fvec[i])+abs(fvcold[i]))
          skip=0;
        else
          w[i]=0.0;
        end;
      end;
      if skip==0
        B = B + w*s' / dot(s,s);
      end;
    end;
    if(length(B)<=20)
      if(cond(B)>1e14)
        error("condition B")
      end
    end
    p = - B\fvec;
    p = p[:]
    g = B'*fvec;
    xold = x;
    fvcold = fvec;
    fold=f;
    (x,f,fvec,check,alam) = lnsrch_br(xold,fold,g,p,stpmax,fname);
    fx = fvec;
    test = maximum(abs(fvec));
    if (test < TOLF)
      return (x, 0)
    end;
    if check==1
      if (restrt>0)
        return (x, 1)
      else
        test=0.0;
        den=max(f,0.5*n);
        test = maximum(abs(g) .* maximum([abs(x');ones(1,n)])') / den;
        if (test < TOLMIN)
          return (x, 0)
        else
          restrt=1;
        end;
      end;
    else
      restrt=0;
      if(alam<1e-3)
        restrt=1;
      end
      test= maximum(abs(x-xold) ./ maximum([abs(x');ones(1,n)])');
      if (test < TOLX)
        return (x, 0)
      end;
    end;
  end;
  # warning("MXITS exceeded in broydn");
  return (x, 3)
end



function lnsrch_br(xold,fold::Float64,g,p,stpmax::Float64,funcname)
  n = size(xold,1);
  ALF = 1.0e-4;
  TOLX = 1.0e-10;
  check=0;
  sum = dot(p,p);
  sum=sqrt(sum);
  if (sum > stpmax)
    p = p * stpmax/sum;
  end;
  slope = dot(g,p);
  test=maximum(abs(p) ./ maximum([abs(xold');ones(1,n)],1)');
  # min() introduced by MR 5.4.2006:
  # otherwise, alamin can be >1 !!
  alamin=min(0.1,TOLX/test[1]);
  alam=1.0;
  f2 = NaN;
  fold2 = NaN;
  alam2 = NaN;
  f = NaN;
  fvec = [];
  for its=1:100
    # print([xold p])
    x=xold + alam*p;
    try
      (f,fvec)=fmin_br(x,funcname);
    catch
      f = 1e100;  #something outrageous;
    end;
    if (alam < alamin)
      x = xold;
      check=1;
      return (x,f,fvec,check,alam)
    else
      if(f[1]>=1e100)
        alam = 0.1*alam;
        continue;
      end
      if f <= fold+ALF*alam*slope
        return (x,f,fvec,check,alam)
      else
        if(~isfinite(f2)) # not the first iteration!
          tmplam = -slope/(2.0*(f-fold-slope));
        else
          rhs1 = f-fold-alam*slope;
          rhs2=f2-fold2-alam2*slope;
          a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
          b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
          if (a == 0.0)
            tmplam = -slope/(2.0*b);
          else
            disc=b*b-3.0*a*slope;
            if (disc<0.0)
              error("Roundoff problem in lnsrch");
            else
              tmplam=(-b+sqrt(disc))/(3.0*a);
            end;
          end;
          if (tmplam>0.5*alam)
            tmplam=0.5*alam;
          end;
        end;
      end;
    end;
    alam2=alam;
    f2 = f;
    fold2=fold;
    alam=max(tmplam,0.1*alam);
  end;
end

function fmin_br(x,funcname);

  fvec = funcname(x);
  f = 0.5*dot(fvec,fvec);
  return (f,fvec)
end

# solve min_x (f(x) + B*dx)'(f(x) + B*dx)
# Normal equations:
# B'*B dx + B'f(x) = 0
# use Tychonov regularization
# then f(x) + B*dx = (-B*inv(B'*B)B' + I)*f(x)
function stepwiseOLS(fname,xold,opts,dfname=x->x);
  # global Parameter broydn_init_B does not interfere with recursive calling of broydn,
  #  since only used in the beginning to initialize B, must be set anew before
  #  each call of broydn

  if(length(opts)>0)
    TOLF = opts[1];
  else  #default
    TOLF = 1e-6;
  end
  if(length(opts)>1)
    iAD = opts[2];
  else  #default
    iAD = 0;
  end
  if(length(opts)>2)
    iprint = opts[3];
  else  #default
    iprint = 0;
  end
  if(length(opts)>3)
    step_jacob = opts[4];
  else  #default
    step_jacob = 1e-6;
  end
  if(length(opts)>=5)
    MAXITS = opts[5];
  else  #default
    MAXITS = 20;
  end

  # initialize outputs:;rid
  check = 0;
  local B;

  xold = xold[:];
  x = xold;

  EPS = 1.0e-7;
  TOLX = EPS;
  STPMX = 100.0;
  TOLMIN = 1.0e-6;
  n=size(x,1);

  (f,fvec)=fmin_br(x,fname);
  fStart = f;
  fvcold = fvec;
  fold=f;
  fx = fvec;
  if abs(fvec[1])>=1e100
    warning("overflow in function given to stepwiseOLS at initial vector");
    check = -1;
    return (x, check)
  end;
  if (maxabs(fvec)<TOLF)
    return (x, 0)
  end;
  stpmax=1000.;
  restrt=1;
  alam = 1;
  for its=1:MAXITS
    txt = @sprintf("its:  %d; f:  %e; step size taken: %e\n",its,getval(f),alam);
    if(iprint>0)
      print(txt)
    end
    # determine whether restart (new computation of Jacobian):
    if(restrt==1) #compute Jacobian at 3n iterations since last restart
      if(iAD==1)
        B = dfname(x);
        if(!all(isfinite(B)))
          ii = find(!all(isfinite(sum(B,2)))); # columns not all finite
          FF = open("bnotfinite.dat","w");
          println(FF,ii);
          println(FF,fvec[ii]);
          println(FF,B[ii,:]);
          close(FF);
        end
        # fvec = getval(fvec); # not needed??
      else
        println("start Jacob")
        B = jacob(fname,x,step_jacob);
        #println("end Jacob")
      end;
    end;
    p = -ols_ridge(B,fvec,1e-6);
    fvecPredict = fvec + B*p
    fPredict = 0.5*dot(fvecPredict,fvecPredict);
    @printf("predicted sumsqr = %e\n",fPredict)
    xold = x;
    fvcold = fvec;
    fold=f;
    alam = 1;
    #println(sqrt(p'*p))
    while(alam>1e-5)
      x = xold + alam*p;
      (f,fvec)=fmin_br(x,fname);
      #println(f)
      if(f<fold)
        restrt=0;
        break;
      end
      alam *= 0.75;
    end
    if(f<fStart*1e-3 || maxabs(fvec)<1e-6)
      @printf("alam and f at exit: %f, %e\n",alam,f);
      return (x,0);
    end
    if(f<fStart+0.9*(fPredict-fStart))
      @printf("alam and f at exit: %f, %e\n",alam,f);
      return (x,0);
    end
    if(alam<1e-5)
      if(restrt==1)
        return (xold,1);
      else
        x = xold;
        f = fold;
        fvec = fvcold;
        restrt=1;
      end
    end
    fx = fvec;
  end;
  return (x,2); # maximum iteration exceeded
end

_jacobSimpleNewton = Array{Float64,2}()
_jacobSimpleNewtonInv = lufact(eye(2))
function simpleNewton(func,dfunc,x0,maxits,ftol,doPrint=0)
  x = copy(x0);
  f = func(x);
  if(maxabs(f).<ftol)
    return (x,0);
  end

  local B,Binv;
  restart = false
  if(isempty(_jacobSimpleNewton) || size(_jacobSimpleNewton,1)!=length(f))
    restart = true;
  else
    # B = copy(_jacobSimpleNewton)
    Binv = _jacobSimpleNewtonInv
  end
  alam = NaN
  for its = 1:maxits
    if(restart)
      (f,B) = dfunc(x);
      Binv = lufact(B);
      global _jacobSimpleNewton,_jacobSimpleNewtonInv;
      _jacobSimpleNewton = copy(B) # no copy for LU-factorization!! fix it?
      _jacobSimpleNewtonInv = Base.LinAlg.LU(copy(Binv.factors),copy(Binv.ipiv),Binv.info)
    else
      f = func(x);
    end
    fsqr = dot(f,f);
    fcrit = maxabs(f)
    if(doPrint>0)
      @printf("%e; %e, %f\n",fsqr,fcrit,alam);
    end

    if(maxabs(f).<ftol)
      return (x,0);
    end
    p = -(Binv\f);
    alam = 1;
    # if(doPrint>=3)
    #   saveT(@sprintf("res0.bmt"),x,f,B)
    # end
    for i=1:5
      xtry = x + alam*p;
      # if(doPrint>=2)
      #   println(xtry[1:12])
      # end
      try
        ftry = func(xtry);
        # saveT(@sprintf("res%d.bmt",i),xtry,ftry)
        f2 = dot(ftry,ftry)
        if(f2<fsqr)
          x = xtry;
          break;
        end
        alam = alam/10;
        restart = true;
      catch
        alam = alam/10;
        restart = true;
      end
    end
  end
  return (x,1); # signal failure
end

function wrapBroydn(f::Function,fTol,iPrint,xStart::Array{Float64},args...)
  nX = length(xStart);
  return broydn(x->f(x,args...),xStart,[fTol,1,iPrint],
  x -> getjac(f(indep(Deriv1a{nX},x),args...)));
end


type Fsolvepars
  pNewton::Array{Float64,1}
  pc::Array{Float64,1}
  rho::Float64
  tau::Float64
  Delta::Float64
end
sumsqr(x) = dot(x,x)
mulMM(a,b) = a*b
mulMTM(a,b) = a'*b
function mycopy!(x,y)
  for i=1:length(x)
    x[i] = y[i];
  end
end
function mycopy!(x,y,n::Int)
  for i=1:n
    x[i] = y[i];
  end
end
function mycopy!(x,offs::Int,n::Int,y)
  for i=1:n
    x[offs+i] = y[i];
  end
end

function fsolve_trust(f,df,x0::Array{Float64,1}, DeltaBar::Float64, TOLF::Float64=1e-6,  iprint::Int=0, maxiter = 200)

  fp = Fsolvepars(copy(x0),copy(x0),0.0,1.0,1.0);
  fp.Delta = DeltaBar/2.;
  M=0.0; Mtry=0.0;Mj=0.0;
  eta = 0.01;
  x = copy(x0)
  new_x = 1;
  local r,J
  for its=1:maxiter
    if(new_x>0)
      (r::Array{Float64,1},J) = df(x);
      if(!all(isfinite(r)))
        runtime_error("value not finite at iter %d in fsolve_trust",its);
      end
      Jscaled = copy(J)
      # Jscaled = J  + 1e-20*eye(size(J,1));
      # maxr = vec(1.0 ./ maximum(abs(Jscaled),2));
      # Jscaled = broadcast(*,Jscaled,maxr);
      # maxc = vec(1.0 ./ maximum(abs(Jscaled),1));
      # Jscaled = broadcast(*,Jscaled,maxc);
      M = sumsqr(r);
      if(iprint>0)
        @printf("fsolve,its:  %d; sumsqu:  %e; bounds: %e\n",its,M,fp.Delta);
      end
      if(maxabs(r)<TOLF)
        return (x,0);
      end
      # saveT("Jfsolve.bmt",Jscaled)
      _tmp = - (Jscaled\r)
      # _tmp = - vec((Jscaled\(r.*maxr)) .* maxc);
      mycopy!(fp.pNewton,_tmp)
    end
    at_bound = true;
    if(norm(fp.pNewton)<=fp.Delta) #take full Newton step:
      p = copy(fp.pNewton);
      at_bound = false;
    else # steepest descent direction:
      Jr = mulMTM(J,r);
      JJr = mulMM(J,Jr);
      rJJJJr = sumsqr(JJr);
      NJr = norm(Jr);
      fp.tau = NJr*NJr*NJr/(fp.Delta*rJJJJr);
      if(fp.tau>=1)
        # take steepest descent to boundary:
        p = -(fp.Delta/NJr) * Jr;
      else
        # compute Cauchy point:
        mycopy!(fp.pc,-(fp.tau*fp.Delta/NJr) * Jr);
        # from there, take Newton step to boundary:
        (fp.tau,flag::Int) = bisect(t-> norm(fp.pc+t*fp.pNewton)-fp.Delta,0.,1.,1e-6);
        assert(flag==0)
        p = fp.pc+fp.tau*fp.pNewton;
      end
    end
    if(its==1)
      pold = copy(p);
    end
    xtry = x + p;
    local rtry
    try
      rtry::Array{Float64,1} = f(xtry);
      Mtry = sumsqr(rtry);
      Mj = sumsqr(r + mulMM(J,p));
      if(Mtry>=M)
        fp.rho = -1;
      else
        fp.rho = (M-Mtry) / (M-Mj);
      end
    catch
      fp.rho = -1;
    end
    if(!isfinite(fp.rho))
      fp.rho = -1;
    end
    if(fp.rho<0.25)
      fp.Delta = 0.25*fp.Delta;
      if(fp.Delta<1e-12*DeltaBar)
        return (x,1);
      end
    elseif(fp.rho>0.75 && at_bound)
      fp.Delta = min(DeltaBar,2.0*fp.Delta);
    end
    if(fp.rho>eta || (fp.rho>-0.5 && DeltaBar<1e-7))
      new_x = 1;
      x = xtry;
      r = rtry;
      M = Mtry;
      pold = copy(p);
    else
      new_x = 0;
    end
  end
  # error("not converged in fsolve");
  return (x,1); # never get here
end
