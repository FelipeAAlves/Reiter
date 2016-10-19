function xout=subsref(x,s)
  nx = size(x);
  switch s.type
    case '()'
      % xout.v = x.v(s.subs{:});
      xout.v = subsref(x.v,s);
      Is = reshape(1:prod(size(x.v)),size(x.v));
      ii = subsref(Is,s);
      xout.d = x.d(ii(:),:);
      
        if (nnz(xout.d)==0)
            xout = xout.v;
        else
            xout = deriv1(xout);
        end
    case '.'
        switch s(1).subs
            case 'v'
                xout = x.v;
            case 'd'
                xout = x.d;
        end
  end


