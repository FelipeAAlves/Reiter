function D = invdistr(Pi)

    %== invariant distribution ==%
    D.v = invdistr(Pi.v);

    %== derivative ?? ==%
    n = length(Pi.v);
    p = nindep(Pi);
    for i=1:p
        rhs = reshape(Pi.d(:,i),n,n)*D.v;
    end

    D=deriv1(D);
