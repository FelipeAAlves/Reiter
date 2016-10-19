function [Y] = test_equ2(X)

    Y = initsize(X);
    Y(1) = exp( X(1) ) * X(2);
    Y(2) = log(X(1) + X(2));
    
end

