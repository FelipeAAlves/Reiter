function [Y] = test_equ1(X)

    Y = initsize(X);
    Y(1) = 2*X(1) + 3*X(2);
    Y(2) = X(1) + X(2);
    
end

