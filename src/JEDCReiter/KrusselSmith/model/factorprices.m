function [R, wage] = factorprices(K, dz)

    global MP;

    %== Prices ==%
    R      = 1 + netintr(K,1+dz);
    wage   = wagefunc(K,1+dz);
