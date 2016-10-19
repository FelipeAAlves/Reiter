function r = netintr(K,z);
    global MP;

    % adjust for differences when L ≠ 1
    KLratio = K / MP.aggEmployment;

    r = MP.A * z * MP.alpha * (KLratio) .^ (MP.alpha-1) - MP.delta;
