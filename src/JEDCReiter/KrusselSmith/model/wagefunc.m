
function wage = wagefunc(K,z);

    global MP;

    % adjust for differences when L ≠ 1
    KLratio = K / MP.aggEmployment;

    wage = MP.A * z * ( 1 - MP.alpha ) * (KLratio) .^ MP.alpha;
