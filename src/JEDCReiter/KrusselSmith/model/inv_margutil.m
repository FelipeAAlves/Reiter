
function c = inv_margutil(mu)

    global MP;
    assert(all(mu>0),'mu not positive');
    c = mu.^(-1/MP.gam);
