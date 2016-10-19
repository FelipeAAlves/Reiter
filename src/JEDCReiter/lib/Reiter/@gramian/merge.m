function g = merge(g1,g2)
g.U = blockdiag(g1.U,g2.U);
g.s = [g1.s g2.s];

