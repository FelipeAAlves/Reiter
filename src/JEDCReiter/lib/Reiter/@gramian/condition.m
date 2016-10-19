function S2 = condition(S,eps)
  S2 = S;
  S2.s = max(S.s,S.s(1)*eps);
