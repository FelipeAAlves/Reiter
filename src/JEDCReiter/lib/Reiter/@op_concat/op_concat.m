function S = op_concat(a,b)
  S.a = a;
  S.b = b;
  assert(size(a,2)==size(b,1));
  S = class(S,'op_concat');
