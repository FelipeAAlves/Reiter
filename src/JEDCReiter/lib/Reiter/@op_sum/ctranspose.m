function ST = ctranspose(S)
  ST.a = ctranspose(S.a);
  ST.b = ctranspose(S.b);
  ST = class(ST,'op_sum');
