function ST = ctranspose(S)
  ST.a = ctranspose(S.b);
  ST.b = ctranspose(S.a);
  ST = class(ST,'op_concat');
