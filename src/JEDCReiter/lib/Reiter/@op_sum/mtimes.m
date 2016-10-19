function z = mtimes(S,v)

if(isa(S,'op_sum'))
  z = S.a*v + S.b*v;
else
  z = S*v.a + S*v.b;
end


