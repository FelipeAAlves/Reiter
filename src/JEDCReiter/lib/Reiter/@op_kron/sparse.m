function A = sparse(S)
  A = kron(sparse(S.a),sparse(S.b));
