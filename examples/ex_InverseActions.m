%% Applying operators, their adjoints, and their inverses
A = [4+1i,1;2,3-1i];
op = MatrixInvLinearOp(A);
b = [1;2-1i];

x = op.inv(b);          % Solve A*x = b.
z = op.invAdj(b);       % Solve A'*z = b.
B = op.inverse();       % Operator wrapper for A^(-1); no solve yet.
H = op.adjoint();       % Operator wrapper for A'; no matrix assembly.

fprintf('Inverse residual: %.3e\n',norm(A*x-b));
fprintf('Inverse-adjoint residual: %.3e\n',norm(A'*z-b));
fprintf('Inverse-wrapper agreement: %.3e\n',norm(B.fwd(b)-x));
fprintf('Adjoint inverse agreement: %.3e\n',norm(H.inv(b)-z));

% Inverse operators can be composed using the existing CompositeLinearOp:
% identityOp = CompositeLinearOp({B,op});
% fprintf('Composition error: %.3e\n',norm(identityOp.fwd(b)-b));

% Explicit inverse matrix is computed only on request:
Ainv = B.matrixRep();
fprintf('Explicit inverse residual: %.3e\n',norm(A*Ainv-eye(2)));
