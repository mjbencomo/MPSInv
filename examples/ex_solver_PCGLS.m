%% PCGLS with a column-scaling preconditioner
n = 30;
d = logspace(0,-3,n).';
D = diag(d);
B = [D;0.5*D];
A = MatrixLinearOp(B);
xTrue = sin((1:n).');
b = A.fwd(xTrue);

% Jacobi preconditioner for the normal equations: M = diag(B'*B).
% Here the columns are orthogonal, so this equals B'*B exactly.
% General problems will not normally converge in a single PCGLS step.
M = MatrixInvLinearOp(diag(sum(abs(B).^2,1)));
[x,info] = solver_PCGLS(A,b,M, ...
    TrueSolution=xTrue,RelativeTolerance=1e-12, ...
    AbsoluteTolerance=0,MaxIterations=100);

fprintf('PCGLS: %s\n',info.message);
fprintf('Iterations: %d; true error: %.3e\n', ...
    info.iterations,info.trueErrorNorm);

figure;
tiledlayout(1,2);
nexttile;
plot(1:n,xTrue,'k-',1:n,x,'ro');
xlabel('Component'); ylabel('Solution');
legend('True','PCGLS','Location','best'); grid on;
nexttile;
semilogy(0:info.iterations,info.trueErrorHistory,'o-');
xlabel('Iteration'); ylabel('True error norm'); grid on;
