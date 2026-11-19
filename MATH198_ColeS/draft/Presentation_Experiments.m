%% Generating the Toeplitz Matrix
% The Setup
rcv_loc = [0.3, 0.7];
src_location0 = [0.5,0.4];
src_loc = [src_location0];
derivativeOrders = [0,0];
tc0 = 0.35; tc1 = 0.35; tc2 = 0.35;
w = {@(t) exp(-((t-tc0)/0.1).^2)};
Experiment1 = runGreenCGLS2D(rcv_loc,src_loc,derivativeOrders,w);

%The plot
A_G = Experiment1.greenOperator.matrixRep();
nSources = size(Experiment1.sourceLocations,1);
nSourceTimes = size(Experiment1.trueSourceSamples,1);
figure('Name','Green Operator Blocks')
for m = 1:nSources
    columns = (m-1)*nSourceTimes + (1:nSourceTimes);
    subplot(1,nSources,m)
    imagesc(A_G(:,columns))
    axis image
    colorbar
    xlabel('Source-time index')
    ylabel('Data index')
    title(sprintf('A_{%d}',m-1))
end

%% Monopole Source
%Creating receiver grid
x_rcv = linspace(0.25,0.75,10)';
y_rcv = linspace(0.25,0.65,8)';

rcv_loc = [
x_rcv, 0.25*ones(numel(x_rcv),1);
x_rcv, 0.65*ones(numel(x_rcv),1);
0.25*ones(numel(y_rcv)-2,1), y_rcv(2:end-1);
0.75*ones(numel(y_rcv)-2,1), y_rcv(2:end-1)];
src_location0 = [0.5,0.45];
src_loc = [src_location0];
derivativeOrders = [0,0];
tc0 = 0.35;
w = {@(t) exp(-((t-tc0)/0.1).^2)};

Experiment2 = runGreenCGLS2D(rcv_loc,src_loc,derivativeOrders,w,StopAtNoiseLevel=true);
plotGreenCGLS2D(Experiment2);

%% Dipole Source
%Creating receiver grid
x_rcv = linspace(0.25,0.75,10)';
y_rcv = linspace(0.25,0.65,8)';

rcv_loc = [
    x_rcv, 0.25*ones(numel(x_rcv),1);
    x_rcv, 0.65*ones(numel(x_rcv),1);
    0.25*ones(numel(y_rcv)-2,1), y_rcv(2:end-1);
    0.75*ones(numel(y_rcv)-2,1), y_rcv(2:end-1)];
src_location0 = [0.5,0.45];
src_loc = [src_location0; src_location0];
derivativeOrders = [1,0; 0,1];
tc1 = 0.35; tc2 = 0.35;
w = {@(t) exp(-((t-tc1)/0.1).^2), ...
    @(t) exp(-((t-tc2)/0.1).^2)};

Experiment3 = runGreenCGLS2D(rcv_loc,src_loc,derivativeOrders,w,StopAtNoiseLevel=true);
plotGreenCGLS2D(Experiment3);

%% Multipole Source

%Creating receiver grid
x_rcv = linspace(0.25,0.75,10)';
y_rcv = linspace(0.25,0.65,8)';

rcv_loc = [
    x_rcv, 0.25*ones(numel(x_rcv),1);
    x_rcv, 0.65*ones(numel(x_rcv),1);
    0.25*ones(numel(y_rcv)-2,1), y_rcv(2:end-1);
    0.75*ones(numel(y_rcv)-2,1), y_rcv(2:end-1)];
src_location0 = [0.5,0.45];
src_loc = [src_location0; src_location0; src_location0];
derivativeOrders = [0,0; 1,0; 0,1];
tc0 = 0.35; tc1 = 0.35; tc2 = 0.35;
w = {@(t) exp(-((t-tc0)/0.1).^2), ...
    @(t) exp(-((t-tc1)/0.1).^2), ...
    @(t) exp(-((t-tc2)/0.1).^2)};

Experiment4 = runGreenCGLS2D(rcv_loc,src_loc,derivativeOrders,w,StopAtNoiseLevel=true);
plotGreenCGLS2D(Experiment4);

%% Attempts at Preconditioning

A = Experiment4.forwardOperator;
B = A.matrixRep();
xTrue = Experiment4.trueSourceSamples(:);
b = Experiment4.noiseFreeData;

% Jacobi preconditioner for the normal equations: M = diag(B'*B).
% Here the columns are orthogonal, so this equals B'*B exactly.
% General problems will not normally converge in a single PCGLS step.
M = MatrixInvLinearOp(diag(sum(abs(B).^2,1)));
[x,info] = solver_PCGLS(A,b,M, ...
    TrueSolution=xTrue,RelativeTolerance=1e-12, ...
    AbsoluteTolerance=0,MaxIterations=1000);

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
