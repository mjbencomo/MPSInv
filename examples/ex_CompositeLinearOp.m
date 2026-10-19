%% Compose interpolation and multi-channel convolution
% Requires LinearOp, ConvMultiOp, and the corrected TimeInterpOp on the path.
clear
K = 101;
tConv = linspace(0,1,K).';
dt = tConv(2)-tConv(1);
lags = (-(K-1):(K-1)).'*dt;
M = 2;
N = 3;
kernels = zeros(2*K-1,N,M);
for n = 1:N
    for m = 1:M
        kernels(:,n,m) = (lags >= 0).*exp(-(n+m)*max(lags,0));
    end
end
G = ConvMultiOp(kernels,Scale=dt);
tSource = linspace(0,1,31).';
tReceiver = linspace(0,1,46).';
Pin = TimeInterpOp(tSource,tConv,NumChannels=G.N_input);
Pout = TimeInterpOp(tConv,tReceiver,NumChannels=G.N_output);

% Mathematical product order: C = Pout o G o Pin.
C = CompositeLinearOp({Pout,G,Pin});
X = [exp(-((tSource-0.25)/0.08).^2), ...
     exp(-((tSource-0.55)/0.10).^2)];
x = X(:);
y = C.fwd(x);
Y = reshape(y,numel(tReceiver),N);
figure
plot(tReceiver,Y,'o-')
xlabel('Time')
ylabel('Receiver data')
grid on

z = sin((1:C.dimRange).');
lhs = dot(C.fwd(x),z);
rhs = dot(x,C.adj(z));
fprintf('Relative adjoint error: %.3e\n', ...
    abs(lhs-rhs)/max([1,abs(lhs),abs(rhs)]));

% Assemble only for small examples or verification.
A = C.matrixRep();
fprintf('Forward/matrix discrepancy: %.3e\n',norm(y-A*x));
