%% Interpolate source samples, convolve, then sample receiver times
% Run after adding the project's src tree to the MATLAB path.
clear
K = 101;
t = linspace(0,1,K).';
dt = t(2)-t(1);
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

% Each interpolation uses the same times for every channel.
tSource = linspace(0,1,31).';
tReceiver = linspace(0,1,46).';
Pin = TimeInterpOp(tSource,t,NumChannels=G.N_input);
Pout = TimeInterpOp(t,tReceiver,NumChannels=G.N_output);

% Existing LinearOp has no composition overload. Function handles suffice
% until a general composition class is added. Adjoint order is reversed.
forward = @(x) Pout.fwd(G.fwd(Pin.fwd(x)));
adjoint = @(y) Pin.adj(G.adj(Pout.adj(y)));
X = [exp(-((tSource-0.25)/0.08).^2), ...
     exp(-((tSource-0.55)/0.10).^2)];
x = X(:);
y = forward(x);
Y = reshape(y,numel(tReceiver),N);
figure
plot(tReceiver,Y,'o-')
xlabel('Time')
ylabel('Receiver data')
grid on

% Verify the composed adjoint using deterministic data.
z = sin((1:Pout.dimRange).');
lhs = dot(forward(x),z);
rhs = dot(x,adjoint(z));
fprintf('Relative adjoint error: %.3e\n', ...
    abs(lhs-rhs)/max([1,abs(lhs),abs(rhs)]));
