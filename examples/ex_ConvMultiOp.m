%% Multi-channel Convolution as a Linear Operator
% This example constructs a ConvMultiOp with two input channels and three output channels. It illustrates the channel-block vector convention, the forward and adjoint actions, and the matrix representation.

clear
close all

%% Time grid and channel dimensions
% Each channel has K samples on the same uniform time grid. A kernel has 2K-1 lag samples, ordered from -(K-1)dt through (K-1)dt. The zero-lag value is stored at index K.

T = 1;
K = 101;
times = linspace(0,T,K).';
dt = times(2)-times(1);
lags = (-(K-1):(K-1)).'*dt;

M = 2;  % input channels
N = 3;  % output channels

%% Construct convolution kernels
% kernels(:,n,m) maps input channel m to output channel n. These illustrative kernels are causal, so their negative-lag samples vanish.

kernels = zeros(2*K-1,N,M);
causal = lags >= 0;

for n = 1:N
    for m = 1:M
        decayRate = 2+n+m;
        frequency = n+0.5*m;
        kernels(causal,n,m) = ...
            exp(-decayRate*lags(causal)).* ...
            cos(2*pi*frequency*lags(causal));
    end
end

figure
tiledlayout(N,M)
for n = 1:N
    for m = 1:M
        nexttile
        plot(lags,kernels(:,n,m),'LineWidth',1.2)
        xline(0,':k')
        title(sprintf('Output %d, input %d',n,m))
        xlabel('Lag')
        grid on
    end
end

%% Construct the operator
% convMulti computes the discrete convolution sum. Scale=dt supplies the quadrature factor required by the time-convolution approximation.

A = ConvMultiOp(kernels,Scale=dt);

fprintf('Domain dimension: %d = K*M\n',A.dimDomain)
fprintf('Range dimension:  %d = K*N\n',A.dimRange)

%% Apply the forward operator
% The K-by-M matrix W stores one time function per column. The LinearOp interface uses w=W(:), which stacks the complete input channels. The output vector similarly contains N complete output channels.

W = zeros(K,M);
W(:,1) = exp(-((times-0.25)/0.06).^2);
W(:,2) = (1-2*(pi*8*(times-0.45)).^2).* ...
    exp(-(pi*8*(times-0.45)).^2);

w = W(:);
d = A.fwd(w);
D = reshape(d,K,N);

figure
tiledlayout(2,1)
nexttile
plot(times,W,'LineWidth',1.2)
xlabel('Time')
title('Input channels')
legend('w_1','w_2')
grid on

nexttile
plot(times,D,'LineWidth',1.2)
xlabel('Time')
title('Output channels')
legend('d_1','d_2','d_3')
grid on

%% Verify the matrix representation
% matrixRep constructs the block-Toeplitz matrix in the standard bases. It is useful for tests and small examples, while the matrix-free actions are preferable for large inverse problems.

Amat = A.matrixRep();
matrixError = norm(A.fwd(w)-Amat*w)/norm(Amat*w);
fprintf('Relative matrix/forward error: %.3e\n',matrixError)

figure
imagesc(Amat)
axis image
colorbar
xlabel('Domain coordinate')
ylabel('Range coordinate')
title('Block-Toeplitz matrix representation')

%% Verify the adjoint
% For Euclidean inner products, the defining identity is dot(A.fwd(w),y) = dot(w,A.adj(y)). This is the key verification before using the operator in CGLS or another iterative inversion method.

rng(4)
y = randn(A.dimRange,1);
lhs = dot(A.fwd(w),y);
rhs = dot(w,A.adj(y));
adjointError = abs(lhs-rhs)/max([1,abs(lhs),abs(rhs)]);

fprintf('Relative adjoint error: %.3e\n',adjointError)
