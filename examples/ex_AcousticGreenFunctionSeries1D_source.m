%% Computing a 1D acoustic Green-function series
% This example computes the pressure Green functions for two multipole
% components at three receiver locations. The resulting ConvMultiOp maps
% two source-time-function channels to three receiver-data channels.

clear
close all

%% Configure the acoustic solver
% Every component in a GreenFunctionSeries uses the same configured solver,
% output grid, time grid, and observed field.

domain = [0,1];
computationalGrid = GridSpace1D( ...
    Grid1DUnifPrimal(domain,201));

receiverLocations = [0.2,0.5,0.8];
outputGrid = ReceiverSpace1D(receiverLocations);

medium = AcousticMedium1D.constant(KappaValue=1,BetaValue=1);
pml = PML1D(WidthLeft=0.2,WidthRight=0.2);

[~,~,wspeed,cMax] = sampleMedium1D(computationalGrid,medium);
T = 1;
dtEstimate = 0.8*(6/7)*computationalGrid.h/cMax;
K = ceil(T/dtEstimate)+1;
times = linspace(0,T,K);
solverTimes = [-times(2),times];

solver = AcousticSolver1D_PML( ...
    computationalGrid,solverTimes,outputGrid,solverTimes,medium,pml, ...
    SpatialOrder=4);

%% Define the input channels
% The component order determines the input-channel order. Here channel 1 is
% a pressure monopole and channel 2 is a first spatial derivative source.
sourceLocations = [0.4;0.6];
components = { ...
    MultipoleComponent1D(sourceLocations(1),0,TargetField="pressure"), ...
    MultipoleComponent1D(sourceLocations(2),1,TargetField="pressure")};

G_series = AcousticGreenFunctionSeries1D( ...
    components,solver,OutputField="pressure");

%% Plotting source-receiver geometry
plotAcqGeo(wspeed,sourceLocations,receiverLocations);
title("Source-Receiver Locations")

%% Compute and inspect the kernels
% A scalar approximation order is applied to every component. A vector with
% one entry per input channel may be supplied instead.

G_series.compute(4);

fprintf('Kernel array size: %d-by-%d-by-%d\n',size(G_series.values))
fprintf('Input channels: %d\n',G_series.numberOfInputChannels)
fprintf('Output channels: %d\n',G_series.numberOfOutputChannels)

figure
tiledlayout(G_series.numberOfInputChannels,1)
for m = 1:G_series.numberOfInputChannels
    nexttile
    plot(G_series.lags,G_series.values(:,:,m),'LineWidth',1.1)
    xline(0,':k')
    xlabel('Time lag')
    ylabel('Pressure kernel')
    title(sprintf('Input channel %d: %s',m, ...
        components{m}.description()))
    grid on
end

%% Construct and apply the multichannel operator
% Channel-block vectorization is used. If W has size K-by-M, then W(:)
% stores all samples of input channel 1, followed by input channel 2.

A = G_series.makeLinearOp();

W = zeros(K,G_series.numberOfInputChannels);
W(:,1) = exp(-((times(:)-0.25)/0.06).^2);
W(:,2) = 0.4*exp(-((times(:)-0.55)/0.08).^2);

d = A.fwd(W(:));
D = reshape(d,K,G_series.numberOfOutputChannels);

Amat = A.matrixRep();

figure
imagesc(Amat)
axis xy
axis equal
title('Convolution operator matrix')
colorbar
xlabel('column index')
ylabel('row index')
xlim([0,size(Amat,2)])
ylim([0,size(Amat,1)])

figure
tiledlayout(2,1)
nexttile
plot(times,W,'LineWidth',1.2)
xlabel('Time')
ylabel('Source amplitude')
legend('input 1','input 2',Location='best')
title('Input channels')
grid on

nexttile
plot(times,D,'LineWidth',1.2)
xlabel('Time')
ylabel('Pressure')
legend(compose('x_r = %.1f',receiverLocations),Location='best')
title('Receiver data')
grid on

%% Verify the discrete adjoint
% rng(5)
% y = randn(A.dimRange,1);
% lhs = dot(A.fwd(W(:)),y);
% rhs = dot(W(:),A.adj(y));
% relativeAdjointError = abs(lhs-rhs)/max([1,abs(lhs),abs(rhs)]);
% fprintf('Relative adjoint error: %.3e\n',relativeAdjointError)
