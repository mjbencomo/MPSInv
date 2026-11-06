function results = runGreenCGLS2D( ...
        receiverLocations,sourceLocations,derivativeOrders, ...
        sourceTimeFunctions,options)
%RUNGREENCGLS2D Simulate and invert 2D acoustic multipole data.
%   RESULTS = RUNGREENCGLS2D(RCVLOC,SRCLOC,ORDERS,W) uses the direct
%   acoustic solver to generate pressure data, optionally adds Gaussian
%   noise, builds a Green-function forward operator, and recovers the
%   source time functions with CGLS.
%
%   RCVLOC is an N-by-2 array of [x,y] receiver coordinates. The receiver
%   points may have any arrangement; they need not form a Cartesian grid.
%
%   SRCLOC is an M-by-2 array of [x,y] source coordinates. ORDERS is an
%   M-by-2 array of nonnegative integer derivative orders. Row j defines
%   the multipole
%
%       d^(ORDERS(j,1)+ORDERS(j,2)) /
%       (dx^ORDERS(j,1) dy^ORDERS(j,2)) delta(x-SRCLOC(j,:)).
%
%   For example, [0,0] is a monopole, [1,0] is an x-dipole, and [0,1] is a
%   y-dipole. 
% 
%   W is a cell array containing one time-function handle per 
%   source component. Each handle must accept a vector of times and 
%   return a vector of the same length. The component order in ORDERS, 
%   SRCLOC, and W is preserved in the forward operator and recovered 
%   source matrix.
%
%   Name-value options:
%       BetaValue           Buoyancy used (default 1)
%       KappaValue          Bulk Modulus (default 1)
%       TFinal              Final time (default 1).
%       DomainX, DomainY    Spatial intervals (defaults [0,1], [0,1.2]).
%       NxDirect, NyDirect  Direct solver grid sizes (defaults 201, 211).
%       NxGreen, NyGreen    Green solver grid sizes (default to direct grid).
%       ReceiverTimes       Receiver sample times (default: 101 uniform
%                           times from 0 to TFinal).
%       NumReceiverTimes    Number of default receiver times (default 101).
%       CFLFactor           CFL multiplier (default 0.8).
%       PMLWidth            Width on each PML side (default 0.2).
%       SpatialOrder        Spatial approximation order (default 4).
%       Noise_Percent       Gaussian noise scale (default 0.01). The
%                           component standard deviation is this fraction
%                           of max(abs(noise-free data)). Set to zero to
%                           use noise-free data.
%       RandomSeed          Optional nonnegative integer seed. When
%                           provided, a local random stream is used.
%       MaxIterations       CGLS iteration limit (default 1000).
%       RelativeTolerance   CGLS relative normal-residual tolerance
%                           (default 1e-6).
%       MakePlots           Make diagnostic plots (default false).
%
%   RESULTS contains the direct and noisy data, Green kernels, forward
%   operators, true and recovered source samples, receiver and source time
%   grids, and CGLS diagnostics. The receiver data vectors are stacked in
%   receiver-channel order, matching reshape(P.values.',[],1).
%
%   Example: three source components and an arbitrary receiver layout.
%
%       loc = [0.33,0.43; 0.62,0.70; 0.48,0.88];
%       ord = [0,0; 1,0; 0,1];
%       w = {@(t) exp(-((t-0.35)/0.1).^2), ...
%            @(t) exp(-((t-0.40)/0.1).^2), ...
%            @(t) exp(-((t-0.45)/0.1).^2)};
%       rcv = [0.2,0.5; 0.5,0.5; 0.8,0.5; 0.5,0.9];
%       result = runGreenCGLS2D(rcv,loc,ord,w,Noise_Percent=0.01);

arguments
    receiverLocations (:,2) double {mustBeReal,mustBeFinite}
    sourceLocations (:,2) double {mustBeReal,mustBeFinite}
    derivativeOrders (:,2) double ...
        {mustBeInteger,mustBeNonnegative,mustBeFinite}
    sourceTimeFunctions (1,:) cell
    options.BetaValue (1,1) double ...
        {mustBeReal,mustBePositive,mustBeFinite} = 1
    options.KappaValue (1,1) double ...
        {mustBeReal,mustBePositive,mustBeFinite} = 1
    options.TFinal (1,1) double ...
        {mustBeReal,mustBePositive,mustBeFinite} = 1
    options.DomainX (1,2) double {mustBeReal,mustBeFinite} = [0,1]
    options.DomainY (1,2) double {mustBeReal,mustBeFinite} = [0,1.2]
    options.NxDirect (1,1) double {mustBeInteger,mustBePositive} = 201
    options.NyDirect (1,1) double {mustBeInteger,mustBePositive} = 211
    options.NxGreen double = []
    options.NyGreen double = []
    options.ReceiverTimes double = []
    options.NumReceiverTimes (1,1) double ...
        {mustBeInteger,mustBePositive} = 101
    options.CFLFactor (1,1) double {mustBePositive,mustBeFinite} = 0.8
    options.PMLWidth (1,1) double {mustBeNonnegative,mustBeFinite} = 0.2
    options.SpatialOrder (1,1) double ...
        {mustBeInteger,mustBePositive} = 4
    options.Noise_Percent (1,1) double ...
        {mustBeNonnegative,mustBeFinite} = 0.01
    options.RandomSeed double = []
    options.MaxIterations (1,1) double ...
        {mustBeInteger,mustBePositive} = 1000
    options.RelativeTolerance (1,1) double ...
        {mustBeNonnegative,mustBeFinite} = 1e-6
    options.MakePlots (1,1) logical = false
end

validateExperimentInputs(receiverLocations,sourceLocations, ...
    derivativeOrders,sourceTimeFunctions,options);

numberOfSources = size(sourceLocations,1);
numberOfReceivers = size(receiverLocations,1);
spatialOrder = options.SpatialOrder;
domainX = options.DomainX;
domainY = options.DomainY;
tFinal = options.TFinal;

if isempty(options.NxGreen)
    NxGreen = options.NxDirect;
else
    NxGreen = options.NxGreen;
end
if isempty(options.NyGreen)
    NyGreen = options.NyDirect;
else
    NyGreen = options.NyGreen;
end

% Initial conditions, PML, receivers, and constant acoustic medium.
p0 = @(x,y) zeros(size(x));
vx0 = @(x,y) zeros(size(x));
vy0 = @(x,y) zeros(size(x));
pml = PML2D( ...
    WidthLeft=options.PMLWidth,WidthRight=options.PMLWidth, ...
    WidthBottom=options.PMLWidth,WidthTop=options.PMLWidth);
gridSp_rcv = ReceiverSpace2D(receiverLocations);
medium = AcousticMedium2D.constant(...
    BetaValue=options.BetaValue,KappaValue=options.KappaValue);

% Construct one multipole component and source term per input row.
sourceComponents = cell(1,numberOfSources);
mps_terms = cell(1,numberOfSources);
for j = 1:numberOfSources
    sourceComponents{j} = MultipoleComponent2D( ...
        sourceLocations(j,:),derivativeOrders(j,:), ...
        TargetField="pressure");
    mps_terms{j} = MultipoleSrcTerm2D( ...
        sourceComponents{j},sourceTimeFunctions{j});
end
mps = MultipoleSrcSeries2D(mps_terms);

% Direct acoustic simulation on its own computational grid.
[gridSp_P,gridSp_Vx,gridSp_Vy,h] = makeAcousticGrids( ...
    domainX,domainY,options.NxDirect,options.NyDirect);
[~,~,~,waveSpeed,cMax] = sampleMedium2D(gridSp_P,medium);
[t_cmp,~] = makeComputationalTimes(tFinal,cMax,h,options.CFLFactor);

if isempty(options.ReceiverTimes)
    t_rcv = linspace(0,tFinal,options.NumReceiverTimes).';
else
    t_rcv = options.ReceiverTimes(:);
end

source = mps.discretize( ...
    gridSp_P,gridSp_Vx,gridSp_Vy,spatialOrder);
solver = AcousticSolver2D_PML( ...
    gridSp_P,t_cmp,gridSp_rcv,t_rcv,medium,pml, ...
    SpatialOrder=spatialOrder);
[P,~,~] = solver.solve(p0,vx0,vy0,Source=source);
data_nf = reshape(P.values.',[],1);

% Add independent, zero-mean Gaussian noise using the stated 1%-of-peak
% convention: noise = randn(size(data_nf))*norm(data_nf,inf)*fraction.
if isempty(options.RandomSeed)
    standardNormalNoise = randn(size(data_nf));
else
    stream = RandStream('mt19937ar','Seed',options.RandomSeed);
    standardNormalNoise = randn(stream,size(data_nf));
end
noise = standardNormalNoise*norm(data_nf,inf)*options.Noise_Percent;
data_noisy = data_nf + noise;

% Build Green responses on an independently selectable spatial grid.
[gridSp_PG,~,~,hG] = makeAcousticGrids( ...
    domainX,domainY,NxGreen,NyGreen);
[~,~,~,~,cMaxG] = sampleMedium2D(gridSp_PG,medium);
[t_cmpG,dt_cmpG] = makeComputationalTimes( ...
    tFinal,cMaxG,hG,options.CFLFactor);
t_cmpG_solver = (-1:(numel(t_cmpG)-1))*dt_cmpG;

solver_green = AcousticSolver2D_PML( ...
    gridSp_PG,t_cmpG_solver,gridSp_rcv,t_cmpG_solver,medium,pml, ...
    SpatialOrder=spatialOrder);
G = AcousticGreenFunctionSeries2D( ...
    sourceComponents,solver_green,OutputField="pressure");
G.compute(spatialOrder);

% Compose source interpolation, multichannel convolution, and receiver
% interpolation into one operator mapping source samples to data samples.
t_src = t_cmpG(:);
T_op_in = TimeInterpOp( ...
    t_src,t_cmpG,NumChannels=numberOfSources);
T_op_out = TimeInterpOp( ...
    t_cmpG,t_rcv,NumChannels=numberOfReceivers);
A_G_op = G.makeLinearOp();
A_op = CompositeLinearOp({T_op_out,A_G_op,T_op_in});

% Sample every input source function on the convolution source grid.
W_src = zeros(numel(t_src),numberOfSources);
for j = 1:numberOfSources
    sourceSamples = sourceTimeFunctions{j}(t_src);
    if ~isvector(sourceSamples) || ~isreal(sourceSamples) || ...
            numel(sourceSamples) ~= numel(t_src) || ...
            any(~isfinite(sourceSamples(:)))
        error('runGreenCGLS2D:InvalidSourceSamples', ...
            ['Each source time function must return one finite value ' ...
             'for each input time.']);
    end
    W_src(:,j) = sourceSamples(:);
end
w_data = W_src(:);
w_initial = zeros(size(w_data));

[w_cgls,info_cgls] = solver_CGLS( ...
    A_op,data_noisy, ...
    InitialGuess=w_initial, ...
    MaxIterations=options.MaxIterations, ...
    RelativeTolerance=options.RelativeTolerance, ...
    TrueSolution=w_data);
data_recovered = A_op.fwd(w_cgls);
W_cgls = reshape(w_cgls,size(W_src));

% Return results and diagnostic quantities without opening figures by
% default, so this function is suitable for experiment loops.
results = struct();
results.receiverLocations = receiverLocations;
results.sourceLocations = sourceLocations;
results.derivativeOrders = derivativeOrders;
results.sourceComponents = sourceComponents;
results.receiverTimes = t_rcv;
results.directComputationalTimes = t_cmp;
results.greenSourceTimes = t_src;
results.greenSolverTimes = t_cmpG_solver;
results.greenTimeStep = dt_cmpG;
results.noiseFraction = options.Noise_Percent;
results.sourceTimeFunctions = sourceTimeFunctions;
results.trueSourceSamples = W_src;
results.recoveredSourceSamples = W_cgls;
results.directPressureValues = P.values;
results.noiseFreeData = data_nf;
results.noise = noise;
results.noisyData = data_noisy;
results.noiseNorm = norm(noise,2);
results.relativeNoiseNorm = safeRelativeNorm(norm(noise,2),norm(data_nf,2));
results.greenKernels = G.values;
results.greenLags = G.lags;
results.greenOperator = A_G_op;
results.forwardOperator = A_op;
results.recoveredData = data_recovered;
results.recoveredDataByReceiver = reshape( ...
    data_recovered,numel(t_rcv),numberOfReceivers);
results.cglsInfo = info_cgls;
results.relativeSourceError = safeRelativeNorm( ...
    norm(w_cgls-w_data,2),norm(w_data,2));
results.relativeErrorAgainstNoisyData = safeRelativeNorm( ...
    norm(data_recovered-data_noisy,2),norm(data_noisy,2));
results.relativeErrorAgainstNoiseFreeData = safeRelativeNorm( ...
    norm(data_recovered-data_nf,2),norm(data_nf,2));
results.wSpeed = waveSpeed;

if options.MakePlots
    plotExperimentResults(t_rcv,data_nf,data_noisy,noise, ...
        W_src,W_cgls,t_src,numberOfReceivers,numberOfSources);
end
end

function validateExperimentInputs( ...
        receiverLocations,sourceLocations,derivativeOrders, ...
        sourceTimeFunctions,options)
if isempty(receiverLocations)
    error('runGreenCGLS2D:NoReceivers', ...
        'Provide at least one receiver location.');
end
if isempty(sourceLocations)
    error('runGreenCGLS2D:NoSources', ...
        'Provide at least one source component.');
end
if size(derivativeOrders,1) ~= size(sourceLocations,1)
    error('runGreenCGLS2D:SourceCountMismatch', ...
        'sourceLocations and derivativeOrders must have the same row count.');
end
if numel(sourceTimeFunctions) ~= size(sourceLocations,1)
    error('runGreenCGLS2D:TimeFunctionCountMismatch', ...
        'Provide exactly one time-function handle for each source component.');
end
if any(~cellfun(@(f) isa(f,'function_handle'),sourceTimeFunctions))
    error('runGreenCGLS2D:InvalidTimeFunction', ...
        'Every entry in sourceTimeFunctions must be a function handle.');
end
if any(options.DomainX(2) <= options.DomainX(1)) || ...
        any(options.DomainY(2) <= options.DomainY(1))
    error('runGreenCGLS2D:InvalidDomain', ...
        'DomainX and DomainY must each be increasing intervals.');
end
if options.CFLFactor > 1
    error('runGreenCGLS2D:InvalidCFLFactor', ...
        'CFLFactor must be positive and no greater than 1.');
end
if any(receiverLocations(:,1) < options.DomainX(1) | ...
        receiverLocations(:,1) > options.DomainX(2) | ...
        receiverLocations(:,2) < options.DomainY(1) | ...
        receiverLocations(:,2) > options.DomainY(2))
    error('runGreenCGLS2D:ReceiverOutsideDomain', ...
        'All receiver locations must lie within DomainX and DomainY.');
end
if any(sourceLocations(:,1) < options.DomainX(1) | ...
        sourceLocations(:,1) > options.DomainX(2) | ...
        sourceLocations(:,2) < options.DomainY(1) | ...
        sourceLocations(:,2) > options.DomainY(2))
    error('runGreenCGLS2D:SourceOutsideDomain', ...
        'All source locations must lie within DomainX and DomainY.');
end
if options.NxDirect < 5 || options.NyDirect < 5
    error('runGreenCGLS2D:GridTooSmall', ...
        'Direct grid sizes must each be at least 5.');
end
if ~isempty(options.NxGreen) && ...
        (~isreal(options.NxGreen) || ~isscalar(options.NxGreen) || ...
         ~isfinite(options.NxGreen) || options.NxGreen < 5 || ...
         options.NxGreen ~= fix(options.NxGreen))
    error('runGreenCGLS2D:InvalidGreenGridSize', ...
        'NxGreen must be empty or an integer of at least 5.');
end
if ~isempty(options.NyGreen) && ...
        (~isreal(options.NyGreen) || ~isscalar(options.NyGreen) || ...
         ~isfinite(options.NyGreen) || options.NyGreen < 5 || ...
         options.NyGreen ~= fix(options.NyGreen))
    error('runGreenCGLS2D:InvalidGreenGridSize', ...
        'NyGreen must be empty or an integer of at least 5.');
end
if options.NumReceiverTimes < 2 && isempty(options.ReceiverTimes)
    error('runGreenCGLS2D:TooFewReceiverTimes', ...
        'NumReceiverTimes must be at least 2.');
end
if ~isempty(options.ReceiverTimes)
    tr = options.ReceiverTimes(:);
    if ~isvector(options.ReceiverTimes) || ~isreal(tr) || ...
            any(~isfinite(tr)) || any(diff(tr) <= 0) || ...
            tr(1) < 0 || tr(end) > options.TFinal
        error('runGreenCGLS2D:InvalidReceiverTimes', ...
            'ReceiverTimes must increase strictly within [0,TFinal].');
    end
end
if options.PMLWidth*2 >= min(diff(options.DomainX),diff(options.DomainY))
    error('runGreenCGLS2D:InvalidPMLWidth', ...
        'PMLWidth must be less than half of each domain dimension.');
end
if ~isempty(options.RandomSeed) && ...
        (~isreal(options.RandomSeed) || ~isscalar(options.RandomSeed) || ...
         ~isfinite(options.RandomSeed) || ...
         options.RandomSeed < 0 || options.RandomSeed ~= fix(options.RandomSeed))
    error('runGreenCGLS2D:InvalidRandomSeed', ...
        'RandomSeed must be empty or a nonnegative integer scalar.');
end
end

function [gridSp_P,gridSp_Vx,gridSp_Vy,h] = ...
        makeAcousticGrids(domainX,domainY,Nx,Ny)
grid_x_primal = Grid1DUnifPrimal(domainX,Nx);
grid_y_primal = Grid1DUnifPrimal(domainY,Ny);
grid_x_dual = Grid1DUnifDual(grid_x_primal);
grid_y_dual = Grid1DUnifDual(grid_y_primal);
gridSp_P = GridSpace2D(grid_x_primal,grid_y_primal);
gridSp_Vx = GridSpace2D(grid_x_dual,grid_y_primal);
gridSp_Vy = GridSpace2D(grid_x_primal,grid_y_dual);
h = gridSp_P.h;
end

function [t,dt] = makeComputationalTimes(tFinal,cMax,h,CFLFactor)
dt = CFLFactor*(6/7)/(cMax*sqrt(sum(1./h.^2)));
numberOfTimes = ceil(tFinal/dt)+1;
dt = tFinal/(numberOfTimes-1);
t = (0:numberOfTimes-1)'*dt;
end

function value = safeRelativeNorm(numerator,denominator)
if denominator == 0
    value = NaN;
else
    value = numerator/denominator;
end
end

function plotExperimentResults( ...
        t_rcv,data_nf,data_noisy,noise,W_src,W_cgls,t_src, ...
        numberOfReceivers,numberOfSources)
data_nf_by_receiver = reshape(data_nf,numel(t_rcv),numberOfReceivers);
data_noisy_by_receiver = reshape( ...
    data_noisy,numel(t_rcv),numberOfReceivers);
noise_by_receiver = reshape(noise,numel(t_rcv),numberOfReceivers);

figure
for r = 1:numberOfReceivers
    subplot(numberOfReceivers,2,2*r-1)
    plot(t_rcv,data_nf_by_receiver(:,r));
    hold on
    plot(t_rcv,data_noisy_by_receiver(:,r));
    hold off
    legend('Noise-free data','Noisy data')
    xlabel('Time')
    ylabel('Pressure')
    title(sprintf('Receiver %d data',r))

    subplot(numberOfReceivers,2,2*r)
    plot(t_rcv,noise_by_receiver(:,r));
    xlabel('Time')
    ylabel('Noise')
    title(sprintf('Receiver %d added noise',r))
end

figure
for j = 1:numberOfSources
    subplot(numberOfSources,1,j)
    plot(t_src,W_src(:,j));
    hold on
    plot(t_src,W_cgls(:,j));
    hold off
    legend('True source','CGLS estimate')
    xlabel('Time')
    ylabel('Source amplitude')
    title(sprintf('Source component %d',j))
end
end
