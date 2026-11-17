function plotGreenCGLS2D(results)
%PLOTGREENCGLS2D Plot experiment diagnostics returned by runGreenCGLS2D.
%   PLOTGREENCGLS2D(RESULTS) creates the source, geometry, pressure,
%   Green-kernel, operator, and CGLS convergence figures for one experiment.

arguments
    results (1,1) struct
end

requiredFields = { ...
    'receiverTimes','sourceTimeFunctions','sourceLocations', ...
    'receiverLocations','wSpeed','greenLags','greenSourceTimes', ...
    'greenKernels', ...
    'derivativeOrders','directPressureValues','trueSourceSamples', ...
    'recoveredSourceSamples','recoveredDataByReceiver','recoveredData', ...
    'noisyData','greenOperator','cglsInfo'};
missingFields = requiredFields(~isfield(results,requiredFields));
if ~isempty(missingFields)
    error('plotGreenCGLS2D:MissingResultFields', ...
        'The results structure is missing required field(s): %s.', ...
        strjoin(missingFields,', '));
end

t_rcv = results.receiverTimes(:);
t_src = results.greenSourceTimes(:);
t_ext = results.greenLags(:);
sourceTimeFunctions = results.sourceTimeFunctions;
greenKernels = results.greenKernels;
derivativeOrders = results.derivativeOrders;
W_src = results.trueSourceSamples;
W_cgls = results.recoveredSourceSamples;
P_direct = results.directPressureValues.';       % time x receiver
P_recovered = results.recoveredDataByReceiver;   % time x receiver
nSources = size(results.sourceLocations,1);
nReceivers = size(results.receiverLocations,1);

if numel(sourceTimeFunctions) ~= nSources || ...
        size(greenKernels,2) ~= nReceivers || ...
        size(greenKernels,3) ~= nSources || ...
        size(derivativeOrders,1) ~= nSources || ...
        size(P_direct,1) ~= numel(t_rcv) || ...
        size(P_direct,2) ~= nReceivers || ...
        ~isequal(size(W_src),size(W_cgls)) || ...
        size(W_src,1) ~= numel(t_src) || size(W_src,2) ~= nSources || ...
        ~isequal(size(P_direct),size(P_recovered))
    error('plotGreenCGLS2D:InconsistentResults', ...
        'The source, receiver, time, and Green-kernel dimensions do not agree.');
end

sourceNames = cell(nSources,1);
for m = 1:nSources
    sourceNames{m} = multipoleName(derivativeOrders(m,:));
end
receiverNames = arrayfun(@(r) sprintf('Receiver %d',r),1:nReceivers, ...
    'UniformOutput',false);

%% Source time functions on the receiver time grid
W_rcv = zeros(numel(t_rcv),nSources);
for m = 1:nSources
    values = sourceTimeFunctions{m}(t_rcv);
    if ~isvector(values) || numel(values) ~= numel(t_rcv) || ...
            ~isreal(values) || any(~isfinite(values(:)))
        error('plotGreenCGLS2D:InvalidSourceFunctionOutput', ...
            'Source function %d must return one finite real value per time.',m);
    end
    W_rcv(:,m) = values(:);
end

figure('Name','Source Time Functions')
plot(t_rcv,W_rcv,'LineWidth',1.2)
xlabel('Time')
ylabel('Source amplitude')
title('Source Time Functions')
legend(sourceNames,'Location','best')
grid on

%% Acquisition geometry
figGeo = figure('Name','Source-Receiver Locations');
axGeo = axes(figGeo);
plotAcqGeo(results.wSpeed,results.sourceLocations, ...
    results.receiverLocations,Parent=axGeo)
title('Source-Receiver Locations')

%% Pressure heatmap: time down the vertical axis, receiver index horizontally
figure('Name','Direct Pressure Heatmap')
imagesc(1:nReceivers,t_rcv,P_direct)
set(gca,'YDir','reverse')
xticks(1:nReceivers)
xlabel('Receiver Index')
ylabel('Time (s)')
title('Pressure')
cb = colorbar;
cb.Label.String = 'Pressure (Pa)';
nColors = 128;
blueToWhite = [linspace(0,1,nColors)', ...
    linspace(0,1,nColors)',ones(nColors,1)];
whiteToRed = [ones(nColors,1), ...
    linspace(1,0,nColors)',linspace(1,0,nColors)'];
colormap([blueToWhite;whiteToRed])
clim([-25 25])

%% Noise-free direct pressure time series, one curve per receiver
figure('Name','Direct Pressure Time Series')
plot(t_rcv,P_direct,'LineWidth',1.2)
xlabel('Time')
ylabel('Pressure')
title('Noise-Free Pressure at Each Receiver')
legend(receiverNames,'Location','best')
grid on

%% Green functions, showing every receiver for each source component
figure('Name','Pressure Green Functions')
for m = 1:nSources
    subplot(nSources,1,m)
    hold on
    for r = 1:nReceivers
        plot(t_ext,greenKernels(:,r,m),'LineWidth',1.1)
    end
    hold off
    xlabel('Time lag')
    ylabel('Pressure')
    title(['Pressure Green''s Function: ',sourceNames{m}])
    legend(receiverNames,'Location','best')
    grid on
end

%% Green operator matrix, separated into one block per source
A_G = results.greenOperator.matrixRep();
nSourceTimes = size(W_src,1);
if size(A_G,2) ~= nSourceTimes*nSources
    error('plotGreenCGLS2D:UnexpectedGreenOperatorSize', ...
        'The Green operator columns do not match source-time components.');
end
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

%% True and recovered source signals
figure('Name','Source Recovery')
for m = 1:nSources
    subplot(nSources,1,m)
    plot(t_src,W_src(:,m))
    hold on
    plot(t_src,W_cgls(:,m),'--')
    hold off
    legend('True Source','Source from CGLS','Location','best')
    xlabel('Time')
    ylabel('Source')
    title(['Source Recovery Comparison: ',sourceNames{m}])
    grid on
end

%% Source reconstruction errors
figure('Name','Source Errors')
for m = 1:nSources
    subplot(nSources,1,m)
    plot(t_src,W_src(:,m)-W_cgls(:,m))
    xlabel('Time')
    ylabel('Source error')
    title(['Source Error: ',sourceNames{m}])
    grid on
end

%% CGLS convergence histories
info = results.cglsInfo;
iterations = 0:info.iterations;
dataScale = norm(results.noisyData,2);
dataResidual = info.dataResidualHistory(:) / dataScale;
noiseLimit = results.noiseNorm / dataScale;
normalResidual = info.normalResidualHistory(:);
trueError = info.trueErrorHistory(:);
historyLength = min([numel(iterations),numel(dataResidual), ...
    numel(normalResidual),numel(trueError)]);
iterations = iterations(1:historyLength);


figure('Name','CGLS Relative Pressure Residual')
semilogy(iterations,dataResidual(1:historyLength))
hold on
if noiseLimit > 0
    yline(noiseLimit,'r--')
end
hold off
xlabel('Iterations')
ylabel('Relative Pressure Residual')
title('Convergence of CGLS')
grid on

figure('Name','CGLS Normal Residual')
semilogy(iterations,normalResidual(1:historyLength))
xlabel('Iterations')
ylabel('Normal Residual')
title('Convergence of CGLS')
grid on

figure('Name','CGLS True Source Error')
semilogy(iterations,trueError(1:historyLength))
xlabel('Iterations')
ylabel('True Error')
title('Convergence of CGLS')
grid on

%% Normalized pressure fit error at the receiver times
pressureScale = max(abs(results.noisyData));
if pressureScale == 0
    pressureError = zeros(size(P_direct));
else
    pressureError = (results.recoveredDataByReceiver - ...
        reshape(results.noisyData,numel(t_rcv),nReceivers))/pressureScale;
end

figure('Name','Receiver Pressure Error')
plot(t_rcv,pressureError)
xlabel('Time')
ylabel('Normalized pressure error')
title('Receiver Pressure Error (Receiver Time Grid)')
legend(receiverNames,'Location','best')
grid on

%% Direct and CGLS-recovered pressure at each receiver
figure('Name','Direct and Recovered Pressure')
for r = 1:nReceivers
    subplot(nReceivers,1,r)
    plot(t_rcv,P_direct(:,r))
    hold on
    plot(t_rcv,P_recovered(:,r),'--')
    hold off
    legend('Direct','Recovered','Location','best')
    xlabel('Time')
    ylabel('Pressure')
    title(sprintf('Pressure at Receiver %d',r))
    grid on
end
end

function name = multipoleName(order)
if isequal(order,[0,0])
    name = 'Monopole';
elseif isequal(order,[1,0])
    name = 'Dipole in X';
elseif isequal(order,[0,1])
    name = 'Dipole in Y';
else
    name = sprintf('Multipole [%d,%d]',order(1),order(2));
end
end
