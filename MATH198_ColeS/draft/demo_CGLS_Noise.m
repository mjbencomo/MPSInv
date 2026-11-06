% Demo: Green's Functions in 2D Acoustics Using CGLS with Interesting Multipoles
% Converted from demo_CGLS_TimeInterp.mlx to a regular MATLAB script.
% This version adds 1% Gaussian noise to the direct pressure data.

%% Initial set up

clear

t_final = 1;
domain_t = [0,t_final];
domain_x = [0,1];
domain_y = [0,1.2];


% Initial conditions
p0 = @(x,y) zeros(size(x));
vx0 = @(x,y) zeros(size(x));
vy0 = @(x,y) zeros(size(x));

% Initializing PML
pml = PML2D(WidthLeft=0.2,WidthRight=0.2,WidthBottom=0.2,WidthTop=0.2);

% Receiver Setup
receiverLocations = [
    0.2, 0.5;
    0.5, 0.5;
    0.8, 0.5];
gridSp_rcv = ReceiverSpace2D(receiverLocations);
nReceivers = gridSp_rcv.numPts;

% Continuous acoustic medium object function
medium = AcousticMedium2D.constant(BetaValue=1,KappaValue=1);

%% Computing pressure data directly

% Spatial computational grids and grid spaces
Nx=201;
Ny=211;
grid_x_primal = Grid1DUnifPrimal(domain_x,Nx);
grid_y_primal = Grid1DUnifPrimal(domain_y,Ny);
grid_x_dual   = Grid1DUnifDual(grid_x_primal);
grid_y_dual   = Grid1DUnifDual(grid_y_primal);

gridSp_P = GridSpace2D(grid_x_primal, grid_y_primal);   %grid space for pressure
gridSp_Vx = GridSpace2D(grid_x_dual,  grid_y_primal);   %grid spaces for velocity in x
gridSp_Vy = GridSpace2D(grid_x_primal,  grid_y_dual);   %grid spaces for velocity in y

h = gridSp_P.h; %spatial grid sizes, h = [dx,dy]

% Sampling medium parameters to compute max wave speed 
[~,~,~,wspeed,cMax] = sampleMedium2D(gridSp_P,medium);

% Determining computational time grid size
CFL_factor = 0.8;
dt_cmp = CFL_factor*(6/7)/(cMax*sqrt(sum(1./h.^2))); %CFL condition for 2-4 FDM
Nt_cmp = ceil(diff(domain_t)/dt_cmp)+1;
dt_cmp = diff(domain_t)/(Nt_cmp-1); %recomputing to ensure time domain divisible by dt

% Computational time points
t_cmp = domain_t(1):dt_cmp:domain_t(2);

% Receiver time points
Nt_rcv = 101;
t_rcv  = linspace(domain_t(1),domain_t(2),Nt_rcv);
dt_rcv = diff(domain_t)/(Nt_rcv-1);

location0 = [0.33,0.43];
sourceLocation = [location0; location0; location0];
tc0 = 0.35;
tc1 = 0.35;
tc2 = 0.35;
s0 = [0,0];
s1 = [1,0];
s2 = [0,1];
w0 = @(t) exp(-((t-tc0)/0.1).^2);
w1 = @(t) exp(-((t-tc1)/0.1).^2);
w2 = @(t) exp(-((t-tc2)/0.1).^2);

monopole = MultipoleComponent2D(location0,s0,TargetField="pressure");
dipole_x = MultipoleComponent2D(location0,s1,TargetField="pressure");
dipole_y = MultipoleComponent2D(location0,s2,TargetField="pressure");
sourceComponents = {monopole,dipole_x,dipole_y};

% We spit the source into 3 terms
term0 = MultipoleSrcTerm2D(monopole,w0);
term1 = MultipoleSrcTerm2D(dipole_x,w1);
term2 = MultipoleSrcTerm2D(dipole_y,w2);
mps = MultipoleSrcSeries2D({term0,term1,term2});
spatialOrder = 4;  
source = mps.discretize(gridSp_P,gridSp_Vx,gridSp_Vy,spatialOrder);

W_rcv = [w0(t_rcv(:)), w1(t_rcv(:)), w2(t_rcv(:))];

figure
plot(t_rcv,W_rcv)
xlabel("time")
title("Source time function");
figure
plotAcqGeo(wspeed,sourceLocation,receiverLocations);
title("Source-Receiver Locations")

solver = AcousticSolver2D_PML( ...
    gridSp_P,   t_cmp, ...
    gridSp_rcv, t_rcv, ...
    medium, pml, ...
    SpatialOrder=spatialOrder);

[P,~,~] = solver.solve(p0,vx0,vy0,Source=source);
data_nf = reshape(P.values.',[],1);

P.plotTimeSeries();

% Add 1% Gaussian noise:
% d_noisy = d_noise_free + randn(size(d_noise_free))*||d_noise_free||_inf*0.01
noise_percent = 0.01;
noise = randn(size(data_nf))*norm(data_nf,inf)*noise_percent;
data_noisy = data_nf + noise;

% % Compare the clean data, added noise, and noisy data by receiver.
% data_nf_by_receiver = reshape(data_nf,Nt_rcv,nReceivers);
% noise_by_receiver = reshape(noise,Nt_rcv,nReceivers);
% data_noisy_by_receiver = reshape(data_noisy,Nt_rcv,nReceivers);
% figure
% for r = 1:nReceivers
%     subplot(nReceivers,2,2*r-1)
%     plot(t_rcv,data_nf_by_receiver(:,r));
%     hold on
%     plot(t_rcv,data_noisy_by_receiver(:,r));
%     hold off
%     legend('Noise-free data','Noisy data')
%     xlabel('Time')
%     ylabel('Pressure')
%     title(sprintf('Data: receiver %d',r))
% 
%     subplot(nReceivers,2,2*r)
%     plot(t_rcv,noise_by_receiver(:,r));
%     xlabel('Time')
%     ylabel('Noise')
%     title(sprintf('Added Gaussian noise: receiver %d',r))
% end

%% Computing pressure data via convolution

% Spatial computational grids and grid spaces
NxG=Nx;
NyG=Ny;

grid_x_primalG = Grid1DUnifPrimal(domain_x,NxG);
grid_y_primalG = Grid1DUnifPrimal(domain_y,NyG);
grid_x_dualG   = Grid1DUnifDual(grid_x_primalG);
grid_y_dualG   = Grid1DUnifDual(grid_y_primalG);

gridSp_PG = GridSpace2D(grid_x_primalG, grid_y_primalG);    %grid space for pressure
gridSp_VxG = GridSpace2D(grid_x_dualG,  grid_y_primalG);    %grid spaces for velocity in x
gridSp_VyG = GridSpace2D(grid_x_primalG,  grid_y_dualG);    %grid spaces for velocity in y
hG = gridSp_PG.h;    %spatial grid sizes, h = [dx,dy]

% Sampling medium parameters to compute max wave speed 
[~,~,~,~,cMaxG] = sampleMedium2D(gridSp_PG,medium);

% Determining computational time grid size
CFL_factor = 0.8;
dt_cmpG = CFL_factor*(6/7)/(cMaxG*sqrt(sum(1./hG.^2))); %CFL condition for 2-4 FDM
Nt_cmpG = ceil(diff(domain_t)/dt_cmpG)+1;
dt_cmpG = diff(domain_t)/(Nt_cmpG-1); %recomputing to ensure time domain divisible by dt

% Computational time points
t_cmpG = 0:dt_cmpG:t_final;
% Green solver times: one step before zero through t_final
t_cmpG_solver = (-1:Nt_cmpG-1)*dt_cmpG;

% Initialize solver for direct data.
solver_green = AcousticSolver2D_PML( ...
    gridSp_PG,  t_cmpG_solver, ...
    gridSp_rcv, t_cmpG_solver, ...
    medium, pml, ...
    SpatialOrder=spatialOrder);

G = AcousticGreenFunctionSeries2D( ...
    sourceComponents,solver_green,OutputField="pressure");

% Compute the impulse (Instead of the delta_t line we had previously)
G.compute(spatialOrder);

data_green = G.values;
t_ext = G.lags;

figure
subplot(3,1,1)
plot(t_ext,data_green(:,1,1))
xlabel("time")
title("Pressure Green's Function: Monopole")
subplot(3,1,2)
plot(t_ext,data_green(:,1,2))
xlabel("time")
title("Pressure Green's Function: x-Dipole")
subplot(3,1,3)
plot(t_ext,data_green(:,1,3))
xlabel("time")
title("Pressure Green's Function: y-Dipole")

% Constructing convolution matrix (We do this differently now that we have
% functions that can manipulate linear operators)
A_G_op = G.makeLinearOp();
A_G = A_G_op.matrixRep();

figure
K = size(A_G,2)/3;
for m = 1:3
    columns = (m-1)*K+(1:K);
    A_component = A_G(:,columns);

    subplot(1,3,m)
    imagesc(A_component)
    axis image
    colorbar
    xlabel('column index')
    ylabel('row index')
    title(sprintf('A_%d',m-1))
end

% t_src defines the time grid for source samples used by convolution
t_src = t_cmpG;

% Here we interpolate pressure from t_cmpG to t_rcv using TimeInterpOp
T_op_in = TimeInterpOp(t_src,t_cmpG,NumChannels=G.numberOfInputChannels);
T_op_out = TimeInterpOp(t_cmpG,t_rcv,NumChannels=G.numberOfOutputChannels);

% Composition of Linear Operators: 
A_op = CompositeLinearOp({T_op_out,A_G_op,T_op_in});

% Computing pressure data via convolution
W_src = [w0(t_src(:)), w1(t_src(:)), w2(t_src(:))];
w_data = W_src(:);

% Build a zero source estimate of the true solution
w_initial = zeros(size(w_data));

% Run the CGLS Solver
[w_cgls,info_cgls] = solver_CGLS(...
    A_op,data_noisy,...
    InitialGuess=w_initial,...
    MaxIterations=1000, ...
    TrueSolution=w_data);

% Computing pressure from the newly recovered source the time grid
data_recovered = A_op.fwd(w_cgls);

% Computing error of calculated source
source_error = norm(w_cgls-w_data)/norm(w_data);
fprintf('Relative source error: %.3e\n',source_error);

% Computing error of calculated pressure data (receiver time grid)
data_error = norm(data_recovered-data_noisy)/norm(data_noisy);
fprintf('Relative error against noisy pressure data: %.3e\n',data_error);

data_error_clean = norm(data_recovered-data_nf)/norm(data_nf);
fprintf('Relative error against noise-free pressure data: %.3e\n', ...
    data_error_clean);

% Checking the number of iterations:
fprintf('CGLS iterations: %d\n',info_cgls.iterations);
disp(info_cgls.message)

% Now let's plot some graphs to better show these results
% Comparing the calculated source against the true source
W_cgls = reshape(w_cgls,size(W_src));
graph_names = {'Monopole','Dipole in X','Dipole in Y'};

figure
for m = 1:3
    subplot(3,1,m)
    plot(t_src,W_src(:,m));
    hold on
    plot(t_src,W_cgls(:,m));
    hold off

    legend('True Source','Source from CGLS')
    xlabel('Time')
    ylabel('Source')
    title(['Source Recovery Comparison: ',graph_names{m}])
end

figure
for m = 1:3
    subplot(3,1,m)
    plot(t_src,W_src(:,m) - W_cgls(:,m));

    xlabel('Time')
    ylabel('Source Error')
    title(['Source Error: ',graph_names{m}])
end

% Convergence Speed
iterations_cgls = 0:info_cgls.iterations;
relative_residual = info_cgls.dataResidualHistory/norm(data_noisy);
abs_normalResidual = info_cgls.normalResidualHistory;
trueErrorHistory = info_cgls.trueErrorHistory;


figure
semilogy(iterations_cgls,relative_residual)
xlabel('Iterations')
ylabel('Relative Pressure Residual')
title('Convergence of CGLS')

figure
semilogy(iterations_cgls,abs_normalResidual)
xlabel('Iterations')
ylabel('Normal Residual')
title('Convergence of CGLS')

figure
semilogy(iterations_cgls,trueErrorHistory)
xlabel('Iterations')
ylabel('True Error')
title('Convergence of CGLS')
data_res = reshape( ...
    (data_recovered-data_noisy) / max(abs(data_noisy)), ...
    numel(t_rcv),nReceivers);

% Plotting Pressure Data Error (Receiver Time Grid)
figure
plot(t_rcv,data_res)
hold on
xlabel('Time')
title('Receiver Pressure Error (Receiver Time Grid)')
P_recovered = reshape(data_recovered,size(P.values));

figure
for m = 1:3
    subplot(3,1,m)
    plot(t_rcv,P.values(m,:));
    hold on
    plot(t_rcv,P_recovered(m,:),'--')
    hold off
    legend('Direct','Recovered')
    xlabel('Time')
    ylabel('Pressure')
end
