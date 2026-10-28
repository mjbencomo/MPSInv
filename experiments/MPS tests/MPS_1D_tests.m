% Testing the discretization of MPS in 1D

clear

domain_x = [0,2];
domain_t = [0,1];

% Setting up spatial grid
Nx=201;
grid_primal = Grid1DUnifPrimal(domain_x,Nx);
grid_dual = Grid1DUnifDual(grid_primal);
gsp_P = GridSpace1D(grid_primal);
gsp_V = GridSpace1D(grid_dual);
dx = gsp_P.h;
% Acoustic medium parameters
medium = AcousticMedium1D.constant(BetaValue=1,KappaValue=1);
[kappa,beta,wspeed,cMax] = sampleMedium1D(gsp_P,medium);

% Setting up time grid
CFL_factor = 0.9;
dt = CFL_factor*(6/7)*dx/cMax;
Nt = ceil(diff(domain_t)/dt)+1;
dt = diff(domain_t)/(Nt-1); %recomputing to ensure time domain divisible by dt
t_pts = domain_t(1):dt:domain_t(2);

% Output spatial grid and time points
gsp_out   = GridSpace1D(Grid1DUnifPrimal(domain_x,101));
t_pts_out = linspace(domain_t(1),domain_t(2),51);

% Initial conditions
p0 = @(x) zeros(size(x));
v0 = @(x) zeros(size(x));

% Initializing PML
pml = PML1D(WidthLeft=0.2,WidthRight=0);

% Initialize solver
solver = AcousticSolver1D_PML( ...
    gsp_P,t_pts, ...
    gsp_out,t_pts_out, ...
    medium, pml, ...
    SpatialOrder=4);

%% Setting up MPS

x_src = 0.4;
d_src = 0;
mps_comp = MultipoleComponent1D(x_src,d_src,TargetField="pressure");

w_src = @(t) exp(-((t-0.2)/0.1).^2);
mps_term = MultipoleSrcTerm1D(mps_comp,w_src);

appxOrder = 4;
mps_src = mps_term.discretize(gsp_P,gsp_V,appxOrder);

%% Running solver
[P,~] = solver.solve(p0,v0,Source=mps_src);

%% Plotting

P.animate();

tau = @(x,t) t-abs(x)/cMax;
true_soln = @(x,t) w_src(tau(x,t))/(2*cMax);
P_true = GridFunctionTimeSeries(gsp_out,t_pts_out,)