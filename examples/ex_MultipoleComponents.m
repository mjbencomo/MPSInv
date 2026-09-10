%% Multipole components
% A multipole component contains only the time-independent description of
% one source channel. A time function is attached when a MultipoleTerm is
% needed for a simulation.

clear

%% One-dimensional components
xc = 0.4; %source location
monopole1D = MultipoleComponent1D( ...
    xc,0, ...
    ApproximationOrder=4, ...
    Label="Pressure monopole");

dipole1D = MultipoleComponent1D( ...
    xc,1, ...
    TargetField="velocity", ...
    ApproximationOrder=4, ...
    Label="Velocity dipole");

% The same component can be paired with different time functions to
% generate a MultipoleTerm
gaussian = @(t) exp(-((t-0.35)/0.1).^2);
term1D = monopole1D.withTimeFunction(gaussian,Amplitude=2);

fprintf("%s: value at t=0.35 is %.6g\n", ...
    monopole1D.label,term1D.evaluateTime(0.35));

%% Constructing an impulse source for a Green-function calculation
dt = 1e-3;
deltaApproximation = @(t) deltaTime(t,0,dt);
greenTerm1D = monopole1D.withTimeFunction(deltaApproximation);

% The existing source construction remains unchanged.
xPrimal = Grid1DUnifPrimal([0,1],201);
pressureGrid1D = GridSpace1D(xPrimal);
velocityGrid1D = GridSpace1D(Grid1DUnifDual(xPrimal));
greenSource1D = greenTerm1D.discretize( ...
    pressureGrid1D,velocityGrid1D, ...
    ApproximationOrder=monopole1D.approximationOrder);

fprintf("The Green source has %d pressure term.\n", ...
    greenSource1D.numberOfPressureTerms);

%% Two-dimensional components
mixedDerivative2D = MultipoleComponent2D( ...
    [0.45,0.55],[1,1], ...
    TargetField="pressure", ...
    ApproximationOrder=4, ...
    Label="Mixed pressure derivative");

ricker = @(t) (1-2*(pi*20*(t-0.1)).^2).* ...
    exp(-(pi*20*(t-0.1)).^2);
term2D = mixedDerivative2D.withTimeFunction(ricker,Amplitude=-1);

fprintf("%s targets the %s equation.\n", ...
    mixedDerivative2D.label,term2D.targetField);

%% Collections represent input channels
% Components in one collection must have the same concrete class.
components1D = [monopole1D,dipole1D];
fprintf("Number of one-dimensional input channels: %d\n", ...
    numel(components1D));
