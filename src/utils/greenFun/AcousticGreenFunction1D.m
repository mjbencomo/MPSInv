classdef AcousticGreenFunction1D < GreenFunction
    % AcousticGreenFunction1D computes a 1D acoustic impulse response.

    methods
        function obj = AcousticGreenFunction1D(component,solver,options)
            arguments
                component (1,1) MultipoleComponent1D
                solver
                options.OutputField (1,1) string ...
                    {mustBeMember(options.OutputField, ...
                    ["pressure","velocity"])} = "pressure"
            end

            if ~isa(solver,'AcousticSolver') || ~isscalar(solver)
                error('AcousticGreenFunction1D:InvalidSolver', ...
                    'solver must be a scalar AcousticSolver object.');
            end
            if solver.computationalGrid.dim ~= 1 || solver.outputGrid.dim ~= 1
                error('AcousticGreenFunction1D:SolverDimensionMismatch', ...
                    'The injected solver must use one-dimensional grids.');
            end
            obj@GreenFunction(component,solver,options.OutputField);
        end

        function compute(obj,approximationOrder)
            arguments
                obj
                approximationOrder (1,1) double ...
                    {mustBeInteger,mustBePositive}
            end

            computationalGrid = obj.solver.computationalGrid;
            if ~isa(computationalGrid.x,'Grid1DUnifPrimal')
                error('AcousticGreenFunction1D:PrimalGridRequired', ...
                    ['computationalGrid must contain a uniform primal ' ...
                     'grid.']);
            end

            velocityGrid = GridSpace1D( ...
                Grid1DUnifDual(computationalGrid.x));
            impulse = @(t) deltaTime(t,0,obj.dt);
            term = MultipoleSrcTerm1D(obj.component,impulse);
            source = term.discretize( ...
                computationalGrid,velocityGrid,approximationOrder);

            zero = @(x) zeros(size(x));
            [pressure,velocity] = obj.solver.solve(zero,zero,Source=source);

            if obj.outputField == "pressure"
                selectedOutput = pressure;
            else
                selectedOutput = velocity;
            end

            obj.storeComputedValues( ...
                obj.extendSolverOutput(selectedOutput));
        end
    end
end
