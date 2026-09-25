classdef AcousticGreenFunction2D < GreenFunction
    % AcousticGreenFunction2D computes a 2D acoustic impulse response.

    methods
        function obj = AcousticGreenFunction2D(component,solver,options)
            arguments
                component (1,1) MultipoleComponent2D
                solver
                options.OutputField (1,1) string ...
                    {mustBeMember(options.OutputField, ...
                    ["pressure","velocityX","velocityY"])} = "pressure"
            end

            if ~isa(solver,'AcousticSolver') || ~isscalar(solver)
                error('AcousticGreenFunction2D:InvalidSolver', ...
                    'solver must be a scalar AcousticSolver object.');
            end
            if solver.computationalGrid.dim ~= 2 || solver.outputGrid.dim ~= 2
                error('AcousticGreenFunction2D:SolverDimensionMismatch', ...
                    'The injected solver must use two-dimensional grids.');
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
            if ~isa(computationalGrid.x,'Grid1DUnifPrimal') || ...
                    ~isa(computationalGrid.y,'Grid1DUnifPrimal')
                error('AcousticGreenFunction2D:PrimalGridRequired', ...
                    ['computationalGrid must contain uniform primal ' ...
                     'grids in both directions.']);
            end

            velocityXGrid = GridSpace2D( ...
                Grid1DUnifDual(computationalGrid.x), ...
                computationalGrid.y);
            velocityYGrid = GridSpace2D( ...
                computationalGrid.x, ...
                Grid1DUnifDual(computationalGrid.y));

            impulse = @(t) deltaTime(t,0,obj.dt);
            term = MultipoleSrcTerm2D(obj.component,impulse);
            source = term.discretize( ...
                computationalGrid,velocityXGrid,velocityYGrid, ...
                approximationOrder);

            zero = @(x,y) zeros(size(x));
            [pressure,velocityX,velocityY] = obj.solver.solve( ...
                zero,zero,zero,Source=source);

            switch obj.outputField
                case "pressure"
                    selectedOutput = pressure;
                case "velocityX"
                    selectedOutput = velocityX;
                case "velocityY"
                    selectedOutput = velocityY;
            end

            obj.storeComputedValues( ...
                obj.extendSolverOutput(selectedOutput));
        end
    end
end
