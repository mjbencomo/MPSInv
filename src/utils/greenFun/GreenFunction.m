classdef (Abstract) GreenFunction < handle
    % GreenFunction stores one multipole impulse response.
    %
    % One GreenFunction corresponds to one MultipoleComponent, one injected
    % solver, and one observed acoustic field. The solver owns the
    % computational/output grids, medium, time grid, and numerical options.
    % Values are stored on
    % the full lag grid -(K-1)*dt:dt:(K-1)*dt with size (2*K-1)-by-N,
    % where N is the number of output spatial points.

    properties (SetAccess = protected)
        component
        solver
        outputGrid
        outputField = ""
        times = zeros(1,0)
        lags = zeros(1,0)
        values
        isComputed (1,1) logical = false
    end

    properties (Dependent)
        dt
        numberOfTimes
        numberOfOutputChannels
        zeroLagIndex
    end

    methods
        function obj = GreenFunction(component,solver,outputField)
            arguments
                component (1,1) MultipoleComponent
                solver
                outputField (1,1) string
            end

            GreenFunction.validateSolverTimeGrid(solver);
            solverTimes = solver.outputTimes(:).';
            times = solverTimes(2:end);
            times(1) = 0;
            K = numel(times);
            N = prod(solver.outputGrid.N);

            obj.component = component;
            obj.solver = solver;
            obj.outputGrid = solver.outputGrid;
            obj.outputField = outputField;
            obj.times = times;
            obj.lags = (-(K-1):(K-1))*solver.dt;
            obj.values = zeros(2*K-1,N);
        end

        function value = get.dt(obj)
            value = obj.times(2)-obj.times(1);
        end

        function value = get.numberOfTimes(obj)
            value = numel(obj.times);
        end

        function value = get.numberOfOutputChannels(obj)
            value = prod(obj.outputGrid.N);
        end

        function value = get.zeroLagIndex(obj)
            value = obj.numberOfTimes;
        end

        function kernel = kernelAt(obj,outputIndex)
            % Return the kernel for one output spatial point.
            arguments
                obj
                outputIndex (1,1) double ...
                    {mustBeInteger,mustBePositive}
            end

            if outputIndex > obj.numberOfOutputChannels
                error('GreenFunction:OutputIndexOutOfRange', ...
                    'outputIndex exceeds the number of output channels.');
            end
            if ~obj.isComputed
                error('GreenFunction:NotComputed', ...
                    'Compute the Green function before accessing its kernel.');
            end

            kernel = obj.values(:,outputIndex);
        end

        function op = makeLinearOp(obj)
            % Construct the one-input multi-output convolution operator.
            if ~obj.isComputed
                error('GreenFunction:NotComputed', ...
                    ['Compute the Green function before constructing ' ...
                     'its linear operator.']);
            end

            kernels = reshape(obj.values, ...
                2*obj.numberOfTimes-1, ...
                obj.numberOfOutputChannels,1);
            op = ConvMultiOp(kernels,Scale=obj.dt);
        end
    end

    methods (Access = protected)
        function storeComputedValues(obj,newValues)
            expectedSize = [2*obj.numberOfTimes-1, ...
                obj.numberOfOutputChannels];
            if ~isnumeric(newValues) || ...
                    ~isequal(size(newValues),expectedSize) || ...
                    any(~isfinite(newValues),'all')
                error('GreenFunction:InvalidComputedValues', ...
                    ['Computed values must be a finite numeric array of ' ...
                     'size (2*K-1)-by-N.']);
            end

            obj.values = newValues;
            obj.isComputed = true;
        end

        function fullValues = extendSolverOutput(obj,timeSeries)
            % Embed output on [-dt,T] into the full lag grid [-T,T].
            K = obj.numberOfTimes;
            expectedNt = K+1;

            if timeSeries.Nt ~= expectedNt
                error('GreenFunction:UnexpectedSolverOutput', ...
                    'The solver output must contain K+1 time levels.');
            end

            shortValues = reshape(timeSeries.values,[],expectedNt).';
            if size(shortValues,2) ~= obj.numberOfOutputChannels
                error('GreenFunction:UnexpectedSolverOutput', ...
                    ['The solver output spatial size does not match the ' ...
                     'Green-function output grid.']);
            end

            fullValues = zeros(2*K-1,obj.numberOfOutputChannels, ...
                'like',shortValues);
            fullValues(K-1:end,:) = shortValues;
        end
    end

    methods (Static, Access = private)
        function validateSolverTimeGrid(solver)
            computationalTimes = solver.computationalTimes(:).';
            outputTimes = solver.outputTimes(:).';
            scale = max(1,max(abs(computationalTimes)));
            tolerance = 100*eps(scale);

            if numel(outputTimes) < 3
                error('GreenFunction:InvalidSolverTimes', ...
                    ['The solver must output at least the time levels ' ...
                     '[-dt,0,dt].']);
            end
            if numel(outputTimes) ~= numel(computationalTimes) || ...
                    any(abs(outputTimes-computationalTimes) > tolerance)
                error('GreenFunction:OutputTimesMustMatch', ...
                    ['Solver outputTimes must equal computationalTimes so ' ...
                     'that every computed time level is stored.']);
            end
            expectedStart = -solver.dt;
            if abs(outputTimes(1)-expectedStart) > tolerance || ...
                    abs(outputTimes(2)) > tolerance
                error('GreenFunction:InvalidSolverTimes', ...
                    ['The solver time grid must begin [-dt,0] and then ' ...
                     'continue with its uniform time step.']);
            end
        end
    end

    methods (Abstract)
        compute(obj,varargin)
        % Numerically compute and store the impulse response.
    end
end
