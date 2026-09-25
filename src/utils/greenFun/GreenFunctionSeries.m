classdef GreenFunctionSeries < handle
    % GreenFunctionSeries stores compatible component Green functions.
    %
    % If the series contains M Green functions sampled at N output points
    % and K convolution times, values has size (2*K-1)-by-N-by-M. The third
    % dimension follows the component/input-channel ordering.

    properties (SetAccess = private)
        greenFunctions (1,:) cell
    end

    properties (Dependent)
        components
        solver
        outputGrid
        outputField
        times
        lags
        values
        dt
        numberOfTimes
        numberOfInputChannels
        numberOfOutputChannels
        zeroLagIndex
        isComputed
    end

    methods
        function obj = GreenFunctionSeries(greenFunctions)
            greenFunctions = GreenFunctionSeries.normalizeGreenFunctions( ...
                greenFunctions);
            GreenFunctionSeries.validateCompatibility(greenFunctions);
            obj.greenFunctions = greenFunctions;
        end

        function value = get.components(obj)
            value = cellfun(@(G) G.component,obj.greenFunctions, ...
                UniformOutput=false);
        end

        function value = get.solver(obj)
            value = obj.greenFunctions{1}.solver;
        end

        function value = get.outputGrid(obj)
            value = obj.greenFunctions{1}.outputGrid;
        end

        function value = get.outputField(obj)
            value = obj.greenFunctions{1}.outputField;
        end

        function value = get.times(obj)
            value = obj.greenFunctions{1}.times;
        end

        function value = get.lags(obj)
            value = obj.greenFunctions{1}.lags;
        end

        function value = get.values(obj)
            value = zeros(2*obj.numberOfTimes-1, ...
                obj.numberOfOutputChannels,obj.numberOfInputChannels, ...
                'like',obj.greenFunctions{1}.values);
            for m = 1:obj.numberOfInputChannels
                value(:,:,m) = obj.greenFunctions{m}.values;
            end
        end

        function value = get.dt(obj)
            value = obj.greenFunctions{1}.dt;
        end

        function value = get.numberOfTimes(obj)
            value = obj.greenFunctions{1}.numberOfTimes;
        end

        function value = get.numberOfInputChannels(obj)
            value = numel(obj.greenFunctions);
        end

        function value = get.numberOfOutputChannels(obj)
            value = obj.greenFunctions{1}.numberOfOutputChannels;
        end

        function value = get.zeroLagIndex(obj)
            value = obj.greenFunctions{1}.zeroLagIndex;
        end

        function value = get.isComputed(obj)
            value = all(cellfun(@(G) G.isComputed,obj.greenFunctions));
        end

        function compute(obj,approximationOrders)
            % Compute every component Green function.
            arguments
                obj
                approximationOrders (1,:) double ...
                    {mustBeInteger,mustBePositive}
            end

            orders = obj.normalizeApproximationOrders(approximationOrders);
            for m = 1:obj.numberOfInputChannels
                obj.greenFunctions{m}.compute(orders(m));
            end
        end

        function G = greenFunctionAt(obj,inputIndex)
            % Return the GreenFunction for one input channel.
            obj.validateInputIndex(inputIndex);
            G = obj.greenFunctions{inputIndex};
        end

        function kernel = kernelAt(obj,outputIndex,inputIndex)
            % Return G_m(x_n,:) for output n and input m.
            obj.validateOutputIndex(outputIndex);
            obj.validateInputIndex(inputIndex);
            if ~obj.greenFunctions{inputIndex}.isComputed
                error('GreenFunctionSeries:NotComputed', ...
                    'Compute the requested Green function first.');
            end
            kernel = obj.greenFunctions{inputIndex}.kernelAt(outputIndex);
        end

        function op = makeLinearOp(obj)
            % Create the M-input, N-output convolution operator.
            if ~obj.isComputed
                error('GreenFunctionSeries:NotComputed', ...
                    ['Compute every Green function before constructing ' ...
                     'the linear operator.']);
            end
            op = ConvMultiOp(obj.values,Scale=obj.dt);
        end
    end

    methods (Access = private)
        function orders = normalizeApproximationOrders( ...
                obj,approximationOrders)
            if isscalar(approximationOrders)
                orders = repmat(approximationOrders,1, ...
                    obj.numberOfInputChannels);
            elseif numel(approximationOrders) == ...
                    obj.numberOfInputChannels
                orders = reshape(approximationOrders,1,[]);
            else
                error( ...
                    'GreenFunctionSeries:ApproximationOrderSizeMismatch', ...
                    ['approximationOrders must be scalar or contain one ' ...
                     'entry per input channel.']);
            end
        end

        function validateInputIndex(obj,index)
            arguments
                obj
                index (1,1) double {mustBeInteger,mustBePositive}
            end
            if index > obj.numberOfInputChannels
                error('GreenFunctionSeries:InputIndexOutOfRange', ...
                    'inputIndex exceeds the number of input channels.');
            end
        end

        function validateOutputIndex(obj,index)
            arguments
                obj
                index (1,1) double {mustBeInteger,mustBePositive}
            end
            if index > obj.numberOfOutputChannels
                error('GreenFunctionSeries:OutputIndexOutOfRange', ...
                    'outputIndex exceeds the number of output channels.');
            end
        end
    end

    methods (Static, Access = private)
        function greenFunctions = normalizeGreenFunctions(greenFunctions)
            if isempty(greenFunctions)
                error('GreenFunctionSeries:EmptySeries', ...
                    'A Green-function series must contain at least one entry.');
            elseif iscell(greenFunctions)
                if any(~cellfun(@(G) isa(G,'GreenFunction'),greenFunctions))
                    error('GreenFunctionSeries:InvalidGreenFunctions', ...
                        'Every entry must inherit from GreenFunction.');
                end
                greenFunctions = reshape(greenFunctions,1,[]);
            elseif isa(greenFunctions,'GreenFunction')
                greenFunctions = num2cell(reshape(greenFunctions,1,[]));
            else
                error('GreenFunctionSeries:InvalidGreenFunctions', ...
                    'greenFunctions must contain GreenFunction objects.');
            end
        end

        function validateCompatibility(greenFunctions)
            reference = greenFunctions{1};
            for m = 2:numel(greenFunctions)
                candidate = greenFunctions{m};
                if ~isequal(candidate.solver,reference.solver)
                    error('GreenFunctionSeries:SolverMismatch', ...
                        'All Green functions must use the same solver object.');
                end
                if candidate.outputField ~= reference.outputField
                    error('GreenFunctionSeries:OutputFieldMismatch', ...
                        'All Green functions must use the same output field.');
                end
                if ~isequal(candidate.times,reference.times) || ...
                        ~isequal(candidate.lags,reference.lags)
                    error('GreenFunctionSeries:TimeGridMismatch', ...
                        'All Green functions must use the same time grid.');
                end
            end
        end
    end
end
