classdef MultipoleTerm2D
    % MultipoleTerm2D represents a(t) Dx^sx Dy^sy delta(x-xc,y-yc).

    properties (SetAccess = private)
        location (1,2) double = [0,0]
        derivativeOrder (1,2) double = [0,0]
        timeFunction (1,1) function_handle = @(t) ones(size(t))
        amplitude (1,1) double = 1
        targetField (1,1) string = "pressure"
    end

    methods
        function obj = MultipoleTerm2D( ...
                location,derivativeOrder,timeFunction,options)
            arguments
                location (1,2) double {mustBeFinite}
                derivativeOrder (1,2) double ...
                    {mustBeInteger,mustBeNonnegative}
                timeFunction (1,1) function_handle
                options.Amplitude (1,1) double {mustBeFinite} = 1
                options.TargetField (1,1) string ...
                    {mustBeMember(options.TargetField, ...
                    ["pressure","velocityX","velocityY"])} = "pressure"
            end

            obj.location = location;
            obj.derivativeOrder = derivativeOrder;
            obj.timeFunction = timeFunction;
            obj.amplitude = options.Amplitude;
            obj.targetField = options.TargetField;
        end

        function value = evaluateTime(obj,t)
            % Evaluate the amplitude-scaled temporal factor.
            arguments
                obj
                t (1,1) double {mustBeFinite}
            end

            value = obj.timeFunction(t);
            if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value)
                error('MultipoleTerm2D:InvalidTimeOutput', ...
                    ['The time function must return one finite numeric ' ...
                     'scalar for a scalar time input.']);
            end
            value = obj.amplitude*value;
        end
    end
end
