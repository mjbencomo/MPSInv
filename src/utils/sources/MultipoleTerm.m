classdef (Abstract) MultipoleTerm
    % MultipoleTerm is the common base class for one multipole term.

    properties (SetAccess = private)
        timeFunction (1,1) function_handle
        amplitude (1,1) double
    end

    properties (Abstract, SetAccess = private)
        location
        derivativeOrder
        targetField
    end

    methods
        function obj = MultipoleTerm(timeFunction,amplitude)
            arguments
                timeFunction (1,1) function_handle
                amplitude (1,1) double {mustBeFinite} = 1
            end

            obj.timeFunction = timeFunction;
            obj.amplitude = amplitude;
        end

        function value = evaluateTime(obj,t)
            % Evaluate the amplitude-scaled temporal factor.
            arguments
                obj
                t (1,1) double {mustBeFinite}
            end

            value = obj.timeFunction(t);
            if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value)
                error('MultipoleTerm:InvalidTimeOutput', ...
                    ['The time function must return one finite numeric ' ...
                     'scalar for a scalar time input.']);
            end
            value = obj.amplitude*value;
        end
    end

    methods (Abstract)
        source = discretize(obj,varargin)
        % Construct an AcousticSource containing only this term.

        addTo(obj,source,varargin)
        % Discretize this term and append it to an AcousticSource.
    end
end
