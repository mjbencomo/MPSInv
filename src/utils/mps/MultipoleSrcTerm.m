classdef (Abstract) MultipoleSrcTerm
    % MultipoleSrcTerm combines one spatial multipole component with its
    % time-dependent coefficient.

    properties (SetAccess = private)
        component 
        timeFunction
    end

    properties (Dependent)
        location
        derivativeOrder
        targetField
    end

    methods
        function obj = MultipoleSrcTerm(component,timeFunction)
            arguments
                component (1,1) MultipoleComponent
                timeFunction (1,1) function_handle
            end

            obj.component = component;
            obj.timeFunction = timeFunction;
        end

        function value = get.location(obj)
            value = obj.component.location;
        end

        function value = get.derivativeOrder(obj)
            value = obj.component.derivativeOrder;
        end

        function value = get.targetField(obj)
            value = obj.component.targetField;
        end

        function value = evaluateTime(obj,t)
            % Evaluate the complete temporal coefficient w(t).
            arguments
                obj
                t (1,1) double {mustBeFinite}
            end

            value = obj.timeFunction(t);
            if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value)
                error('MultipoleSrcTerm:InvalidTimeOutput', ...
                    ['The time function must return one finite numeric ' ...
                     'scalar for a scalar time input.']);
            end
        end
    end

    methods (Abstract)
        source = discretize(obj,varargin)
        % Construct an AcousticSource containing only this term.

        addTo(obj,source,varargin)
        % Discretize this term and append it to an AcousticSource.
    end
end
