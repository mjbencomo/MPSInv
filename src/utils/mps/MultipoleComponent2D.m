classdef MultipoleComponent2D < MultipoleComponent
    % MultipoleComponent2D describes
    % Dx^sx Dy^sy delta(x-xc,y-yc) without its time function.

    properties (SetAccess = private)
        location (1,2) double = [0,0]
        derivativeOrder (1,2) double = [0,0]
        targetField (1,1) string = "pressure"
    end

    methods
        function obj = MultipoleComponent2D( ...
                location,derivativeOrder,options)
            arguments
                location (1,2) double {mustBeFinite}
                derivativeOrder (1,2) double ...
                    {mustBeInteger,mustBeNonnegative}
                options.TargetField (1,1) string ...
                    {mustBeMember(options.TargetField, ...
                    ["pressure","velocityX","velocityY"])} = "pressure"
                options.ApproximationOrder (1,1) double ...
                    {mustBeInteger,mustBePositive} = 4
            end

            obj@MultipoleComponent(options.ApproximationOrder);
            obj.location = location;
            obj.derivativeOrder = derivativeOrder;
            obj.targetField = options.TargetField;
        end

        function text = description(obj)
            % Return a readable description generated from component data.
            text = sprintf([ ...
                '%s multipole, derivative order [%d,%d], ' ...
                'at [x,y] = [%.6g,%.6g]'], ...
                obj.targetField, ...
                obj.derivativeOrder(1),obj.derivativeOrder(2), ...
                obj.location(1),obj.location(2));
            text = string(text);
        end

        function term = withTimeFunction(obj,timeFunction,options)
            % Construct a MultipoleTerm2D from this spatial component.
            arguments
                obj
                timeFunction (1,1) function_handle
                options.Amplitude (1,1) double {mustBeFinite} = 1
            end

            term = MultipoleTerm2D( ...
                obj.location,obj.derivativeOrder,timeFunction, ...
                Amplitude=options.Amplitude, ...
                TargetField=obj.targetField);
        end
    end
end
