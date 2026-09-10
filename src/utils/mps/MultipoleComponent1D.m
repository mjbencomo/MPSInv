classdef MultipoleComponent1D < MultipoleComponent
    % MultipoleComponent1D describes D^s delta(x-xc) without its time
    % function.

    properties (SetAccess = private)
        location (1,1) double = 0
        derivativeOrder (1,1) double = 0
        targetField (1,1) string = "pressure"
    end

    methods
        function obj = MultipoleComponent1D( ...
                location,derivativeOrder,options)
            arguments
                location (1,1) double {mustBeFinite}
                derivativeOrder (1,1) double ...
                    {mustBeInteger,mustBeNonnegative}
                options.TargetField (1,1) string ...
                    {mustBeMember(options.TargetField, ...
                    ["pressure","velocity"])} = "pressure"
                options.ApproximationOrder (1,1) double ...
                    {mustBeInteger,mustBePositive} = 4
                options.Label (1,1) string = ""
            end

            obj@MultipoleComponent( ...
                options.ApproximationOrder,options.Label);
            obj.location = location;
            obj.derivativeOrder = derivativeOrder;
            obj.targetField = options.TargetField;
        end

        function term = withTimeFunction(obj,timeFunction,options)
            % Construct a MultipoleTerm1D from this spatial component.
            arguments
                obj
                timeFunction (1,1) function_handle
                options.Amplitude (1,1) double {mustBeFinite} = 1
            end

            term = MultipoleTerm1D( ...
                obj.location,obj.derivativeOrder,timeFunction, ...
                Amplitude=options.Amplitude, ...
                TargetField=obj.targetField);
        end
    end
end
