classdef MultipoleComponent1D < MultipoleComponent
    % MultipoleComponent1D describes D^s delta(x-xc).

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
            end

            obj.location = location;
            obj.derivativeOrder = derivativeOrder;
            obj.targetField = options.TargetField;
        end

        function text = description(obj)
            % Return a readable description generated from component data.
            text = sprintf( ...
                '%s multipole, derivative order %d, at x = %.6g', ...
                obj.targetField,obj.derivativeOrder,obj.location);
            text = string(text);
        end

        function [indices,weights,values] = createStencil( ...
                obj,grid,approximationOrder)
            % Discretize the component on a one-dimensional grid.
            arguments
                obj
                grid (1,1) GridSpace1D
                approximationOrder (1,1) double ...
                    {mustBeInteger,mustBePositive}
            end

            [values,indices,weights] = MPSappx( ...
                grid,obj.location,approximationOrder, ...
                obj.derivativeOrder);
        end
    end
end
