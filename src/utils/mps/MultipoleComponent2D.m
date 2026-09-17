classdef MultipoleComponent2D < MultipoleComponent
    % MultipoleComponent2D describes
    % Dx^sx Dy^sy delta(x-xc,y-yc).

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
            end

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

        function [indices,weights,xIndices,yIndices] = createStencil( ...
                obj,grid,approximationOrder)
            % Discretize the component using a tensor-product stencil.
            arguments
                obj
                grid (1,1) GridSpace2D
                approximationOrder (1,1) double ...
                    {mustBeInteger,mustBePositive}
            end

            [~,xIndices,xWeights] = MPSappx( ...
                grid.x.pts,obj.location(1),approximationOrder, ...
                obj.derivativeOrder(1));
            [~,yIndices,yWeights] = MPSappx( ...
                grid.y.pts,obj.location(2),approximationOrder, ...
                obj.derivativeOrder(2));

            [IX,IY] = ndgrid(xIndices,yIndices);
            weights = xWeights(:)*yWeights(:).';
            indices = sub2ind(grid.N,IX,IY);

            indices = indices(:);
            weights = weights(:);
        end
    end
end
