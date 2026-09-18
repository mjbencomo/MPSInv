classdef MultipoleComponent2D < MultipoleComponent
    % MultipoleComponent2D describes
<<<<<<< HEAD
    % Dx^sx Dy^sy delta(x-xc,y-yc) without its time function.
=======
    % Dx^sx Dy^sy delta(x-xc,y-yc).
>>>>>>> b0d5b1dcea97a019fdcf42ee218d06a571246d9f

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
<<<<<<< HEAD
                options.ApproximationOrder (1,1) double ...
                    {mustBeInteger,mustBePositive} = 4
            end

            obj@MultipoleComponent(options.ApproximationOrder);
=======
            end

>>>>>>> b0d5b1dcea97a019fdcf42ee218d06a571246d9f
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

<<<<<<< HEAD
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
=======
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
>>>>>>> b0d5b1dcea97a019fdcf42ee218d06a571246d9f
        end
    end
end
