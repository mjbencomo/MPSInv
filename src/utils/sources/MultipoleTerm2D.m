classdef MultipoleTerm2D < MultipoleTerm
    % MultipoleTerm2D represents a(t) Dx^sx Dy^sy delta(x-xc,y-yc).

    properties (SetAccess = private)
        location (1,2) double = [0,0]
        derivativeOrder (1,2) double = [0,0]
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

            obj@MultipoleTerm(timeFunction,options.Amplitude);
            obj.location = location;
            obj.derivativeOrder = derivativeOrder;
            obj.targetField = options.TargetField;
        end

        function source = discretize( ...
                obj,pressureGrid,velocityXGrid,velocityYGrid,options)
            % Construct an AcousticSource2D containing this term.
            arguments
                obj
                pressureGrid (1,1) GridSpace2D
                velocityXGrid (1,1) GridSpace2D
                velocityYGrid (1,1) GridSpace2D
                options.ApproximationOrder (1,1) double ...
                    {mustBeInteger,mustBePositive} = 4
            end

            source = AcousticSource2D.zero();
            obj.addTo(source,pressureGrid,velocityXGrid,velocityYGrid, ...
                options.ApproximationOrder);
        end

        function addTo( ...
                obj,source,pressureGrid,velocityXGrid,velocityYGrid,q)
            % Discretize this term and append it to source.
            arguments
                obj
                source (1,1) AcousticSource2D
                pressureGrid (1,1) GridSpace2D
                velocityXGrid (1,1) GridSpace2D
                velocityYGrid (1,1) GridSpace2D
                q (1,1) double {mustBeInteger,mustBePositive}
            end

            timeFactor = @(t) obj.evaluateTime(t);
            switch obj.targetField
                case "pressure"
                    [indices,weights,xIndices,yIndices] = ...
                        obj.createStencil(pressureGrid,q);

                    if any(xIndices == 1 | ...
                            xIndices == pressureGrid.N(1)) || ...
                            any(yIndices == 1 | ...
                            yIndices == pressureGrid.N(2))
                        error( ...
                            'MultipoleTerm2D:PressureStencilAtBoundary', ...
                            ['The pressure multipole stencil must exclude ' ...
                             'all pressure-grid edges.']);
                    end

                    source.addLocalizedPressureTerm( ...
                        pressureGrid,indices,weights,timeFactor);

                case "velocityX"
                    [indices,weights] = ...
                        obj.createStencil(velocityXGrid,q);
                    source.addLocalizedVelocityXTerm( ...
                        velocityXGrid,indices,weights,timeFactor);

                case "velocityY"
                    [indices,weights] = ...
                        obj.createStencil(velocityYGrid,q);
                    source.addLocalizedVelocityYTerm( ...
                        velocityYGrid,indices,weights,timeFactor);
            end
        end
    end

    methods (Access = private)
        function [indices,weights,xIndices,yIndices] = ...
                createStencil(obj,grid,q)
            [~,xIndices,xWeights] = MPSappx( ...
                grid.x.pts,obj.location(1),q,obj.derivativeOrder(1));
            [~,yIndices,yWeights] = MPSappx( ...
                grid.y.pts,obj.location(2),q,obj.derivativeOrder(2));

            [IX,IY] = ndgrid(xIndices,yIndices);
            weights = xWeights*yWeights.';
            indices = sub2ind(grid.N,IX,IY);
            indices = indices(:);
            weights = weights(:);
        end
    end
end
