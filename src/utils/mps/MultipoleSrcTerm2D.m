classdef MultipoleSrcTerm2D < MultipoleSrcTerm
    % MultipoleSrcTerm2D represents
    % w(t) Dx^sx Dy^sy delta(x-xc,y-yc).

    methods
        function obj = MultipoleSrcTerm2D(component,timeFunction)
            arguments
                component (1,1) MultipoleComponent2D
                timeFunction (1,1) function_handle
            end

            obj@MultipoleSrcTerm(component,timeFunction);
        end

        function source = discretize(obj,pressureGrid,velocityXGrid, ...
                velocityYGrid,approximationOrder)
            % Construct an AcousticSource2D containing this term.
            arguments
                obj
                pressureGrid (1,1) GridSpace2D
                velocityXGrid (1,1) GridSpace2D
                velocityYGrid (1,1) GridSpace2D
                approximationOrder (1,1) double ...
                    {mustBeInteger,mustBePositive}
            end

            source = AcousticSource2D.zero();
            obj.addTo(source,pressureGrid,velocityXGrid,velocityYGrid, ...
                approximationOrder);
        end

        function addTo(obj,source,pressureGrid,velocityXGrid, ...
                velocityYGrid,approximationOrder)
            % Discretize this term and append it to source.
            arguments
                obj
                source (1,1) AcousticSource2D
                pressureGrid (1,1) GridSpace2D
                velocityXGrid (1,1) GridSpace2D
                velocityYGrid (1,1) GridSpace2D
                approximationOrder (1,1) double ...
                    {mustBeInteger,mustBePositive}
            end

            timeFactor = @(t) obj.evaluateTime(t);
            switch obj.targetField
                case "pressure"
                    [indices,weights,xIndices,yIndices] = ...
                        obj.component.createStencil( ...
                        pressureGrid,approximationOrder);

                    if any(xIndices == 1 | ...
                            xIndices == pressureGrid.N(1)) || ...
                            any(yIndices == 1 | ...
                            yIndices == pressureGrid.N(2))
                        error( ...
                            'MultipoleSrcTerm2D:PressureStencilAtBoundary', ...
                            ['The pressure multipole stencil must exclude ' ...
                             'all pressure-grid edges.']);
                    end

                    source.addLocalizedPressureTerm( ...
                        pressureGrid,indices,weights,timeFactor);

                case "velocityX"
                    [indices,weights] = obj.component.createStencil( ...
                        velocityXGrid,approximationOrder);
                    source.addLocalizedVelocityXTerm( ...
                        velocityXGrid,indices,weights,timeFactor);

                case "velocityY"
                    [indices,weights] = obj.component.createStencil( ...
                        velocityYGrid,approximationOrder);
                    source.addLocalizedVelocityYTerm( ...
                        velocityYGrid,indices,weights,timeFactor);
            end
        end
    end
end
