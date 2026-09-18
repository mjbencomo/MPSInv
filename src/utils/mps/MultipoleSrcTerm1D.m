classdef MultipoleSrcTerm1D < MultipoleSrcTerm
    % MultipoleSrcTerm1D represents w(t) D^s delta(x-xc).

    methods
        function obj = MultipoleSrcTerm1D(component,timeFunction)
            arguments
                component (1,1) MultipoleComponent1D
                timeFunction (1,1) function_handle
            end

            obj@MultipoleSrcTerm(component,timeFunction);
        end

        function source = discretize( ...
                obj,pressureGrid,velocityGrid,approximationOrder)
            % Construct an AcousticSource1D containing this term.
            arguments
                obj
                pressureGrid (1,1) GridSpace1D
                velocityGrid (1,1) GridSpace1D
                approximationOrder (1,1) double ...
                    {mustBeInteger,mustBePositive}
            end

            source = AcousticSource1D.zero();
            obj.addTo(source,pressureGrid,velocityGrid, ...
                approximationOrder);
        end

        function addTo( ...
                obj,source,pressureGrid,velocityGrid,approximationOrder)
            % Discretize this term and append it to source.
            arguments
                obj
                source (1,1) AcousticSource1D
                pressureGrid (1,1) GridSpace1D
                velocityGrid (1,1) GridSpace1D
                approximationOrder (1,1) double ...
                    {mustBeInteger,mustBePositive}
            end

            timeFactor = @(t) obj.evaluateTime(t);
            switch obj.targetField
                case "pressure"
                    [indices,weights] = obj.component.createStencil( ...
                        pressureGrid,approximationOrder);

                    if any(indices == 1 | indices == pressureGrid.N)
                        error( ...
                            'MultipoleSrcTerm1D:PressureStencilAtBoundary', ...
                            ['The pressure multipole stencil must exclude ' ...
                             'the pressure-grid endpoints.']);
                    end

                    source.addLocalizedPressureTerm( ...
                        pressureGrid,indices,weights,timeFactor);

                case "velocity"
                    [indices,weights] = obj.component.createStencil( ...
                        velocityGrid,approximationOrder);
                    source.addLocalizedVelocityTerm( ...
                        velocityGrid,indices,weights,timeFactor);
            end
        end
    end
end
