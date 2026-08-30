classdef MultipoleTerm1D < MultipoleTerm
    % MultipoleTerm1D represents a(t) D^s delta(x-xc).

    properties (SetAccess = private)
        location (1,1) double = 0
        derivativeOrder (1,1) double = 0
        targetField (1,1) string = "pressure"
    end

    methods
        function obj = MultipoleTerm1D( ...
                location,derivativeOrder,timeFunction,options)
            arguments
                location (1,1) double {mustBeFinite}
                derivativeOrder (1,1) double ...
                    {mustBeInteger,mustBeNonnegative}
                timeFunction (1,1) function_handle
                options.Amplitude (1,1) double {mustBeFinite} = 1
                options.TargetField (1,1) string ...
                    {mustBeMember(options.TargetField, ...
                    ["pressure","velocity"])} = "pressure"
            end

            obj@MultipoleTerm(timeFunction,options.Amplitude);
            obj.location = location;
            obj.derivativeOrder = derivativeOrder;
            obj.targetField = options.TargetField;
        end

        function source = discretize( ...
                obj,pressureGrid,velocityGrid,options)
            % Construct an AcousticSource1D containing this term.
            arguments
                obj
                pressureGrid (1,1) GridSpace1D
                velocityGrid (1,1) GridSpace1D
                options.ApproximationOrder (1,1) double ...
                    {mustBeInteger,mustBePositive} = 4
            end

            source = AcousticSource1D.zero();
            obj.addTo(source,pressureGrid,velocityGrid, ...
                options.ApproximationOrder);
        end

        function addTo(obj,source,pressureGrid,velocityGrid,q)
            % Discretize this term and append it to source.
            arguments
                obj
                source (1,1) AcousticSource1D
                pressureGrid (1,1) GridSpace1D
                velocityGrid (1,1) GridSpace1D
                q (1,1) double {mustBeInteger,mustBePositive}
            end

            timeFactor = @(t) obj.evaluateTime(t);
            switch obj.targetField
                case "pressure"
                    [~,indices,weights] = MPSappx(pressureGrid, ...
                        obj.location,q,obj.derivativeOrder);

                    if any(indices == 1 | indices == pressureGrid.N)
                        error( ...
                            'MultipoleTerm1D:PressureStencilAtBoundary', ...
                            ['The pressure multipole stencil must exclude ' ...
                             'the pressure-grid endpoints.']);
                    end

                    source.addLocalizedPressureTerm( ...
                        pressureGrid,indices,weights,timeFactor);

                case "velocity"
                    [~,indices,weights] = MPSappx(velocityGrid, ...
                        obj.location,q,obj.derivativeOrder);
                    source.addLocalizedVelocityTerm( ...
                        velocityGrid,indices,weights,timeFactor);
            end
        end
    end
end
