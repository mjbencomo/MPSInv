classdef MultipoleSource1D
    % MultipoleSource1D discretizes sums of 1D multipole terms.
    %
    % The result of discretize is an AcousticSource1D containing localized
    % pressure and velocity terms on their respective staggered grids.

    properties (SetAccess = private)
        terms
        approximationOrder (1,1) double
    end

    methods
        function obj = MultipoleSource1D(terms,options)
            arguments
                terms
                options.ApproximationOrder (1,1) double ...
                    {mustBeInteger,mustBePositive} = 4
            end

            if isempty(terms)
                terms = MultipoleTerm1D.empty(1,0);
            elseif iscell(terms)
                if any(~cellfun(@(term) isa(term,'MultipoleTerm1D'),terms))
                    error('MultipoleSource1D:InvalidTerms', ...
                        'Every term must be a MultipoleTerm1D object.');
                end
                terms = [terms{:}];
            elseif ~isa(terms,'MultipoleTerm1D')
                error('MultipoleSource1D:InvalidTerms', ...
                    'terms must contain MultipoleTerm1D objects.');
            end

            obj.terms = reshape(terms,1,[]);
            obj.approximationOrder = options.ApproximationOrder;
        end

        function source = discretize(obj,pressureGrid,velocityGrid)
            % Discretize all terms on the appropriate staggered grids.
            arguments
                obj
                pressureGrid (1,1) GridSpace1D
                velocityGrid (1,1) GridSpace1D
            end

            source = AcousticSource1D.zero();
            q = obj.approximationOrder;

            for k = 1:numel(obj.terms)
                term = obj.terms(k);
                if term.targetField == "pressure"
                    [~,indices,weights] = MPSappx(pressureGrid, ...
                        term.location,q,term.derivativeOrder);

                    % The non-PML solver overwrites pressure at its outer
                    % endpoints, so require pressure support to be interior.
                    if any(indices == 1 | indices == pressureGrid.N)
                        error('MultipoleSource1D:PressureStencilAtBoundary', ...
                            ['The pressure multipole stencil must exclude ' ...
                             'the pressure-grid endpoints.']);
                    end

                    source.addLocalizedPressureTerm( ...
                        pressureGrid,indices,weights, ...
                        @(t) term.evaluateTime(t));
                else
                    [~,indices,weights] = MPSappx(velocityGrid, ...
                        term.location,q,term.derivativeOrder);
                    source.addLocalizedVelocityTerm( ...
                        velocityGrid,indices,weights, ...
                        @(t) term.evaluateTime(t));
                end
            end
        end

        function value = numberOfTerms(obj)
            value = numel(obj.terms);
        end
    end
end
