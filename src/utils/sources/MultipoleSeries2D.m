classdef MultipoleSeries2D < MultipoleSeries
    % MultipoleSeries2D discretizes a finite collection of 2D terms.

    properties (Dependent)
        locations
    end

    methods
        function obj = MultipoleSeries2D(terms,options)
            arguments
                terms
                options.ApproximationOrders (1,:) double ...
                    {mustBeInteger,mustBePositive} = 4
            end

            terms = MultipoleSeries2D.normalizeTerms(terms);
            obj@MultipoleSeries(terms,options.ApproximationOrders);
        end

        function value = get.locations(obj)
            if isempty(obj.terms)
                value = zeros(0,2);
            else
                value = vertcat(obj.terms.location);
            end
        end

        function source = discretize( ...
                obj,pressureGrid,velocityXGrid,velocityYGrid)
            arguments
                obj
                pressureGrid (1,1) GridSpace2D
                velocityXGrid (1,1) GridSpace2D
                velocityYGrid (1,1) GridSpace2D
            end

            source = AcousticSource2D.zero();
            for k = 1:obj.numberOfTerms
                obj.terms(k).addTo( ...
                    source,pressureGrid,velocityXGrid,velocityYGrid, ...
                    obj.approximationOrders(k));
            end
        end
    end

    methods (Static, Access = private)
        function terms = normalizeTerms(terms)
            if isempty(terms)
                terms = MultipoleTerm2D.empty(1,0);
            elseif iscell(terms)
                if any(~cellfun(@(term) isa(term,'MultipoleTerm2D'),terms))
                    error('MultipoleSeries2D:InvalidTerms', ...
                        'Every term must be a MultipoleTerm2D object.');
                end
                terms = [terms{:}];
            elseif ~isa(terms,'MultipoleTerm2D')
                error('MultipoleSeries2D:InvalidTerms', ...
                    'terms must contain MultipoleTerm2D objects.');
            end
            terms = reshape(terms,1,[]);
        end
    end
end
