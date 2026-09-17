classdef MultipoleSrcSeries2D < MultipoleSrcSeries
    % MultipoleSrcSeries2D represents a finite collection of 2D terms.

    properties (Dependent)
        components
        locations
    end

    methods
        function obj = MultipoleSrcSeries2D(terms)
            terms = MultipoleSrcSeries2D.normalizeTerms(terms);
            obj@MultipoleSrcSeries(terms);
        end

        function value = get.components(obj)
            if isempty(obj.terms)
                value = MultipoleComponent2D.empty(1,0);
            else
                value = [obj.terms.component];
            end
        end

        function value = get.locations(obj)
            if isempty(obj.terms)
                value = zeros(0,2);
            else
                value = vertcat(obj.terms.location);
            end
        end

        function source = discretize(obj,pressureGrid,velocityXGrid, ...
                velocityYGrid,approximationOrders)
            arguments
                obj
                pressureGrid (1,1) GridSpace2D
                velocityXGrid (1,1) GridSpace2D
                velocityYGrid (1,1) GridSpace2D
                approximationOrders (1,:) double ...
                    {mustBeInteger,mustBePositive}
            end

            orders = obj.normalizeApproximationOrders( ...
                approximationOrders);
            source = AcousticSource2D.zero();
            for k = 1:obj.numberOfTerms
                obj.terms(k).addTo(source,pressureGrid,velocityXGrid, ...
                    velocityYGrid,orders(k));
            end
        end
    end

    methods (Static, Access = private)
        function terms = normalizeTerms(terms)
            if isempty(terms)
                terms = MultipoleSrcTerm2D.empty(1,0);
            elseif iscell(terms)
                if any(~cellfun( ...
                        @(term) isa(term,'MultipoleSrcTerm2D'),terms))
                    error('MultipoleSrcSeries2D:InvalidTerms', ...
                        ['Every term must be a ' ...
                         'MultipoleSrcTerm2D object.']);
                end
                terms = [terms{:}];
            elseif ~isa(terms,'MultipoleSrcTerm2D')
                error('MultipoleSrcSeries2D:InvalidTerms', ...
                    'terms must contain MultipoleSrcTerm2D objects.');
            end
            terms = reshape(terms,1,[]);
        end
    end
end
