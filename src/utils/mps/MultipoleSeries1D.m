classdef MultipoleSeries1D < MultipoleSeries
    % MultipoleSeries1D discretizes a finite collection of 1D terms.

    properties (Dependent)
        locations
    end

    methods
        function obj = MultipoleSeries1D(terms,options)
            arguments
                terms
                options.ApproximationOrders (1,:) double ...
                    {mustBeInteger,mustBePositive} = 4
            end

            terms = MultipoleSeries1D.normalizeTerms(terms);
            obj@MultipoleSeries(terms,options.ApproximationOrders);
        end

        function value = get.locations(obj)
            if isempty(obj.terms)
                value = zeros(0,1);
            else
                value = reshape([obj.terms.location],[],1);
            end
        end

        function source = discretize(obj,pressureGrid,velocityGrid)
            arguments
                obj
                pressureGrid (1,1) GridSpace1D
                velocityGrid (1,1) GridSpace1D
            end

            source = AcousticSource1D.zero();
            for k = 1:obj.numberOfTerms
                obj.terms(k).addTo( ...
                    source,pressureGrid,velocityGrid, ...
                    obj.approximationOrders(k));
            end
        end
    end

    methods (Static, Access = private)
        function terms = normalizeTerms(terms)
            if isempty(terms)
                terms = MultipoleTerm1D.empty(1,0);
            elseif iscell(terms)
                if any(~cellfun(@(term) isa(term,'MultipoleTerm1D'),terms))
                    error('MultipoleSeries1D:InvalidTerms', ...
                        'Every term must be a MultipoleTerm1D object.');
                end
                terms = [terms{:}];
            elseif ~isa(terms,'MultipoleTerm1D')
                error('MultipoleSeries1D:InvalidTerms', ...
                    'terms must contain MultipoleTerm1D objects.');
            end
            terms = reshape(terms,1,[]);
        end
    end
end
