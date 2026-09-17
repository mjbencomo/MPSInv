classdef MultipoleSrcSeries1D < MultipoleSrcSeries
    % MultipoleSrcSeries1D represents a finite collection of 1D terms.

    properties (Dependent)
        components
        locations
    end

    methods
        function obj = MultipoleSrcSeries1D(terms)
            terms = MultipoleSrcSeries1D.normalizeTerms(terms);
            obj@MultipoleSrcSeries(terms);
        end

        function value = get.components(obj)
            if isempty(obj.terms)
                value = MultipoleComponent1D.empty(1,0);
            else
                value = [obj.terms.component];
            end
        end

        function value = get.locations(obj)
            if isempty(obj.terms)
                value = zeros(0,1);
            else
                value = reshape([obj.terms.location],[],1);
            end
        end

        function source = discretize(obj,pressureGrid,velocityGrid, ...
                approximationOrders)
            arguments
                obj
                pressureGrid (1,1) GridSpace1D
                velocityGrid (1,1) GridSpace1D
                approximationOrders (1,:) double ...
                    {mustBeInteger,mustBePositive}
            end

            orders = obj.normalizeApproximationOrders( ...
                approximationOrders);
            source = AcousticSource1D.zero();
            for k = 1:obj.numberOfTerms
                obj.terms(k).addTo(source,pressureGrid,velocityGrid, ...
                    orders(k));
            end
        end
    end

    methods (Static, Access = private)
        function terms = normalizeTerms(terms)
            if isempty(terms)
                terms = MultipoleSrcTerm1D.empty(1,0);
            elseif iscell(terms)
                if any(~cellfun( ...
                        @(term) isa(term,'MultipoleSrcTerm1D'),terms))
                    error('MultipoleSrcSeries1D:InvalidTerms', ...
                        ['Every term must be a ' ...
                         'MultipoleSrcTerm1D object.']);
                end
                terms = [terms{:}];
            elseif ~isa(terms,'MultipoleSrcTerm1D')
                error('MultipoleSrcSeries1D:InvalidTerms', ...
                    'terms must contain MultipoleSrcTerm1D objects.');
            end
            terms = reshape(terms,1,[]);
        end
    end
end
