classdef (Abstract) MultipoleSrcSeries
    % MultipoleSrcSeries is the common base for finite multipole sources.

    properties (SetAccess = private)
        terms
    end

    properties (Dependent)
        numberOfTerms
    end

    properties (Abstract, Dependent)
        components
        locations
    end

    methods
        function obj = MultipoleSrcSeries(terms)
            if ~isempty(terms) && ~isa(terms,'MultipoleSrcTerm')
                error('MultipoleSrcSeries:InvalidTerms', ...
                    'Every term must inherit from MultipoleSrcTerm.');
            end

            obj.terms = reshape(terms,1,[]);
        end

        function value = get.numberOfTerms(obj)
            value = numel(obj.terms);
        end
    end

    methods (Access = protected)
        function orders = normalizeApproximationOrders( ...
                obj,approximationOrders)
            arguments
                obj
                approximationOrders (1,:) double ...
                    {mustBeInteger,mustBePositive}
            end

            if isscalar(approximationOrders)
                orders = repmat( ...
                    approximationOrders,1,obj.numberOfTerms);
            elseif numel(approximationOrders) == obj.numberOfTerms
                orders = reshape(approximationOrders,1,[]);
            else
                error( ...
                    'MultipoleSrcSeries:ApproximationOrderSizeMismatch', ...
                    ['approximationOrders must be scalar or contain one ' ...
                     'entry per multipole source term.']);
            end
        end
    end

    methods (Abstract)
        source = discretize(obj,varargin)
    end
end
