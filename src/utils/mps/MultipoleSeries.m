classdef (Abstract) MultipoleSeries
    % MultipoleSeries is the common base for finite multipole collections.

    properties (SetAccess = private)
        terms
        approximationOrders (1,:) double
    end

    properties (Dependent)
        numberOfTerms
    end

    properties (Abstract, Dependent)
        locations
    end

    methods
        function obj = MultipoleSeries(terms,approximationOrders)
            arguments
                terms
                approximationOrders (1,:) double ...
                    {mustBeInteger,mustBePositive}
            end

            if ~isempty(terms) && ~isa(terms,'MultipoleTerm')
                error('MultipoleSeries:InvalidTerms', ...
                    'Every term must inherit from MultipoleTerm.');
            end

            terms = reshape(terms,1,[]);
            nTerms = numel(terms);
            if isscalar(approximationOrders)
                approximationOrders = repmat( ...
                    approximationOrders,1,nTerms);
            elseif numel(approximationOrders) ~= nTerms
                error('MultipoleSeries:ApproximationOrderSizeMismatch', ...
                    ['ApproximationOrders must be scalar or contain one ' ...
                     'entry per multipole term.']);
            end

            obj.terms = terms;
            obj.approximationOrders = reshape(approximationOrders,1,[]);
        end

        function value = get.numberOfTerms(obj)
            value = numel(obj.terms);
        end
    end

    methods (Abstract)
        source = discretize(obj,varargin)
    end
end
