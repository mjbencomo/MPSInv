classdef (Abstract) LinearOp
    % LinearOp is a finite-dimensional linear operator.
    %
    % Subclasses implement the forward action, adjoint action, and matrix
    % representation with respect to the standard coordinate bases.

    properties (SetAccess = protected)
        dimDomain (1,1) double
        dimRange (1,1) double
    end

    methods
        function obj = LinearOp(dimDomain,dimRange)
            arguments
                dimDomain (1,1) double ...
                    {mustBeInteger,mustBePositive}
                dimRange (1,1) double ...
                    {mustBeInteger,mustBePositive}
            end

            obj.dimDomain = dimDomain;
            obj.dimRange = dimRange;
        end
    end

    methods (Access = protected)
        function validateDomainVector(obj,x)
            if ~isnumeric(x) || ~iscolumn(x) || ...
                    numel(x) ~= obj.dimDomain
                error('LinearOp:InvalidDomainVector', ...
                    ['The input must be a numeric column vector with ' ...
                     'dimDomain entries.']);
            end
        end

        function validateRangeVector(obj,y)
            if ~isnumeric(y) || ~iscolumn(y) || ...
                    numel(y) ~= obj.dimRange
                error('LinearOp:InvalidRangeVector', ...
                    ['The input must be a numeric column vector with ' ...
                     'dimRange entries.']);
            end
        end
    end

    methods (Abstract)
        y = fwd(obj,x)
        % Apply the operator to a coordinate vector in its domain.

        x = adj(obj,y)
        % Apply the adjoint operator to a coordinate vector in its range.

        A = matrixRep(obj)
        % Return the matrix representation in the standard bases.
    end
end
