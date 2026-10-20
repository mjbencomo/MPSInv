classdef MatrixInvLinearOp < InvLinearOp
    % MatrixInvLinearOp represents an invertible matrix linear operator.
    %
    % If A is an invertible n-by-n matrix, then this operator acts as
    %
    %   y = A*x.
    %
    % The inverse method returns another MatrixInvLinearOp representing
    % A^(-1). The inverse matrix is computed with a linear solve rather
    % than by calling inv.

    properties (SetAccess = private)
        matrix (:,:) double
    end

    methods
        function obj = MatrixInvLinearOp(A)
            arguments
                A (:,:) double {mustBeNonempty}
            end

            [numRows,numCols] = size(A);
            if numRows ~= numCols
                error('MatrixInvLinearOp:NonSquareMatrix', ...
                    'The matrix must be square.');
            end

            if MatrixInvLinearOp.isNumericallySingular(A)
                error('MatrixInvLinearOp:SingularMatrix', ...
                    'The matrix must be nonsingular.');
            end

            obj@InvLinearOp(numCols);
            obj.matrix = A;
        end

        function y = fwd(obj,x)
            % Apply the matrix to a vector in the operator domain.
            obj.validateDomainVector(x);
            y = obj.matrix*x;
        end

        function x = adj(obj,y)
            % Apply the Euclidean/Hermitian adjoint.
            obj.validateRangeVector(y);
            x = obj.matrix'*y;
        end

        function A = matrixRep(obj)
            % Return the matrix representation in the standard basis.
            A = obj.matrix;
        end

        function inverseOp = inverse(obj)
            % Return an operator representing the matrix inverse.
            identity = eye(obj.dimDomain,'like',obj.matrix);
            inverseMatrix = obj.matrix\identity;
            inverseOp = MatrixInvLinearOp(inverseMatrix);
        end
    end

    methods (Static, Access = private)
        function tf = isNumericallySingular(A)
            if issparse(A)
                conditionEstimate = condest(A);
                tf = ~isfinite(conditionEstimate);
            else
                reciprocalCondition = rcond(A);
                tf = ~isfinite(reciprocalCondition) || ...
                    reciprocalCondition == 0;
            end
        end
    end
end
