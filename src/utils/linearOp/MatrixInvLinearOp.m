classdef MatrixInvLinearOp < InvLinearOp
    % Invertible matrix operator with direct inverse actions.
    % inv(b) uses A\b; invAdj(b) uses A'\b. inverse() is inherited and
    % returns an InverseOfLinearOp wrapper without assembling A^(-1).

    properties (SetAccess = private)
        matrix (:,:) double
    end

    methods
        function obj = MatrixInvLinearOp(A)
            arguments
                A (:,:) double {mustBeNonempty,mustBeFinite}
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
            obj.validateDomainVector(x);
            y = obj.matrix*x;
        end

        function x = adj(obj,y)
            obj.validateRangeVector(y);
            x = obj.matrix'*y;
        end

        function x = inv(obj,b)
            obj.validateRangeVector(b);
            x = obj.matrix\b;
        end

        function x = invAdj(obj,b)
            obj.validateDomainVector(b);
            x = obj.matrix'\b;
        end

        function A = matrixRep(obj)
            A = obj.matrix;
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
