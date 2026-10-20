classdef MatrixLinearOp < LinearOp
    % MatrixLinearOp represents a finite-dimensional linear operator
    % through its matrix representation.
    %
    % If A has size m-by-n, then the operator maps an n-dimensional
    % domain into an m-dimensional range:
    %
    %   y = A*x
    %
    % The adjoint operator is represented by the conjugate transpose A'.

    properties (SetAccess = private)
        matrix (:,:) double
    end

    methods
        function obj = MatrixLinearOp(A)
            arguments
                A (:,:) double {mustBeNonempty}
            end

            [numRows,numCols] = size(A);
            obj@LinearOp(numCols,numRows);
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
            % Return the matrix representation in the standard bases.
            A = obj.matrix;
        end
    end
end