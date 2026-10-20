classdef TestMatrixLinearOp < matlab.unittest.TestCase
    % Unit tests for MatrixLinearOp.

    methods (Test)
        function constructorStoresMatrixAndDimensions(testCase)
            A = [1,2,3;4,5,6];

            op = MatrixLinearOp(A);

            testCase.verifyEqual(op.dimDomain,3);
            testCase.verifyEqual(op.dimRange,2);
            testCase.verifyEqual(op.matrix,A);
        end

        function matrixRepReturnsConstructorMatrix(testCase)
            A = [1,-2,3;4,0,5];
            op = MatrixLinearOp(A);

            testCase.verifyEqual(op.matrixRep(),A);
        end

        function forwardActionMatchesMatrixMultiplication(testCase)
            A = [1,2,-1;0,3,4];
            x = [2;-1;3];
            op = MatrixLinearOp(A);

            testCase.verifyEqual(op.fwd(x),A*x);
        end

        function adjointActionMatchesConjugateTranspose(testCase)
            A = [1+2i,2-1i,-3;4i,5,6+2i];
            y = [2-1i;-3+4i];
            op = MatrixLinearOp(A);

            testCase.verifyEqual(op.adj(y),A'*y);
        end

        function forwardActionSupportsComplexVectors(testCase)
            A = [1,2i;-3i,4];
            x = [2-1i;3+2i];
            op = MatrixLinearOp(A);

            testCase.verifyEqual(op.fwd(x),A*x);
        end

        function forwardAndAdjointSatisfyAdjointIdentity(testCase)
            A = [1+2i,-2,3i;4,5-1i,-6i];
            x = [2-1i;-3+2i;1+4i];
            y = [-1+3i;2-2i];
            op = MatrixLinearOp(A);

            lhs = y'*op.fwd(x);
            rhs = op.adj(y)'*x;

            testCase.verifyEqual(lhs,rhs,'AbsTol',1e-12);
        end

        function rectangularMatrixHasCorrectMappingDimensions(testCase)
            A = reshape(1:20,4,5);
            op = MatrixLinearOp(A);

            x = (1:5).';
            y = (1:4).';

            testCase.verifySize(op.fwd(x),[4,1]);
            testCase.verifySize(op.adj(y),[5,1]);
        end

        function sparseMatrixIsPreserved(testCase)
            A = sparse([1,0,2;0,3,0]);
            x = [1;2;3];
            op = MatrixLinearOp(A);

            testCase.verifyTrue(issparse(op.matrixRep()));
            testCase.verifyEqual(op.fwd(x),A*x);
        end

        function forwardRejectsRowVector(testCase)
            op = MatrixLinearOp(eye(2));

            testCase.verifyError(@() op.fwd([1,2]), ...
                'LinearOp:InvalidDomainVector');
        end

        function forwardRejectsIncorrectVectorLength(testCase)
            op = MatrixLinearOp(ones(2,3));

            testCase.verifyError(@() op.fwd([1;2]), ...
                'LinearOp:InvalidDomainVector');
        end

        function adjointRejectsRowVector(testCase)
            op = MatrixLinearOp(eye(2));

            testCase.verifyError(@() op.adj([1,2]), ...
                'LinearOp:InvalidRangeVector');
        end

        function adjointRejectsIncorrectVectorLength(testCase)
            op = MatrixLinearOp(ones(2,3));

            testCase.verifyError(@() op.adj([1;2;3]), ...
                'LinearOp:InvalidRangeVector');
        end
    end
end
