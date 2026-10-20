classdef TestMatrixInvLinearOp < matlab.unittest.TestCase
    % Unit tests for MatrixInvLinearOp.

    methods (Test)
        function constructorStoresMatrixAndSquareDimensions(testCase)
            A = [2,1;1,3];

            op = MatrixInvLinearOp(A);

            testCase.verifyEqual(op.dimDomain,2);
            testCase.verifyEqual(op.dimRange,2);
            testCase.verifyEqual(op.matrix,A);
        end

        function forwardAndAdjointMatchMatrixActions(testCase)
            A = [1+2i,2;-3i,4-i];
            x = [2-i;3+2i];
            y = [-1+3i;2-i];
            op = MatrixInvLinearOp(A);

            testCase.verifyEqual(op.fwd(x),A*x);
            testCase.verifyEqual(op.adj(y),A'*y);
        end

        function inverseReturnsMatrixInvLinearOp(testCase)
            op = MatrixInvLinearOp([3,1;1,2]);

            inverseOp = op.inverse();

            testCase.verifyClass(inverseOp,'MatrixInvLinearOp');
            testCase.verifyEqual(inverseOp.dimDomain,op.dimRange);
            testCase.verifyEqual(inverseOp.dimRange,op.dimDomain);
        end

        function inverseUndoesForwardAction(testCase)
            A = [4,1,-1;2,5,1;1,-2,6];
            x = [2;-3;1];
            op = MatrixInvLinearOp(A);
            inverseOp = op.inverse();

            recovered = inverseOp.fwd(op.fwd(x));

            testCase.verifyEqual(recovered,x,'AbsTol',1e-12);
        end

        function inverseMatrixRepSatisfiesTwoSidedIdentity(testCase)
            A = [2+1i,1;3,4-2i];
            op = MatrixInvLinearOp(A);
            inverseOp = op.inverse();
            inverseMatrix = inverseOp.matrixRep();
            identity = eye(2);

            testCase.verifyEqual(inverseMatrix*A,identity,'AbsTol',1e-12);
            testCase.verifyEqual(A*inverseMatrix,identity,'AbsTol',1e-12);
        end

        function inverseOfInverseRecoversOriginalOperator(testCase)
            A = [2,-1;3,4];
            op = MatrixInvLinearOp(A);

            inverseOp = op.inverse();
            recoveredOp = inverseOp.inverse();
            recoveredMatrix = recoveredOp.matrixRep();

            testCase.verifyEqual(recoveredMatrix,A,'AbsTol',1e-12);
        end

        function supportsSparseMatrices(testCase)
            A = sparse([3,0,1;0,2,0;1,0,4]);
            x = [1;2;3];
            op = MatrixInvLinearOp(A);
            inverseOp = op.inverse();

            testCase.verifyTrue(issparse(op.matrixRep()));
            testCase.verifyEqual(inverseOp.fwd(op.fwd(x)),x, ...
                'AbsTol',1e-12);
        end

        function rejectsNonsquareMatrix(testCase)
            testCase.verifyError(@() MatrixInvLinearOp(ones(2,3)), ...
                'MatrixInvLinearOp:NonSquareMatrix');
        end

        function rejectsSingularMatrix(testCase)
            testCase.verifyError(@() MatrixInvLinearOp([1,2;2,4]), ...
                'MatrixInvLinearOp:SingularMatrix');
        end

        function forwardRejectsInvalidVector(testCase)
            op = MatrixInvLinearOp(eye(2));

            testCase.verifyError(@() op.fwd([1,2]), ...
                'LinearOp:InvalidDomainVector');
        end

        function adjointRejectsInvalidVector(testCase)
            op = MatrixInvLinearOp(eye(2));

            testCase.verifyError(@() op.adj([1,2]), ...
                'LinearOp:InvalidRangeVector');
        end
    end
end
