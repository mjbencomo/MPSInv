classdef TestLinearOpOperatorMethods < matlab.unittest.TestCase
    methods (Test)
        function rectangularAdjointExchangesDimensions(tc)
            A = [1,2,3;4,5,6];
            op = MatrixLinearOp(A);
            H = op.adjoint();
            tc.verifyTrue(isa(H,'LinearOp'));
            tc.verifyEqual(H.dimDomain,2);
            tc.verifyEqual(H.dimRange,3);
            tc.verifyEqual(H.fwd([2;3]),A'*[2;3]);
            tc.verifyEqual(H.adj([1;2;3]),A*[1;2;3]);
            tc.verifyEqual(H.matrixRep(),A');
        end

        function complexAdjointUsesConjugation(tc)
            A = [1+2i,3,-1i;2,1i,4];
            op = MatrixLinearOp(A);
            H = op.adjoint();
            x = [2+1i;-1+3i];
            y = [1-2i;3+1i;2];
            tc.verifyEqual(H.fwd(x),A'*x,'AbsTol',1e-12);
            tc.verifyEqual(H.adj(y),A*y,'AbsTol',1e-12);
            tc.verifyEqual(H.matrixRep(),A');
            tc.verifyEqual(dot(H.fwd(x),y),dot(x,H.adj(y)), ...
                'AbsTol',1e-12);
        end

        function doubleAdjointReturnsOriginalOperator(tc)
            op = MatrixLinearOp([1,2,3;4,5,6]);
            H = op.adjoint();
            recovered = H.adjoint();
            tc.verifyClass(recovered,'MatrixLinearOp');
            tc.verifyEqual(recovered,op);
        end

        function invertibleAdjointPreservesInterface(tc)
            A = [2+1i,1;3,4-2i];
            op = MatrixInvLinearOp(A);
            H = op.adjoint();
            tc.verifyTrue(isa(H,'InvLinearOp'));
            tc.verifyEqual(H.matrixRep(),A');
            recovered = H.adjoint();
            tc.verifyEqual(recovered,op);
            tc.verifyClass(recovered,'MatrixInvLinearOp');
        end

        function inverseRemainsAnOperator(tc)
            A = [3,1;2,4];
            op = MatrixInvLinearOp(A);
            B = op.inverse();
            tc.verifyTrue(isa(B,'InvLinearOp'));
            x = [2;-3];
            tc.verifyEqual(B.fwd(op.fwd(x)),x,'AbsTol',1e-12);
            tc.verifyEqual(op.fwd(B.fwd(x)),x,'AbsTol',1e-12);
            z = [1;5];
            tc.verifyEqual(B.adj(z),A'\z,'AbsTol',1e-12);
        end

        function inverseAndAdjointCommute(tc)
            A = [2+1i,1;3,4-2i];
            op = MatrixInvLinearOp(A);
            H = op.adjoint();
            B = H.inverse();
            baseInverse = op.inverse();
            expected = baseInverse.adjoint();
            tc.verifyTrue(isa(B,'InvLinearOp'));
            tc.verifyEqual(B.matrixRep(),expected.matrixRep(),'AbsTol',1e-12);
            x = [2+1i;-3+2i];
            tc.verifyEqual(B.fwd(H.fwd(x)),x,'AbsTol',1e-12);
            tc.verifyEqual(H.fwd(B.fwd(x)),x,'AbsTol',1e-12);
            tc.verifyEqual(B.fwd(x),A'\x,'AbsTol',1e-12);
            tc.verifyEqual(B.adj(x),A\x,'AbsTol',1e-12);
            recovered = B.inverse();
            tc.verifyEqual(recovered.matrixRep(),H.matrixRep(),'AbsTol',1e-12);
        end

        function sparseRepresentationRemainsSparse(tc)
            op = MatrixLinearOp(sparse([1,0,2;0,3,0]));
            H = op.adjoint();
            tc.verifyTrue(issparse(H.matrixRep()));
            tc.verifyEqual(H.matrixRep(),op.matrixRep()');
        end

        function wrappersValidateVectors(tc)
            op = MatrixLinearOp(ones(2,3));
            H = op.adjoint();
            tc.verifyError(@() H.fwd([1,2]), ...
                'LinearOp:InvalidDomainVector');
            tc.verifyError(@() H.fwd(ones(3,1)), ...
                'LinearOp:InvalidDomainVector');
            tc.verifyError(@() H.adj([1,2,3]), ...
                'LinearOp:InvalidRangeVector');
            tc.verifyError(@() H.adj(ones(2,1)), ...
                'LinearOp:InvalidRangeVector');
            invertible = MatrixInvLinearOp(eye(2));
            HI = invertible.adjoint();
            tc.verifyError(@() HI.fwd([1,2]), ...
                'LinearOp:InvalidDomainVector');
            tc.verifyError(@() HI.adj([1,2]), ...
                'LinearOp:InvalidRangeVector');
        end

    end
end
