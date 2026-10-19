classdef TestInverseActions < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addHelperPath(tc)
            helperPath = fullfile(fileparts(mfilename('fullpath')),'helpers');
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(helperPath));
        end
    end

    methods (Test)
        function inverseActionsSolveComplexSystems(tc)
            A = [3+1i,1-2i;2,5-1i];
            op = MatrixInvLinearOp(A);
            b = [2+3i;-1+2i];
            x = op.inv(b);
            z = op.invAdj(b);
            tc.verifyEqual(A*x,b,'AbsTol',1e-12);
            tc.verifyEqual(A'*z,b,'AbsTol',1e-12);
            tc.verifyEqual(x,A\b,'AbsTol',1e-12);
            tc.verifyEqual(z,A'\b,'AbsTol',1e-12);
            y = [1-1i;3+2i];
            tc.verifyEqual(dot(op.inv(b),y),dot(b,op.invAdj(y)), ...
                'AbsTol',1e-12);
        end

        function sparseInverseActions(tc)
            A = sparse([4,1i,0;2,5,1;0,1,3-1i]);
            op = MatrixInvLinearOp(A);
            b = [1;2+1i;3];
            tc.verifyEqual(A*op.inv(b),b,'AbsTol',1e-12);
            tc.verifyEqual(A'*op.invAdj(b),b,'AbsTol',1e-12);
        end

        function inverseWrapperDelegatesAllFourActions(tc)
            A = [3+1i,1-2i;2,5-1i];
            op = MatrixInvLinearOp(A);
            B = op.inverse();
            b = [2+3i;-1+2i];
            tc.verifyClass(B,'InverseOfLinearOp');
            tc.verifyTrue(isa(B,'InvLinearOp'));
            tc.verifyEqual(B.fwd(b),op.inv(b),'AbsTol',1e-12);
            tc.verifyEqual(B.adj(b),op.invAdj(b),'AbsTol',1e-12);
            tc.verifyEqual(B.inv(b),op.fwd(b),'AbsTol',1e-12);
            tc.verifyEqual(B.invAdj(b),op.adj(b),'AbsTol',1e-12);
            tc.verifyEqual(B.inverse(),op);
            M = B.matrixRep();
            tc.verifyEqual(M*A,eye(2),'AbsTol',1e-12);
            tc.verifyEqual(A*M,eye(2),'AbsTol',1e-12);
        end

        function adjointWrapperDelegatesInverseActions(tc)
            op = MatrixInvLinearOp([3+1i,1-2i;2,5-1i]);
            H = op.adjoint();
            b = [2+3i;-1+2i];
            tc.verifyEqual(H.inv(b),op.invAdj(b),'AbsTol',1e-12);
            tc.verifyEqual(H.invAdj(b),op.inv(b),'AbsTol',1e-12);
            B = H.inverse();
            D = op.inverse();
            K = D.adjoint();
            tc.verifyEqual(B.fwd(b),K.fwd(b),'AbsTol',1e-12);
            tc.verifyEqual(B.adj(b),K.adj(b),'AbsTol',1e-12);
            tc.verifyEqual(B.matrixRep(),K.matrixRep(),'AbsTol',1e-12);
        end

        function wrappersAndActionsDoNotRequireMatrixAssembly(tc)
            d = [2+1i;3-2i];
            op = MatrixFreeInvOpForTest(d);
            B = op.inverse();
            H = op.adjoint();
            BH = B.adjoint();
            HB = H.inverse();
            b = [1+2i;3-1i];
            tc.verifyEqual(B.fwd(b),b./d);
            tc.verifyEqual(B.adj(b),b./conj(d));
            tc.verifyEqual(B.inv(b),d.*b);
            tc.verifyEqual(B.invAdj(b),conj(d).*b);
            tc.verifyEqual(H.inv(b),b./conj(d));
            tc.verifyEqual(H.invAdj(b),b./d);
            tc.verifyEqual(BH.fwd(b),HB.fwd(b));
            tc.verifyEqual(B.inverse(),op);
            tc.verifyEqual(H.adjoint(),op);
            tc.verifyError(@() B.matrixRep(), ...
                'MatrixFreeInvOpForTest:MatrixForbidden');
        end

        function inverseActionsRejectInvalidVectors(tc)
            op = MatrixInvLinearOp(eye(2));
            B = op.inverse();
            H = op.adjoint();
            ops = {op,B,H};
            for k = 1:numel(ops)
                current = ops{k};
                tc.verifyError(@() current.inv([1,2]), ...
                    'LinearOp:InvalidRangeVector');
                tc.verifyError(@() current.inv(ones(3,1)), ...
                    'LinearOp:InvalidRangeVector');
                tc.verifyError(@() current.invAdj([1,2]), ...
                    'LinearOp:InvalidDomainVector');
                tc.verifyError(@() current.invAdj(ones(3,1)), ...
                    'LinearOp:InvalidDomainVector');
            end
        end
    end
end
