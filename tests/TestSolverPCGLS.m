classdef TestSolverPCGLS < matlab.unittest.TestCase
    methods (Test)
        function solvesConsistentProblem(testCase)
            B = [1,2;3,-1;2,1;1,-3];
            xTrue = [2;-1];
            A = MatrixLinearOp(B);
            M = MatrixInvLinearOp(diag(sum(abs(B).^2,1)));
            [x,info] = solver_PCGLS(A,B*xTrue,M, ...
                TrueSolution=xTrue,RelativeTolerance=1e-11);
            testCase.verifyEqual(x,xTrue,'AbsTol',1e-10);
            testCase.verifyTrue(info.converged);
            testCase.verifyNumElements(info.trueErrorHistory,info.iterations+1);
            testCase.verifyEqual(info.trueErrorHistory(1),norm(xTrue));
            testCase.verifyEqual(info.trueErrorNorm,norm(x-xTrue));
        end

        function solvesInconsistentProblem(testCase)
            B = [1,0;1,1;1,2;1,3];
            b = [1;2;2;5];
            A = MatrixLinearOp(B);
            M = MatrixInvLinearOp([2,0.5;0.5,3]);
            [x,info] = solver_PCGLS(A,b,M,RelativeTolerance=1e-11);
            testCase.verifyEqual(x,B\b,'AbsTol',1e-10);
            testCase.verifyTrue(info.converged);
            testCase.verifyGreaterThan(info.dataResidualNorm,0.1);
            testCase.verifyEmpty(info.trueErrorHistory);
            testCase.verifyEmpty(info.trueErrorNorm);
        end

        function identityMatchesUnpreconditionedCGLS(testCase)
            B = [1,2,0;0,1,3;2,-1,1;1,0,1];
            A = MatrixLinearOp(B);
            b = [1;2;3;4];
            M = MatrixInvLinearOp(eye(3));
            [x,info] = solver_PCGLS(A,b,M, ...
                MaxIterations=2,RelativeTolerance=1e-25);
            [expected,reference] = solver_CGLS(A,b, ...
                MaxIterations=2,RelativeTolerance=1e-25);
            testCase.verifyEqual(x,expected,'AbsTol',1e-13);
            testCase.verifyEqual(info.dataResidualHistory, ...
                reference.dataResidualHistory,'AbsTol',1e-13);
            testCase.verifyEqual(info.normalResidualHistory, ...
                reference.normalResidualHistory,'AbsTol',1e-13);
        end

        function agreesWithSplitPreconditioningForComplexData(testCase)
            % Independent reference: CGLS applied to B/L, y=L*x.
            B = [1+2i,2,0;3i,-1i,1;2-1i,4,2i;1,0,3];
            C = [2,1i,0;0,3,1;0,0,2];
            H = C'*C;
            L = chol(H);
            A = MatrixLinearOp(B);
            M = MatrixInvLinearOp(H);
            b = [1+i;2;-1i;3];
            x0 = [1;-2i;0.5];
            [x,info] = solver_PCGLS(A,b,M,InitialGuess=x0, ...
                RelativeTolerance=1e-25,MaxIterations=2);
            y = solver_CGLS(MatrixLinearOp(B/L),b, ...
                InitialGuess=L*x0,RelativeTolerance=1e-25,MaxIterations=2);
            testCase.verifyEqual(x,L\y,'AbsTol',1e-11);
            testCase.verifyEqual(info.flag,1);
            testCase.verifyTrue(isreal(info.preconditionedNormalResidualHistory));
        end

        function exactPreconditionerSolvesInOneStep(testCase)
            B = [1,2;3,-1;2,1;1,-3];
            A = MatrixLinearOp(B);
            M = MatrixInvLinearOp(B'*B);
            xTrue = [2;-1];
            x0 = [4;3];
            [x,info] = solver_PCGLS(A,B*xTrue,M, ...
                InitialGuess=x0,TrueSolution=xTrue,RelativeTolerance=1e-10);
            testCase.verifyEqual(info.iterations,1);
            testCase.verifyTrue(info.converged);
            testCase.verifyEqual(x,xTrue,'AbsTol',1e-12);
            testCase.verifyEqual(info.trueErrorHistory(1),norm(x0-xTrue));
        end

        function historiesMatchEveryCompletedIterate(testCase)
            B = [1,2,0;0,1,3;2,-1,1;1,0,1];
            A = MatrixLinearOp(B);
            H = [2,0.3,0;0.3,3,0.2;0,0.2,1];
            M = MatrixInvLinearOp(H);
            xTrue = [2;-1;3];
            b = B*xTrue;
            x0 = [1;1;1];
            [~,info] = solver_PCGLS(A,b,M,InitialGuess=x0, ...
                TrueSolution=xTrue,RelativeTolerance=1e-25,MaxIterations=2);
            for k = 0:info.iterations
                if k == 0
                    xk = x0;
                else
                    xk = solver_PCGLS(A,b,M,InitialGuess=x0, ...
                        RelativeTolerance=1e-25,MaxIterations=k);
                end
                r = b-B*xk;
                g = B'*r;
                testCase.verifyEqual(info.dataResidualHistory(k+1), ...
                    norm(r),'AbsTol',1e-12);
                testCase.verifyEqual(info.normalResidualHistory(k+1), ...
                    norm(g),'AbsTol',1e-12);
                testCase.verifyEqual(info.preconditionedNormalResidualHistory(k+1), ...
                    sqrt(real(g'*(H\g))),'AbsTol',1e-12);
                testCase.verifyEqual(info.trueErrorHistory(k+1), ...
                    norm(xk-xTrue),'AbsTol',1e-12);
            end
            testCase.verifyEqual(info.trueErrorNorm,info.trueErrorHistory(end));
        end

        function stationaryInitialGuessHasOneHistoryEntry(testCase)
            A = MatrixLinearOp([1,0;1,1;1,2]);
            M = MatrixInvLinearOp(eye(2));
            x0 = [2;-1];
            [x,info] = solver_PCGLS(A,A.fwd(x0),M, ...
                InitialGuess=x0,TrueSolution=x0);
            testCase.verifyEqual(x,x0);
            testCase.verifyEqual(info.iterations,0);
            testCase.verifyTrue(info.converged);
            testCase.verifyEqual(info.dataResidualHistory,0);
            testCase.verifyEqual(info.normalResidualHistory,0);
            testCase.verifyEqual(info.preconditionedNormalResidualHistory,0);
            testCase.verifyEqual(info.trueErrorHistory,0);
        end

        function smallDataResidualDoesNotTriggerConvergence(testCase)
            % Data residual is small, but the normal residual is not.
            [x,info] = solver_PCGLS(MatrixLinearOp(1e4),1e-4, ...
                MatrixInvLinearOp(1), ...
                RelativeTolerance=0,AbsoluteTolerance=1e-3);
            testCase.verifyEqual(info.iterations,1);
            testCase.verifyTrue(info.converged);
            testCase.verifyEqual(x,1e-8,'AbsTol',1e-20);
        end

        function preconditionerScalingDoesNotTriggerConvergence(testCase)
            % sqrt(g'*M.inv(g)) starts below the absolute tolerance.
            A = MatrixLinearOp(eye(2));
            b = [1;2];
            [x,info] = solver_PCGLS(A,b,MatrixInvLinearOp(1e12*eye(2)), ...
                RelativeTolerance=0,AbsoluteTolerance=1e-3);
            [xRef,reference] = solver_PCGLS(A,b,MatrixInvLinearOp(eye(2)), ...
                RelativeTolerance=0,AbsoluteTolerance=1e-3);
            testCase.verifyEqual(info.iterations,1);
            testCase.verifyTrue(info.converged);
            testCase.verifyEqual(x,b,'AbsTol',1e-12);
            testCase.verifyEqual(x,xRef,'AbsTol',1e-12);
            testCase.verifyEqual(info.stoppingThreshold,reference.stoppingThreshold);
            testCase.verifyEqual(info.iterations,reference.iterations);
            testCase.verifyLessThan(info.preconditionedNormalResidualHistory(1),1e-3);
        end

        function trueSolutionDoesNotChangeIterates(testCase)
            A = MatrixLinearOp(diag([1,2,4]));
            M = MatrixInvLinearOp(diag([1,2,3]));
            b = [1;2;3];
            [x,info] = solver_PCGLS(A,b,M,MaxIterations=1,RelativeTolerance=1e-20);
            [xt,tracked] = solver_PCGLS(A,b,M,TrueSolution=[7;8;9], ...
                MaxIterations=1,RelativeTolerance=1e-20);
            testCase.verifyEqual(xt,x);
            testCase.verifyEqual(tracked.flag,info.flag);
            testCase.verifyEqual(tracked.dataResidualHistory,info.dataResidualHistory);
            testCase.verifyEqual(info.flag,1);
            testCase.verifyNumElements(tracked.trueErrorHistory,2);
        end

        function reportsInvalidPreconditionedInnerProducts(testCase)
            A = MatrixLinearOp(eye(2));
            for H = {-eye(2),1i*eye(2)}
                [~,info] = solver_PCGLS(A,[1;2],MatrixInvLinearOp(H{1}));
                testCase.verifyEqual(info.flag,2);
                testCase.verifyFalse(info.converged);
                testCase.verifyEqual(info.iterations,0);
                testCase.verifyTrue(isnan(info.preconditionedNormalResidualNorm));
            end
        end

        function rejectsIncorrectDimensions(testCase)
            A = MatrixLinearOp(ones(3,2));
            M = MatrixInvLinearOp(eye(2));
            testCase.verifyError(@() solver_PCGLS(A,ones(3,1), ...
                MatrixInvLinearOp(eye(3))), ...
                'solver_PCGLS:InvalidPreconditionerDimension');
            testCase.verifyError(@() solver_PCGLS(A,ones(2,1),M), ...
                'solver_PCGLS:InvalidDataVector');
            testCase.verifyError(@() solver_PCGLS(A,ones(3,1),M, ...
                InitialGuess=[1,2]),'solver_PCGLS:InvalidInitialGuess');
            testCase.verifyError(@() solver_PCGLS(A,ones(3,1),M, ...
                TrueSolution=zeros(3,1)),'solver_PCGLS:InvalidTrueSolution');
        end
    end
end
