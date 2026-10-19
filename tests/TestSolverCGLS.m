classdef TestSolverCGLS < matlab.unittest.TestCase
    % Unit tests for solver_CGLS.

    methods (Test)
        function trueErrorTrackingIsOptional(testCase)
            M = [1,2;3,-1;2,1];
            A = MatrixLinearOp(M);
            xTrue = [2;-1];
            b = M*xTrue;

            [x,info] = solver_CGLS(A,b);
            [xEmpty,infoEmpty] = solver_CGLS(A,b,TrueSolution=[]);
            [xTracked,infoTracked] = solver_CGLS(A,b,TrueSolution=xTrue);

            testCase.verifyEmpty(info.trueErrorNorm);
            testCase.verifyEmpty(info.trueErrorHistory);
            testCase.verifyEqual(xEmpty,x);
            testCase.verifyEqual(infoEmpty,info);
            testCase.verifyEqual(xTracked,x);
            testCase.verifyEqual(infoTracked.flag,info.flag);
            testCase.verifyEqual(infoTracked.iterations,info.iterations);
            testCase.verifyEqual( ...
                infoTracked.dataResidualHistory,info.dataResidualHistory);
            testCase.verifyEqual( ...
                infoTracked.normalResidualHistory,info.normalResidualHistory);
        end

        function trueErrorHistoryMatchesEachIterate(testCase)
            M = diag([1,2,4]);
            A = MatrixLinearOp(M);
            xTrue = [2;-1;3];
            x0 = [-1;2;0];
            b = M*xTrue;

            [x,info] = solver_CGLS(A,b,InitialGuess=x0, ...
                TrueSolution=xTrue,RelativeTolerance=1e-12,MaxIterations=3);

            testCase.verifyNumElements(info.trueErrorHistory,info.iterations+1);
            testCase.verifyEqual(info.trueErrorHistory(1),norm(x0-xTrue));
            testCase.verifyEqual(info.trueErrorNorm,norm(x-xTrue),'AbsTol',1e-14);
            testCase.verifyEqual(info.trueErrorHistory(end),info.trueErrorNorm);
            for k = 1:info.iterations
                xk = solver_CGLS(A,b,InitialGuess=x0, ...
                    RelativeTolerance=1e-12,MaxIterations=k);
                testCase.verifyEqual( ...
                    info.trueErrorHistory(k+1),norm(xk-xTrue),'AbsTol',1e-14);
            end
        end

        function recordsComplexTrueError(testCase)
            M = [1+2i,2;3i,-1i;2-1i,4];
            A = MatrixLinearOp(M);
            xTrue = [2-i;-1+3i];

            [x,info] = solver_CGLS(A,M*xTrue, ...
                TrueSolution=xTrue,RelativeTolerance=1e-11);

            testCase.verifyTrue(isreal(info.trueErrorHistory));
            testCase.verifyEqual(info.trueErrorHistory(1),norm(xTrue));
            testCase.verifyEqual(info.trueErrorNorm,norm(x-xTrue),'AbsTol',1e-14);
        end

        function stationaryGuessCanHaveNonzeroTrueError(testCase)
            % Different solutions can differ by a null-space component.
            A = MatrixLinearOp([1,0;0,0]);
            xTrue = [2;3];
            x0 = [2;0];

            [x,info] = solver_CGLS(A,A.fwd(xTrue), ...
                InitialGuess=x0,TrueSolution=xTrue);

            testCase.verifyEqual(x,x0);
            testCase.verifyEqual(info.iterations,0);
            testCase.verifyTrue(info.converged);
            testCase.verifyEqual(info.trueErrorHistory,3);
            testCase.verifyEqual(info.trueErrorNorm,3);
        end

        function trueErrorHistoryTrimsAtIterationLimit(testCase)
            A = MatrixLinearOp(diag([1,2,4]));
            xTrue = ones(3,1);

            [x,info] = solver_CGLS(A,A.fwd(xTrue), ...
                TrueSolution=xTrue,MaxIterations=1,RelativeTolerance=1e-14);

            testCase.verifyEqual(info.flag,1);
            testCase.verifyNumElements(info.trueErrorHistory,2);
            testCase.verifyEqual(info.trueErrorHistory(end),norm(x-xTrue));
        end

        function rejectsIncorrectTrueSolutionShape(testCase)
            A = MatrixLinearOp(ones(3,2));

            testCase.verifyError( ...
                @() solver_CGLS(A,ones(3,1),TrueSolution=zeros(3,1)), ...
                'solver_CGLS:InvalidTrueSolution');
            testCase.verifyError( ...
                @() solver_CGLS(A,ones(3,1),TrueSolution=[1,2]), ...
                'solver_CGLS:InvalidTrueSolution');
        end

        function solvesConsistentProblemWithDefaults(testCase)
            M = [1,2;3,-1;2,1;1,-3];
            xExact = [2;-1];
            A = MatrixLinearOp(M);

            [x,info] = solver_CGLS(A,M*xExact);

            testCase.verifyEqual(x,xExact,'AbsTol',1e-11);
            testCase.verifyTrue(info.converged);
            testCase.verifyEqual(info.flag,0);
            testCase.verifyGreaterThan(info.iterations,0);
            testCase.verifyNumElements( ...
                info.normalResidualHistory,info.iterations + 1);
        end

        function solvesInconsistentProblem(testCase)
            M = [1,0;1,1;1,2;1,3];
            b = [1;2;2;5];
            expected = M\b;
            A = MatrixLinearOp(M);

            [x,info] = solver_CGLS(A,b,RelativeTolerance=1e-12);

            testCase.verifyEqual(x,expected,'AbsTol',1e-11);
            testCase.verifyTrue(info.converged);
        end

        function supportsNonzeroInitialGuess(testCase)
            M = [2,-1;1,3;4,1];
            xExact = [-2;3];
            A = MatrixLinearOp(M);

            x = solver_CGLS(A,M*xExact, ...
                InitialGuess=[10;-4],RelativeTolerance=1e-12);

            testCase.verifyEqual(x,xExact,'AbsTol',1e-11);
        end

        function supportsComplexProblems(testCase)
            M = [1+2i,2;3i,-1i;2-1i,4];
            xExact = [2-i;-1+3i];
            A = MatrixLinearOp(M);

            [x,info] = solver_CGLS(A,M*xExact,RelativeTolerance=1e-11);

            testCase.verifyEqual(x,xExact,'AbsTol',1e-10);
            testCase.verifyTrue(info.converged);
        end

        function returnsImmediatelyForStationaryInitialGuess(testCase)
            A = MatrixLinearOp([1,2;3,4;5,6]);
            x0 = [2;-1];
            b = A.fwd(x0);

            [x,info] = solver_CGLS(A,b,InitialGuess=x0);

            testCase.verifyEqual(x,x0);
            testCase.verifyEqual(info.iterations,0);
            testCase.verifyTrue(info.converged);
            testCase.verifyNumElements(info.dataResidualHistory,1);
            testCase.verifyNumElements(info.normalResidualHistory,1);
            testCase.verifyEqual(info.dataResidualNorm,0);
            testCase.verifyEqual(info.normalResidualNorm,0);
        end

        function historiesIncludeInitialAndFinalResiduals(testCase)
            M = [1,2;3,-1;2,1];
            b = [1;2;4];
            A = MatrixLinearOp(M);

            [x,info] = solver_CGLS(A,b,RelativeTolerance=1e-12);

            testCase.verifyEqual(info.dataResidualHistory(1),norm(b));
            testCase.verifyEqual( ...
                info.normalResidualHistory(1),norm(A.adj(b)));
            testCase.verifyEqual( ...
                info.dataResidualHistory(end),norm(b-A.fwd(x)), ...
                'AbsTol',1e-14);
            testCase.verifyEqual( ...
                info.normalResidualHistory(end), ...
                norm(A.adj(b-A.fwd(x))),'AbsTol',1e-14);
        end

        function reportsMaximumIterationTermination(testCase)
            M = [1,2,0;0,1,3;2,-1,1;1,0,1];
            A = MatrixLinearOp(M);
            b = [1;2;3;4];

            [~,info] = solver_CGLS(A,b, ...
                RelativeTolerance=1e-14,MaxIterations=1);

            testCase.verifyEqual(info.iterations,1);
            testCase.verifyFalse(info.converged);
            testCase.verifyEqual(info.flag,1);
            testCase.verifyNumElements(info.dataResidualHistory,2);
            testCase.verifyNumElements(info.normalResidualHistory,2);
        end

        function rejectsIncorrectDataLength(testCase)
            A = MatrixLinearOp(ones(3,2));

            testCase.verifyError(@() solver_CGLS(A,[1;2]), ...
                'solver_CGLS:InvalidDataVector');
        end

        function rejectsIncorrectInitialGuessLength(testCase)
            A = MatrixLinearOp(ones(3,2));

            testCase.verifyError( ...
                @() solver_CGLS(A,ones(3,1),InitialGuess=zeros(3,1)), ...
                'solver_CGLS:InvalidInitialGuess');
        end

        function rejectsRowInitialGuess(testCase)
            A = MatrixLinearOp(ones(3,2));

            testCase.verifyError( ...
                @() solver_CGLS(A,ones(3,1),InitialGuess=[1,2]), ...
                'solver_CGLS:InvalidInitialGuess');
        end
    end
end
