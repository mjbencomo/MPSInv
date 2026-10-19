classdef TestSolverTolerances < matlab.unittest.TestCase
    % Shared stopping-rule regression tests for CGLS and PCGLS.
    methods (Test)
        function defaultsUseRelativeNormalResidual(testCase)
            A = MatrixLinearOp(diag([1,2]));
            b = [1;1];
            for preconditioned = [false,true]
                [~,info] = TestSolverTolerances.solve(A,b,preconditioned, ...
                    MaxIterations=1);
                testCase.verifyEqual(info.stoppingThreshold, ...
                    1e-6*norm(A.adj(b)),'AbsTol',1e-20);
                testCase.verifyEqual(info.flag,1);
            end
        end

        function relativeToleranceControlsTermination(testCase)
            A = MatrixLinearOp(diag([1,2]));
            b = [1;1];
            for preconditioned = [false,true]
                [~,loose] = TestSolverTolerances.solve(A,b,preconditioned, ...
                    RelativeTolerance=0.4,AbsoluteTolerance=0);
                [~,tight] = TestSolverTolerances.solve(A,b,preconditioned, ...
                    RelativeTolerance=0.3,AbsoluteTolerance=0);
                testCase.verifyEqual(loose.iterations,1);
                testCase.verifyEqual(tight.iterations,2);
                testCase.verifyTrue(loose.converged);
                testCase.verifyTrue(tight.converged);
                testCase.verifyLessThanOrEqual( ...
                    loose.normalResidualNorm,loose.stoppingThreshold);
            end
        end

        function absoluteToleranceWorksWithoutRelativeTolerance(testCase)
            A = MatrixLinearOp(diag([1,2]));
            b = [1;1];
            for preconditioned = [false,true]
                [~,loose] = TestSolverTolerances.solve(A,b,preconditioned, ...
                    RelativeTolerance=0,AbsoluteTolerance=0.8);
                [~,tight] = TestSolverTolerances.solve(A,b,preconditioned, ...
                    RelativeTolerance=0,AbsoluteTolerance=0.7);
                testCase.verifyEqual(loose.iterations,1);
                testCase.verifyEqual(tight.iterations,2);
                testCase.verifyEqual(loose.stoppingThreshold,0.8);
                testCase.verifyTrue(loose.converged);
            end
        end

        function thresholdUsesMaximumOfTheTwoTerms(testCase)
            A = MatrixLinearOp(diag([1,2]));
            b = [1;1];
            for preconditioned = [false,true]
                [~,absoluteDominates] = TestSolverTolerances.solve(A,b, ...
                    preconditioned,RelativeTolerance=0.1,AbsoluteTolerance=0.8);
                [~,relativeDominates] = TestSolverTolerances.solve(A,b, ...
                    preconditioned,RelativeTolerance=0.4,AbsoluteTolerance=0.5);
                testCase.verifyEqual(absoluteDominates.stoppingThreshold,0.8);
                testCase.verifyEqual(relativeDominates.stoppingThreshold, ...
                    0.4*sqrt(5),'AbsTol',1e-14);
            end
        end

        function absoluteToleranceCanAcceptNonstationaryInitialGuess(testCase)
            A = MatrixLinearOp(eye(2));
            b = [1e-10;-2e-10];
            xTrue = b;
            for preconditioned = [false,true]
                [x,info] = TestSolverTolerances.solve(A,b,preconditioned, ...
                    RelativeTolerance=1e-6,AbsoluteTolerance=1e-8, ...
                    TrueSolution=xTrue);
                testCase.verifyEqual(x,zeros(2,1));
                testCase.verifyEqual(info.iterations,0);
                testCase.verifyTrue(info.converged);
                testCase.verifyGreaterThan(info.normalResidualNorm,0);
                testCase.verifyNumElements(info.normalResidualHistory,1);
                testCase.verifyEqual(info.trueErrorHistory,norm(xTrue));
            end
        end

        function baselineUsesSuppliedInitialGuess(testCase)
            A = MatrixLinearOp(diag([1,2]));
            b = [1;1];
            x0 = [10;-3];
            for preconditioned = [false,true]
                [~,info] = TestSolverTolerances.solve(A,b,preconditioned, ...
                    InitialGuess=x0,RelativeTolerance=0.2,AbsoluteTolerance=0);
                expected = 0.2*norm(A.adj(b-A.fwd(x0)));
                testCase.verifyEqual(info.stoppingThreshold,expected,'AbsTol',1e-14);
            end
        end

        function relativeTestAdaptsToDataScale(testCase)
            A = MatrixLinearOp(diag([1,2]));
            for preconditioned = [false,true]
                for scale = [1e-4,1,1e4]
                    b = scale*[1;1];
                    [~,info] = TestSolverTolerances.solve(A,b,preconditioned, ...
                        RelativeTolerance=0.4,AbsoluteTolerance=0);
                    testCase.verifyEqual(info.iterations,1);
                    testCase.verifyTrue(info.converged);
                    testCase.verifyEqual(info.stoppingThreshold, ...
                        0.4*norm(A.adj(b)),'RelTol',1e-14);
                end
            end
        end

        function zeroTolerancesPermitAnExactSolve(testCase)
            A = MatrixLinearOp(eye(2));
            for preconditioned = [false,true]
                [x,info] = TestSolverTolerances.solve(A,[1;2],preconditioned, ...
                    RelativeTolerance=0,AbsoluteTolerance=0);
                testCase.verifyEqual(x,[1;2]);
                testCase.verifyEqual(info.iterations,1);
                testCase.verifyEqual(info.stoppingThreshold,0);
                testCase.verifyTrue(info.converged);
                [~,stationary] = TestSolverTolerances.solve(A,[0;0],preconditioned, ...
                    RelativeTolerance=0,AbsoluteTolerance=0);
                testCase.verifyEqual(stationary.iterations,0);
                testCase.verifyTrue(stationary.converged);
            end
        end

        function rejectsNegativeTolerances(testCase)
            A = MatrixLinearOp(eye(2));
            for preconditioned = [false,true]
                testCase.verifyError(@() TestSolverTolerances.solve( ...
                    A,[1;2],preconditioned,RelativeTolerance=-1), ...
                    'MATLAB:validators:mustBeNonnegative');
                testCase.verifyError(@() TestSolverTolerances.solve( ...
                    A,[1;2],preconditioned,AbsoluteTolerance=-1), ...
                    'MATLAB:validators:mustBeNonnegative');
            end
        end
    end

    methods (Static, Access = private)
        function [x,info] = solve(A,b,preconditioned,varargin)
            if preconditioned
                M = MatrixInvLinearOp(eye(A.dimDomain));
                [x,info] = solver_PCGLS(A,b,M,varargin{:});
            else
                [x,info] = solver_CGLS(A,b,varargin{:});
            end
        end
    end
end
