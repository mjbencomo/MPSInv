classdef TestAcousticGreenFunction1D < matlab.unittest.TestCase
    methods (Test)
        function rejectsNonAcousticSolver(testCase)
            component = MultipoleComponent1D(0.5,0);
            testCase.verifyError(@() AcousticGreenFunction1D( ...
                component,struct()), ...
                'AcousticGreenFunction1D:InvalidSolver');
        end

        function constructorStoresMetadataAndLagGrid(testCase)
            component = MultipoleComponent1D(0.4,1, ...
                TargetField="velocity");
            solver = TestAcousticGreenFunction1D.makeSolver( ...
                [-0.01,0:0.01:0.04],[0.2,0.7]);

            G = AcousticGreenFunction1D(component,solver, ...
                OutputField="velocity");

            testCase.verifyTrue(isa(G,'GreenFunction'));
            testCase.verifyEqual(G.component,component);
            testCase.verifyEqual(G.solver,solver);
            testCase.verifyEqual(G.outputGrid,solver.outputGrid);
            testCase.verifyEqual(G.outputField,"velocity");
            testCase.verifyEqual(G.dt,0.01,'AbsTol',10*eps);
            testCase.verifyEqual(G.numberOfTimes,5);
            testCase.verifyEqual(G.numberOfOutputChannels,2);
            testCase.verifyEqual(G.zeroLagIndex,5);
            testCase.verifyEqual(G.lags,-0.04:0.01:0.04, ...
                'AbsTol',10*eps);
            testCase.verifySize(G.values,[9,2]);
            testCase.verifyFalse(G.isComputed);
        end

        function rejectsInvalidSolverTimes(testCase)
            component = MultipoleComponent1D(0.5,0);
            solver = TestAcousticGreenFunction1D.makeSolver(0:0.01:0.03,0.5);
            testCase.verifyError(@() AcousticGreenFunction1D(component,solver), ...
                'GreenFunction:InvalidSolverTimes');
        end

        function requiresComputationBeforeUsingKernel(testCase)
            component = MultipoleComponent1D(0.5,0);
            solver = TestAcousticGreenFunction1D.makeSolver( ...
                [-0.01,0:0.01:0.03],0.5);
            G = AcousticGreenFunction1D(component,solver);

            testCase.verifyError(@() G.kernelAt(1), ...
                'GreenFunction:NotComputed');
            testCase.verifyError(@() G.makeLinearOp(), ...
                'GreenFunction:NotComputed');
        end

        function computesPressureGreenFunction(testCase)
            [G,~] = TestAcousticGreenFunction1D.computedExample();

            K = G.numberOfTimes;
            testCase.verifyTrue(G.isComputed);
            testCase.verifySize(G.values,[2*K-1,1]);
            testCase.verifyEqual(G.values(1:K-2,:),zeros(K-2,1));
            testCase.verifyGreaterThan(max(abs(G.values),[],'all'),0);
            testCase.verifyEqual(G.kernelAt(1),G.values(:,1));
        end

        function constructsConsistentLinearOperator(testCase)
            [G,~] = TestAcousticGreenFunction1D.computedExample();
            op = G.makeLinearOp();
            w = (1:G.numberOfTimes).';

            expected = G.dt*convDisc(G.values(:,1),w);

            testCase.verifyTrue(isa(op,'ConvMultiOp'));
            testCase.verifyEqual(op.dimDomain,G.numberOfTimes);
            testCase.verifyEqual(op.dimRange,G.numberOfTimes);
            testCase.verifyEqual(op.fwd(w),expected,'AbsTol',1e-10);
        end
    end

    methods (Static, Access = private)
        function [G,computationalGrid] = computedExample()
            primal = Grid1DUnifPrimal([0,1],21);
            computationalGrid = GridSpace1D(primal);
            outputGrid = GridSpace1D(Grid1D([0,1],0.5));
            component = MultipoleComponent1D(0.5,0);
            medium = AcousticMedium1D.constant( ...
                KappaValue=1,BetaValue=1);
            pml = PML1D();

            solverTimes = [-0.002,0:0.002:0.008];
            solver = AcousticSolver1D_PML( ...
                computationalGrid,solverTimes,outputGrid,solverTimes, ...
                medium,pml,SpatialOrder=2,EnforceCFL=false);
            G = AcousticGreenFunction1D(component,solver);
            G.compute(2);
        end

        function solver = makeSolver(solverTimes,receiverLocations)
            computationalGrid = GridSpace1D( ...
                Grid1DUnifPrimal([0,1],21));
            if isscalar(receiverLocations)
                receiverGrid = Grid1D([0,1],receiverLocations);
            else
                receiverGrid = Grid1DUnifPrimal( ...
                    [receiverLocations(1),receiverLocations(end)], ...
                    numel(receiverLocations));
            end
            outputGrid = GridSpace1D(receiverGrid);
            medium = AcousticMedium1D.constant( ...
                KappaValue=1,BetaValue=1);
            solver = AcousticSolver1D_PML( ...
                computationalGrid,solverTimes,outputGrid,solverTimes, ...
                medium,PML1D(),EnforceCFL=false);
        end
    end
end
