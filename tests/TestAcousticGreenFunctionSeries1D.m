classdef TestAcousticGreenFunctionSeries1D < matlab.unittest.TestCase
    methods (Test)
        function constructorStoresChannelMetadata(testCase)
            solver = TestAcousticGreenFunctionSeries1D.makeSolver();
            components = { ...
                MultipoleComponent1D(0.35,0), ...
                MultipoleComponent1D(0.65,1)};

            series = AcousticGreenFunctionSeries1D( ...
                components,solver,OutputField="velocity");

            testCase.verifyEqual(series.numberOfInputChannels,2);
            testCase.verifyEqual(series.numberOfOutputChannels,2);
            testCase.verifyEqual(series.numberOfTimes,5);
            testCase.verifyEqual(series.outputField,"velocity");
            testCase.verifyEqual(series.solver,solver);
            testCase.verifyEqual(series.components{1},components{1});
            testCase.verifyEqual(series.components{2},components{2});
            testCase.verifySize(series.values,[9,2,2]);
            testCase.verifyFalse(series.isComputed);
        end

        function rejectsInvalidComponents(testCase)
            solver = TestAcousticGreenFunctionSeries1D.makeSolver();
            component2D = MultipoleComponent2D([0.5,0.5],[0,0]);

            testCase.verifyError(@() AcousticGreenFunctionSeries1D( ...
                {},solver), ...
                'AcousticGreenFunctionSeries1D:EmptyComponents');
            testCase.verifyError(@() AcousticGreenFunctionSeries1D( ...
                {component2D},solver), ...
                'AcousticGreenFunctionSeries1D:InvalidComponents');
        end

        function validatesApproximationOrderCount(testCase)
            solver = TestAcousticGreenFunctionSeries1D.makeSolver();
            components = { ...
                MultipoleComponent1D(0.35,0), ...
                MultipoleComponent1D(0.65,1)};
            series = AcousticGreenFunctionSeries1D(components,solver);

            testCase.verifyError(@() series.compute([2,2,2]), ...
                'GreenFunctionSeries:ApproximationOrderSizeMismatch');
        end

        function computesKernelsAndConstructsOperator(testCase)
            solver = TestAcousticGreenFunctionSeries1D.makeSolver();
            components = { ...
                MultipoleComponent1D(0.35,0), ...
                MultipoleComponent1D(0.65,1)};
            series = AcousticGreenFunctionSeries1D(components,solver);

            testCase.verifyError(@() series.makeLinearOp(), ...
                'GreenFunctionSeries:NotComputed');
            series.compute(2);

            K = series.numberOfTimes;
            N = series.numberOfOutputChannels;
            M = series.numberOfInputChannels;
            testCase.verifyTrue(series.isComputed);
            testCase.verifySize(series.values,[2*K-1,N,M]);
            testCase.verifyEqual(series.kernelAt(2,1), ...
                series.values(:,2,1));

            op = series.makeLinearOp();
            W = reshape(1:K*M,K,M);
            expected = series.dt*convMulti(series.values,W);

            testCase.verifyEqual(op.dimDomain,K*M);
            testCase.verifyEqual(op.dimRange,K*N);
            testCase.verifyEqual(op.fwd(W(:)),expected(:), ...
                'AbsTol',1e-10);

            rng(4)
            x = randn(op.dimDomain,1);
            y = randn(op.dimRange,1);
            testCase.verifyEqual(dot(op.fwd(x),y),dot(x,op.adj(y)), ...
                'RelTol',1e-11,'AbsTol',1e-11);
        end
    end

    methods (Static, Access = private)
        function solver = makeSolver()
            computationalGrid = GridSpace1D( ...
                Grid1DUnifPrimal([0,1],31));
            receiverGrid = Grid1DUnifPrimal([0.25,0.75],2);
            outputGrid = GridSpace1D(receiverGrid);
            solverTimes = [-0.002,0:0.002:0.008];
            medium = AcousticMedium1D.constant( ...
                KappaValue=1,BetaValue=1);
            solver = AcousticSolver1D_PML( ...
                computationalGrid,solverTimes,outputGrid,solverTimes, ...
                medium,PML1D(),EnforceCFL=false);
        end
    end
end
