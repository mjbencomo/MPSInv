classdef TestAcousticGreenFunctionSeries2D < matlab.unittest.TestCase
    methods (Test)
        function constructorStoresChannelMetadata(testCase)
            solver = TestAcousticGreenFunctionSeries2D.makeSolver();
            components = { ...
                MultipoleComponent2D([0.4,0.5],[0,0]), ...
                MultipoleComponent2D([0.6,0.5],[1,0])};

            series = AcousticGreenFunctionSeries2D(components,solver);

            testCase.verifyEqual(series.numberOfInputChannels,2);
            testCase.verifyEqual(series.numberOfOutputChannels,1);
            testCase.verifyEqual(series.numberOfTimes,4);
            testCase.verifyEqual(series.components{2},components{2});
            testCase.verifySize(series.values,[7,1,2]);
            testCase.verifyFalse(series.isComputed);
        end

        function rejectsInvalidComponents(testCase)
            solver = TestAcousticGreenFunctionSeries2D.makeSolver();
            component1D = MultipoleComponent1D(0.5,0);

            testCase.verifyError(@() AcousticGreenFunctionSeries2D( ...
                {},solver), ...
                'AcousticGreenFunctionSeries2D:EmptyComponents');
            testCase.verifyError(@() AcousticGreenFunctionSeries2D( ...
                {component1D},solver), ...
                'AcousticGreenFunctionSeries2D:InvalidComponents');
        end

        function computesKernelsAndConstructsOperator(testCase)
            solver = TestAcousticGreenFunctionSeries2D.makeSolver();
            components = { ...
                MultipoleComponent2D([0.4,0.5],[0,0]), ...
                MultipoleComponent2D([0.6,0.5],[1,0])};
            series = AcousticGreenFunctionSeries2D(components,solver);
            series.compute([2,2]);

            K = series.numberOfTimes;
            testCase.verifyTrue(series.isComputed);
            testCase.verifySize(series.values,[2*K-1,1,2]);

            op = series.makeLinearOp();
            testCase.verifyEqual(op.dimDomain,2*K);
            testCase.verifyEqual(op.dimRange,K);
        end
    end

    methods (Static, Access = private)
        function solver = makeSolver()
            computationalGrid = GridSpace2D( ...
                Grid1DUnifPrimal([0,1],11), ...
                Grid1DUnifPrimal([0,1],11));
            outputGrid = GridSpace2D( ...
                Grid1D([0,1],0.5),Grid1D([0,1],0.5));
            solverTimes = [-0.002,0:0.002:0.006];
            medium = AcousticMedium2D.constant( ...
                KappaValue=1,BetaValue=1);
            solver = AcousticSolver2D_PML( ...
                computationalGrid,solverTimes,outputGrid,solverTimes, ...
                medium,PML2D(),EnforceCFL=false);
        end
    end
end
