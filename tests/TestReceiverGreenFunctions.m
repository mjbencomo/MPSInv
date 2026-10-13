classdef TestReceiverGreenFunctions < matlab.unittest.TestCase
    methods (Test)
        function oneDimensionalKernelsFollowReceiverOrder(testCase)
            grid = GridSpace1D(Grid1DUnifPrimal([0,1],21));
            receivers = ReceiverSpace1D([0.6;0.5;0.6]);
            times = [-0.002,0:0.002:0.006];
            medium = AcousticMedium1D.constant( ...
                KappaValue=1,BetaValue=1);
            solver = AcousticSolver1D_PML( ...
                grid,times,receivers,times,medium,PML1D(), ...
                EnforceCFL=false);
            component = MultipoleComponent1D(0.5,0);
            G = AcousticGreenFunction1D(component,solver);
            K = G.numberOfTimes;

            testCase.verifyEqual(G.numberOfOutputChannels,3);
            testCase.verifySize(G.values,[2*K-1,3]);
            G.compute(2);

            testCase.verifyTrue(G.isComputed);
            testCase.verifyGreaterThan(max(abs(G.values),[],'all'),0);
            testCase.verifyEqual(G.values(:,1),G.values(:,3));
            testCase.verifyEqual(G.kernelAt(2),G.values(:,2));

            op = G.makeLinearOp();
            testCase.verifyEqual(op.dimDomain,K);
            testCase.verifyEqual(op.dimRange,3*K);
        end

        function twoDimensionalSeriesHasOneChannelPerPairedReceiver(testCase)
            grid = GridSpace2D( ...
                Grid1DUnifPrimal([0,1],11), ...
                Grid1DUnifPrimal([0,1],11));
            receivers = ReceiverSpace2D( ...
                [0.5 0.5;0.6 0.5;0.5 0.5]);
            times = [-0.002,0:0.002:0.006];
            medium = AcousticMedium2D.constant( ...
                KappaValue=1,BetaValue=1);
            solver = AcousticSolver2D_PML( ...
                grid,times,receivers,times,medium,PML2D(), ...
                EnforceCFL=false);
            components = { ...
                MultipoleComponent2D([0.5,0.5],[0,0]), ...
                MultipoleComponent2D([0.4,0.5],[0,0])};
            series = AcousticGreenFunctionSeries2D(components,solver);
            K = series.numberOfTimes;

            testCase.verifyEqual(series.numberOfInputChannels,2);
            testCase.verifyEqual(series.numberOfOutputChannels,3);
            testCase.verifySize(series.values,[2*K-1,3,2]);
            series.compute(2);

            testCase.verifyTrue(series.isComputed);
            testCase.verifyGreaterThan(max(abs(series.values),[],'all'),0);
            testCase.verifyEqual(series.values(:,1,:),series.values(:,3,:));
            testCase.verifyEqual(series.kernelAt(2,1), ...
                series.values(:,2,1));

            op = series.makeLinearOp();
            testCase.verifyEqual(op.dimDomain,2*K);
            testCase.verifyEqual(op.dimRange,3*K);
        end
    end
end
