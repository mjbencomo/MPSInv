classdef TestAcousticGreenFunction2D < matlab.unittest.TestCase
    methods (Test)
        function rejectsNonAcousticSolver(testCase)
            component = MultipoleComponent2D([0.5,0.5],[0,0]);
            testCase.verifyError(@() AcousticGreenFunction2D( ...
                component,struct()), ...
                'AcousticGreenFunction2D:InvalidSolver');
        end

        function constructorStoresTwoDimensionalMetadata(testCase)
            component = MultipoleComponent2D( ...
                [0.4,0.6],[1,0],TargetField="velocityX");
            outputGrid = GridSpace2D( ...
                Grid1D([0,1],0.3),Grid1D([0,1],0.7));
            computationalGrid = GridSpace2D( ...
                Grid1DUnifPrimal([0,1],11), ...
                Grid1DUnifPrimal([0,1],11));
            solverTimes = [-0.01,0:0.01:0.03];
            solver = AcousticSolver2D_PML( ...
                computationalGrid,solverTimes,outputGrid,solverTimes, ...
                AcousticMedium2D.constant(KappaValue=1,BetaValue=1), ...
                PML2D(),EnforceCFL=false);

            G = AcousticGreenFunction2D(component,solver, ...
                OutputField="velocityY");

            testCase.verifyTrue(isa(G,'GreenFunction'));
            testCase.verifyEqual(G.component,component);
            testCase.verifyEqual(G.outputField,"velocityY");
            testCase.verifyEqual(G.numberOfOutputChannels,1);
            testCase.verifySize(G.values,[7,1]);
        end

        function computesPressureGreenFunction(testCase)
            xPrimal = Grid1DUnifPrimal([0,1],11);
            yPrimal = Grid1DUnifPrimal([0,1],11);
            computationalGrid = GridSpace2D(xPrimal,yPrimal);
            outputGrid = GridSpace2D( ...
                Grid1D([0,1],0.5),Grid1D([0,1],0.5));
            component = MultipoleComponent2D([0.5,0.5],[0,0]);
            medium = AcousticMedium2D.constant( ...
                KappaValue=1,BetaValue=1);
            pml = PML2D();

            solverTimes = [-0.002,0:0.002:0.006];
            solver = AcousticSolver2D_PML( ...
                computationalGrid,solverTimes,outputGrid,solverTimes, ...
                medium,pml,SpatialOrder=2,EnforceCFL=false);
            G = AcousticGreenFunction2D(component,solver);
            G.compute(2);

            K = G.numberOfTimes;
            testCase.verifyTrue(G.isComputed);
            testCase.verifySize(G.values,[2*K-1,1]);
            testCase.verifyEqual(G.values(1:K-2,:),zeros(K-2,1));
            testCase.verifyGreaterThan(max(abs(G.values),[],'all'),0);

            op = G.makeLinearOp();
            testCase.verifyEqual(op.dimDomain,K);
            testCase.verifyEqual(op.dimRange,K);
        end
    end
end
