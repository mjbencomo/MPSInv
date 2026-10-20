classdef TestOtherAcousticSolverReceivers < matlab.unittest.TestCase
    properties
        grid1D
        grid2D
        medium1D
        medium2D
        pml1D
        pml2D
    end

    methods (TestMethodSetup)
        function createFixtures(testCase)
            x = Grid1DUnifPrimal([0,1],21);
            testCase.grid1D = GridSpace1D(x);
            testCase.grid2D = GridSpace2D(x,x);
            testCase.medium1D = AcousticMedium1D.constant( ...
                KappaValue=1,BetaValue=1);
            testCase.medium2D = AcousticMedium2D.constant( ...
                KappaValue=1,BetaValue=1);
            testCase.pml1D = PML1D(WidthLeft=0.2,WidthRight=0.2, ...
                ReflectionCoefficient=1e-6,PolynomialDegree=3);
            testCase.pml2D = PML2D(WidthLeft=0.2,WidthRight=0.2, ...
                WidthBottom=0.2,WidthTop=0.2, ...
                ReflectionCoefficient=1e-6,PolynomialDegree=3);
        end
    end

    methods (Test)
        function oneDimensionalStaggeredSamplesOrderedReceivers(testCase)
            x = [0.63;0.18;0.63];
            receivers = ReceiverSpace1D(x);
            times = [0,0.01,0.02];
            solver = AcousticSolver1DStaggered( ...
                testCase.grid1D,times,receivers,times,testCase.medium1D);
            p0 = @(x) 1+2*x;
            v0 = @(x) 2-3*x;

            [P,V] = solver.solve(p0,v0);
            testCase.verifySize(P.values,[3,3]);
            testCase.verifySize(V.values,[3,3]);
            testCase.verifyEqual(P.values(:,1),p0(x),AbsTol=1e-12);
            testCase.verifyEqual(V.values(:,1),v0(x),AbsTol=1e-12);
            testCase.verifyEqual(P.values(1,:),P.values(3,:));
            testCase.verifyEqual(V.values(1,:),V.values(3,:));
        end

        function oneDimensionalStaggeredChecksVelocityExtent(testCase)
            receivers = ReceiverSpace1D([0.5;0]);
            testCase.verifyError(@() AcousticSolver1DStaggered( ...
                testCase.grid1D,[0,0.01],receivers,0, ...
                testCase.medium1D,VelocityExtrapolation="none"), ...
                'AcousticSolver1DStaggered:VelocityExtrapolationRequired');
            testCase.verifyError(@() AcousticSolver1DStaggered( ...
                testCase.grid1D,[0,0.01],ReceiverSpace1D(1.01),0, ...
                testCase.medium1D), ...
                'AcousticSolver1DStaggered:OutputGridOutsideDomain');
        end

        function oneDimensionalPMLSamplesPhysicalBoundary(testCase)
            x = [1;0.43;0];
            receivers = ReceiverSpace1D(x);
            times = [0,0.01];
            solver = AcousticSolver1D_PML( ...
                testCase.grid1D,times,receivers,times, ...
                testCase.medium1D,testCase.pml1D, ...
                VelocityExtrapolation="none");
            p0 = @(x) 1+2*x;
            v0 = @(x) 2-3*x;

            [P,V] = solver.solve(p0,v0);
            testCase.verifySize(P.values,[3,2]);
            testCase.verifySize(V.values,[3,2]);
            testCase.verifyEqual(P.values(:,1),p0(x),AbsTol=1e-12);
            testCase.verifyEqual(V.values(:,1),v0(x),AbsTol=1e-12);
            testCase.verifyError(@() AcousticSolver1D_PML( ...
                testCase.grid1D,times,ReceiverSpace1D(-0.01),0, ...
                testCase.medium1D,testCase.pml1D), ...
                'AcousticSolver1D_PML:OutputGridOutsideDomain');
        end

        function twoDimensionalPMLSamplesPairedReceivers(testCase)
            xy = [0 0.4;0.61 0.73;0 0.4];
            receivers = ReceiverSpace2D(xy);
            times = [0,0.01];
            solver = AcousticSolver2D_PML( ...
                testCase.grid2D,times,receivers,times, ...
                testCase.medium2D,testCase.pml2D, ...
                VelocityExtrapolation="none");
            p0 = @(x,y) 1+2*x+3*y;
            vx0 = @(x,y) 2*x-3*y;
            vy0 = @(x,y) -x+4*y;

            [P,Vx,Vy] = solver.solve(p0,vx0,vy0);
            testCase.verifySize(P.values,[3,2]);
            testCase.verifySize(Vx.values,[3,2]);
            testCase.verifySize(Vy.values,[3,2]);
            testCase.verifyEqual(P.values(:,1),p0(xy(:,1),xy(:,2)), ...
                AbsTol=1e-12);
            testCase.verifyEqual(Vx.values(:,1),vx0(xy(:,1),xy(:,2)), ...
                AbsTol=1e-12);
            testCase.verifyEqual(Vy.values(:,1),vy0(xy(:,1),xy(:,2)), ...
                AbsTol=1e-12);
            testCase.verifyEqual(P.values(1,:),P.values(3,:));
            testCase.verifyError(@() AcousticSolver2D_PML( ...
                testCase.grid2D,times,ReceiverSpace2D([0.5 1.01]),0, ...
                testCase.medium2D,testCase.pml2D), ...
                'AcousticSolver2D_PML:OutputGridOutsideDomain');
        end
    end
end
