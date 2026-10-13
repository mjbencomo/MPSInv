classdef TestAcousticSolver2DReceivers < matlab.unittest.TestCase
    properties
        grid
        medium
    end

    methods (TestMethodSetup)
        function createFixtures(testCase)
            testCase.grid = GridSpace2D( ...
                Grid1DUnifPrimal([0,1],21), ...
                Grid1DUnifPrimal([0,1],21));
            testCase.medium = AcousticMedium2D.constant( ...
                KappaValue=1,BetaValue=1);
        end
    end

    methods (Test)
        function samplesThreeFieldsAtPairedLocations(testCase)
            xy = [0.18 0.27;0.64 0.73;0.18 0.27];
            receivers = ReceiverSpace2D(xy);
            times = [0,0.01,0.02];
            solver = AcousticSolver2DStaggered( ...
                testCase.grid,times,receivers,times,testCase.medium);

            p0 = @(x,y) 1+2*x+3*y;
            vx0 = @(x,y) 2*x-3*y;
            vy0 = @(x,y) -x+4*y;
            [P,Vx,Vy] = solver.solve(p0,vx0,vy0);

            testCase.verifySize(P.values,[3,3]);
            testCase.verifySize(Vx.values,[3,3]);
            testCase.verifySize(Vy.values,[3,3]);
            testCase.verifyEqual(P.values(:,1),p0(xy(:,1),xy(:,2)), ...
                AbsTol=1e-12);
            testCase.verifyEqual(Vx.values(:,1),vx0(xy(:,1),xy(:,2)), ...
                AbsTol=1e-12);
            testCase.verifyEqual(Vy.values(:,1),vy0(xy(:,1),xy(:,2)), ...
                AbsTol=1e-12);

            testCase.verifyEqual(P.values(1,:),P.values(3,:));
            testCase.verifyEqual(Vx.values(1,:),Vx.values(3,:));
            testCase.verifyEqual(Vy.values(1,:),Vy.values(3,:));
        end

        function acceptsInteriorReceiversWithoutVelocityExtrapolation(testCase)
            receivers = ReceiverSpace2D([0.1 0.2;0.7 0.8]);
            solver = AcousticSolver2DStaggered( ...
                testCase.grid,[0,0.01],receivers,0, ...
                testCase.medium,VelocityExtrapolation="none");
            testCase.verifyEqual(solver.outputGrid,receivers);
        end

        function velocityBoundaryRequiresExtrapolation(testCase)
            receivers = ReceiverSpace2D([0 0.5]);
            testCase.verifyError(@() AcousticSolver2DStaggered( ...
                testCase.grid,[0,0.01],receivers,0, ...
                testCase.medium,VelocityExtrapolation="none"), ...
                'AcousticSolver2DStaggered:VelocityExtrapolationRequired');
        end

        function rejectsReceiversOutsidePhysicalDomain(testCase)
            receivers = ReceiverSpace2D([0.4 0.5;1.01 0.5]);
            testCase.verifyError(@() AcousticSolver2DStaggered( ...
                testCase.grid,[0,0.01],receivers,0,testCase.medium), ...
                'AcousticSolver2DStaggered:OutputGridOutsideDomain');
        end
    end
end
