classdef TestAcousticSolver2D < matlab.unittest.TestCase
    % Order, source, boundary, and storage tests for the 2D solver.

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
        function defaultsToSecondOrder(testCase)
            solver = AcousticSolver2DStaggered( ...
                testCase.grid,[0,0.01],testCase.grid,0, ...
                testCase.medium,EnforceCFL=false);
            testCase.verifyEqual(solver.spatialOrder,2);
        end

        function acceptsFourthOrder(testCase)
            solver = testCase.makeSolver([0,0.01],4,false);
            testCase.verifyEqual(solver.spatialOrder,4);
        end

        function zeroDataRemainZeroForBothOrders(testCase)
            zero = @(x,y) zeros(size(x));
            for order = [2,4]
                solver = testCase.makeSolver(0:0.01:0.04,order,true);
                [P,Vx,Vy] = solver.solve(zero,zero,zero);
                testCase.verifyEqual(P.values,zeros(size(P.values)));
                testCase.verifyEqual(Vx.values,zeros(size(Vx.values)));
                testCase.verifyEqual(Vy.values,zeros(size(Vy.values)));
            end
        end

        function returnsExpectedTwoDimensionalStorage(testCase)
            solver = testCase.makeSolver(0:0.01:0.04,2,true);
            zero = @(x,y) zeros(size(x));
            [P,Vx,Vy] = solver.solve(zero,zero,zero);
            expectedSize = testCase.grid.N;

            testCase.verifySize(P.values,expectedSize);
            testCase.verifySize(Vx.values,expectedSize);
            testCase.verifySize(Vy.values,expectedSize);
        end

        function usesOrderDependentMultidimensionalCFLBound(testCase)
            % dx=dy=0.05 and dt=0.032 give CFL about 0.905.
            times = [0,0.032];
            testCase.verifyWarningFree(@() testCase.makeSolver(times,2,true));
            testCase.verifyError(@() testCase.makeSolver(times,4,true), ...
                'AcousticSolver2DStaggered:CFLViolation');
        end

        function fourthOrderRequiresFivePointsInEachDirection(testCase)
            smallGrid = GridSpace2D( ...
                Grid1DUnifPrimal([0,1],4), ...
                Grid1DUnifPrimal([0,1],9));
            testCase.verifyError(@() AcousticSolver2DStaggered( ...
                smallGrid,[0,0.01],smallGrid,0,testCase.medium, ...
                SpatialOrder=4), ...
                'AcousticSolver2DStaggered:GridTooSmall');
        end

        function pressureBoundaryRemainsZero(testCase)
            solver = testCase.makeSolver(0:0.01:0.04,4,true);
            p0 = @(x,y) sin(pi*x).*sin(pi*y);
            zero = @(x,y) zeros(size(x));
            [P,~,~] = solver.solve(p0,zero,zero);

            testCase.verifyEqual(P.values([1,end],:,:), ...
                zeros(size(P.values([1,end],:,:))),AbsTol=1e-14);
            testCase.verifyEqual(P.values(:,[1,end],:), ...
                zeros(size(P.values(:,[1,end],:))),AbsTol=1e-14);
        end

        function solutionObeysSourceSuperposition(testCase)
            times = 0:0.01:0.04;
            solver = testCase.makeSolver(times,2,true);
            source1 = AcousticSource2D( ...
                @(x,y) sin(pi*x).*sin(pi*y),@(t) 1+t);
            source2 = AcousticSource2D( ...
                @(x,y) x.*(1-x).*y.*(1-y),@(t) cos(t));
            combined = AcousticSource2D.combine({source1,source2});
            zero = @(x,y) zeros(size(x));

            [P1,Vx1,Vy1] = solver.solve(zero,zero,zero,Source=source1);
            [P2,Vx2,Vy2] = solver.solve(zero,zero,zero,Source=source2);
            [P,Vx,Vy] = solver.solve(zero,zero,zero,Source=combined);

            testCase.verifyEqual(P.values,P1.values+P2.values,AbsTol=1e-13);
            testCase.verifyEqual(Vx.values,Vx1.values+Vx2.values,AbsTol=1e-13);
            testCase.verifyEqual(Vy.values,Vy1.values+Vy2.values,AbsTol=1e-13);
        end

        function acceptsLocalizedTermsOnThreeStaggeredGrids(testCase)
            solver = testCase.makeSolver(0:0.01:0.04,2,true);
            [~,Bx,By] = sampleMedium2D(testCase.grid,testCase.medium);
            source = AcousticSource2D.zero();
            center = sub2ind(testCase.grid.N,11,11);
            source.addLocalizedPressureTerm(testCase.grid,center,2,@(t) 1);
            source.addLocalizedVelocityXTerm(Bx.grid,center,-1,@(t) 1);
            source.addLocalizedVelocityYTerm(By.grid,center,1,@(t) 1);
            zero = @(x,y) zeros(size(x));

            [P,Vx,Vy] = solver.solve(zero,zero,zero,Source=source);

            testCase.verifyGreaterThan(max(abs(P.values),[],'all'),0);
            testCase.verifyGreaterThan(max(abs(Vx.values),[],'all'),0);
            testCase.verifyGreaterThan(max(abs(Vy.values),[],'all'),0);
        end
    end

    methods (Access = private)
        function solver = makeSolver(testCase,times,order,enforceCFL)
            solver = AcousticSolver2DStaggered( ...
                testCase.grid,times,testCase.grid,times(end), ...
                testCase.medium,SpatialOrder=order, ...
                EnforceCFL=enforceCFL);
        end
    end
end
