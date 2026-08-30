classdef TestAcousticSolver1D < matlab.unittest.TestCase
    % Order-selection tests for the combined staggered solver.

    properties
        grid
        medium
    end

    methods (TestMethodSetup)
        function createFixtures(testCase)
            testCase.grid = GridSpace1D(Grid1DUnifPrimal([0,1],21));
            testCase.medium = AcousticMedium1D.constant( ...
                KappaValue=1,BetaValue=1);
        end
    end

    methods (Test)
        function defaultsToSecondOrder(testCase)
            solver = AcousticSolver1DStaggered( ...
                testCase.grid,[0,0.025],testCase.grid,0, ...
                testCase.medium,EnforceCFL=false);
            testCase.verifyEqual(solver.spatialOrder,2);
        end

        function acceptsFourthOrder(testCase)
            solver = testCase.makeSolver([0,0.025],4,false);
            testCase.verifyEqual(solver.spatialOrder,4);
        end

        function zeroDataRemainZeroForBothOrders(testCase)
            for order = [2,4]
                solver = testCase.makeSolver(0:0.025:0.1,order,true);
                [P,V] = solver.solve(@(x) 0*x,@(x) 0*x);
                testCase.verifyEqual(P.values,zeros(size(P.values)));
                testCase.verifyEqual(V.values,zeros(size(V.values)));
            end
        end

        function usesOrderDependentCFLBound(testCase)
            % dx=0.05 and dt=0.045 give CFL=0.9: stable for order 2,
            % but above the fourth-order limit 6/7.
            times = [0,0.045];
            testCase.verifyWarningFree(@() testCase.makeSolver( ...
                times,2,true));
            testCase.verifyError(@() testCase.makeSolver(times,4,true), ...
                'AcousticSolver1DStaggered:CFLViolation');
        end

        function fourthOrderRequiresFivePressurePoints(testCase)
            smallGrid = GridSpace1D(Grid1DUnifPrimal([0,1],4));
            testCase.verifyError(@() AcousticSolver1DStaggered( ...
                smallGrid,[0,0.01],smallGrid,0,testCase.medium, ...
                SpatialOrder=4), ...
                'AcousticSolver1DStaggered:GridTooSmall');
        end

        function solutionObeysSourceSuperposition(testCase)
            times = 0:0.01:0.05;
            solver = testCase.makeSolver(times,2,true);
            source1 = AcousticSource1D(@(x) sin(pi*x),@(t) 1+t);
            source2 = AcousticSource1D(@(x) x.*(1-x),@(t) cos(t));
            combined = AcousticSource1D.combine({source1,source2});

            [P1,V1] = solver.solve(@(x) 0*x,@(x) 0*x,Source=source1);
            [P2,V2] = solver.solve(@(x) 0*x,@(x) 0*x,Source=source2);
            [P,V] = solver.solve(@(x) 0*x,@(x) 0*x,Source=combined);

            testCase.verifyEqual(P.values,P1.values+P2.values,AbsTol=1e-13);
            testCase.verifyEqual(V.values,V1.values+V2.values,AbsTol=1e-13);
        end

        function acceptsLocalizedPressureAndVelocityTerms(testCase)
            times = 0:0.01:0.05;
            solver = testCase.makeSolver(times,2,true);
            dualGrid = GridSpace1D(Grid1DUnifDual(testCase.grid.x));
            source = AcousticSource1D.zero();
            source.addLocalizedPressureTerm( ...
                testCase.grid,11,2,@(t) 1);
            source.addLocalizedVelocityTerm( ...
                dualGrid,10,-1,@(t) 1);

            [P,V] = solver.solve(@(x) 0*x,@(x) 0*x,Source=source);

            testCase.verifyGreaterThan(max(abs(P.values),[],'all'),0);
            testCase.verifyGreaterThan(max(abs(V.values),[],'all'),0);
        end
    end

    methods (Access = private)
        function solver = makeSolver(testCase,times,order,enforceCFL)
            solver = AcousticSolver1DStaggered( ...
                testCase.grid,times,testCase.grid,times(end), ...
                testCase.medium,SpatialOrder=order, ...
                EnforceCFL=enforceCFL);
        end
    end
end
