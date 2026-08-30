classdef TestMultipoleSeries1D < matlab.unittest.TestCase
    properties
        PressureGrid
        VelocityGrid
    end

    methods (TestMethodSetup)
        function createGrids(testCase)
            primal = Grid1DUnifPrimal([0,1],41);
            testCase.PressureGrid = GridSpace1D(primal);
            testCase.VelocityGrid = GridSpace1D(Grid1DUnifDual(primal));
        end
    end

    methods (Test)
        function mpsAppxSatisfiesMomentConditions(testCase)
            x = testCase.PressureGrid.x.pts;
            h = testCase.PressureGrid.h;
            xc = 0.437;

            for q = [2,4]
                for s = 0:2
                    [values,indices,weights] = MPSappx( ...
                        testCase.PressureGrid,xc,q,s);

                    testCase.verifyNumElements(indices,q+s);
                    testCase.verifyNumElements(weights,q+s);
                    testCase.verifyEqual(values(indices),weights);

                    for k = 0:q+s-1
                        actual = h*sum(weights.*(x(indices)-xc).^k);
                        expected = 0;
                        if k == s
                            expected = (-1)^s*factorial(s);
                        end
                        tolerance = 1e4*eps*max(1,abs(expected));
                        testCase.verifyEqual(actual,expected, ...
                            AbsTol=tolerance);
                    end
                end
            end
        end

        function mpsAppxPreservesRowOrientation(testCase)
            x = testCase.PressureGrid.x.pts.';
            values = MPSappx(x,0.5,4,0);
            testCase.verifySize(values,size(x));
        end

        function termEvaluatesScaledTimeFunction(testCase)
            term = MultipoleTerm1D(0.4,1,@(t) 1+t, ...
                Amplitude=-2,TargetField="velocity");

            testCase.verifyEqual(term.location,0.4);
            testCase.verifyEqual(term.derivativeOrder,1);
            testCase.verifyEqual(term.targetField,"velocity");
            testCase.verifyEqual(term.evaluateTime(0.5),-3);
            testCase.verifyTrue(isa(term,'MultipoleTerm'));
        end

        function singleTermBuildsAcousticSource(testCase)
            term = MultipoleTerm1D(0.45,1,@(t) 1+t);
            source = term.discretize( ...
                testCase.PressureGrid,testCase.VelocityGrid, ...
                ApproximationOrder=4);

            testCase.verifyTrue(isa(source,'AcousticSource'));
            testCase.verifyTrue(isa(source,'AcousticSource1D'));
            testCase.verifyEqual(source.numberOfPressureTerms,1);
        end

        function discretizesPressureAndVelocityTerms(testCase)
            pressureTerm = MultipoleTerm1D(0.35,0,@(t) 1+t);
            velocityTerm = MultipoleTerm1D(0.65,1,@(t) 2-t, ...
                TargetField="velocity");
            mps = MultipoleSeries1D( ...
                [pressureTerm,velocityTerm],ApproximationOrders=4);

            source = mps.discretize( ...
                testCase.PressureGrid,testCase.VelocityGrid);
            preparedP = source.preparePressure(testCase.PressureGrid);
            preparedV = source.prepareVelocity(testCase.VelocityGrid);

            testCase.verifyTrue(issparse(preparedP.spatialWeights));
            testCase.verifyTrue(issparse(preparedV.spatialWeights));
            testCase.verifyEqual(source.numberOfPressureTerms,1);
            testCase.verifyEqual(source.numberOfVelocityTerms,1);
            testCase.verifyGreaterThan( ...
                norm(source.evaluatePrepared(preparedP,0.2)),0);
            testCase.verifyGreaterThan( ...
                norm(source.evaluatePrepared(preparedV,0.2)),0);
        end

        function storesLocationsAndPerTermOrders(testCase)
            terms = [ ...
                MultipoleTerm1D(0.35,0,@(t) 1), ...
                MultipoleTerm1D(0.65,1,@(t) 1,TargetField="velocity")];
            series = MultipoleSeries1D( ...
                terms,ApproximationOrders=[2,4]);

            testCase.verifyEqual(series.numberOfTerms,2);
            testCase.verifyEqual(series.locations,[0.35;0.65]);
            testCase.verifyEqual(series.approximationOrders,[2,4]);
        end

        function rejectsPressureStencilAtBoundary(testCase)
            term = MultipoleTerm1D(0,0,@(t) 1);
            mps = MultipoleSeries1D(term,ApproximationOrders=4);

            testCase.verifyError(@() mps.discretize( ...
                testCase.PressureGrid,testCase.VelocityGrid), ...
                'MultipoleTerm1D:PressureStencilAtBoundary');
        end

        function multipoleSourceRunsInSolver(testCase)
            medium = AcousticMedium1D.constant( ...
                KappaValue=1,BetaValue=1);
            times = 0:0.005:0.025;
            solver = AcousticSolver1DStaggered( ...
                testCase.PressureGrid,times, ...
                testCase.PressureGrid,times(end),medium, ...
                EnforceCFL=false);

            term = MultipoleTerm1D(0.5,0,@(t) 1);
            mps = MultipoleSeries1D(term,ApproximationOrders=2);
            source = mps.discretize( ...
                testCase.PressureGrid,testCase.VelocityGrid);

            [P,V] = solver.solve(@(x) 0*x,@(x) 0*x,Source=source);
            testCase.verifyGreaterThan(max(abs(P.values),[],'all'),0);
            testCase.verifyGreaterThan(max(abs(V.values),[],'all'),0);
        end
    end
end
