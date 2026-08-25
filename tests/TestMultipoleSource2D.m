classdef TestMultipoleSource2D < matlab.unittest.TestCase
    % Tests for 2D multipole terms and tensor-product discretization.

    properties
        PressureGrid
        VelocityXGrid
        VelocityYGrid
    end

    methods (TestMethodSetup)
        function createGrids(testCase)
            xPrimal = Grid1DUnifPrimal([0,1],41);
            yPrimal = Grid1DUnifPrimal([0,1],41);
            testCase.PressureGrid = GridSpace2D(xPrimal,yPrimal);
            testCase.VelocityXGrid = GridSpace2D( ...
                Grid1DUnifDual(xPrimal),yPrimal);
            testCase.VelocityYGrid = GridSpace2D( ...
                xPrimal,Grid1DUnifDual(yPrimal));
        end
    end

    methods (Test)
        function termEvaluatesScaledTimeFunction(testCase)
            term = MultipoleTerm2D([0.4,0.6],[1,2],@(t) 1+t, ...
                Amplitude=-2,TargetField="velocityY");

            testCase.verifyEqual(term.location,[0.4,0.6]);
            testCase.verifyEqual(term.derivativeOrder,[1,2]);
            testCase.verifyEqual(term.targetField,"velocityY");
            testCase.verifyEqual(term.evaluateTime(0.5),-3);
        end

        function rejectsInvalidTimeOutput(testCase)
            term = MultipoleTerm2D([0.5,0.5],[0,0],@(t) [1,2]);
            testCase.verifyError(@() term.evaluateTime(0), ...
                'MultipoleTerm2D:InvalidTimeOutput');
        end

        function tensorStencilSatisfiesScaledMoments(testCase)
            q = 4;
            derivativeOrder = [1,2];
            location = [0.437,0.463];
            term = MultipoleTerm2D(location,derivativeOrder,@(t) 1);
            mps = MultipoleSource2D(term,ApproximationOrder=q);
            source = mps.discretize(testCase.PressureGrid, ...
                testCase.VelocityXGrid,testCase.VelocityYGrid);
            prepared = source.preparePressure(testCase.PressureGrid);

            weights = full(prepared.spatialWeights(:,1));
            indices = find(weights ~= 0);
            weights = weights(indices);
            [ix,iy] = ind2sub(testCase.PressureGrid.N,indices);
            x = testCase.PressureGrid.x.pts(ix);
            y = testCase.PressureGrid.y.pts(iy);
            hx = testCase.PressureGrid.h(1);
            hy = testCase.PressureGrid.h(2);
            rx = (x-location(1))/hx;
            ry = (y-location(2))/hy;
            scaledWeights = weights* ...
                hx^(derivativeOrder(1)+1)* ...
                hy^(derivativeOrder(2)+1);

            for kx = 0:q+derivativeOrder(1)-1
                for ky = 0:q+derivativeOrder(2)-1
                    actual = sum(scaledWeights.*rx.^kx.*ry.^ky);
                    expectedX = 0;
                    expectedY = 0;
                    if kx == derivativeOrder(1)
                        expectedX = (-1)^derivativeOrder(1)* ...
                            factorial(derivativeOrder(1));
                    end
                    if ky == derivativeOrder(2)
                        expectedY = (-1)^derivativeOrder(2)* ...
                            factorial(derivativeOrder(2));
                    end
                    expected = expectedX*expectedY;
                    testCase.verifyEqual(actual,expected, ...
                        AbsTol=1e5*eps*max(1,abs(expected)));
                end
            end
        end

        function mixedStencilHasExpectedSupportSize(testCase)
            q = 4;
            derivativeOrder = [1,2];
            term = MultipoleTerm2D([0.453,0.557], ...
                derivativeOrder,@(t) 1);
            mps = MultipoleSource2D(term,ApproximationOrder=q);
            source = mps.discretize( ...
                testCase.PressureGrid,testCase.VelocityXGrid, ...
                testCase.VelocityYGrid);
            prepared = source.preparePressure(testCase.PressureGrid);

            expectedSupport = prod(q+derivativeOrder);
            testCase.verifyEqual(nnz(prepared.spatialWeights),expectedSupport);
        end

        function discretizesAllThreeTargetFields(testCase)
            pressureTerm = MultipoleTerm2D( ...
                [0.35,0.45],[0,0],@(t) 1+t);
            velocityXTerm = MultipoleTerm2D( ...
                [0.55,0.45],[1,0],@(t) 2-t, ...
                TargetField="velocityX");
            velocityYTerm = MultipoleTerm2D( ...
                [0.45,0.55],[0,1],@(t) cos(t), ...
                TargetField="velocityY");
            mps = MultipoleSource2D( ...
                [pressureTerm,velocityXTerm,velocityYTerm], ...
                ApproximationOrder=4);

            source = mps.discretize(testCase.PressureGrid, ...
                testCase.VelocityXGrid,testCase.VelocityYGrid);
            preparedP = source.preparePressure(testCase.PressureGrid);
            preparedX = source.prepareVelocityX(testCase.VelocityXGrid);
            preparedY = source.prepareVelocityY(testCase.VelocityYGrid);

            testCase.verifyTrue(issparse(preparedP.spatialWeights));
            testCase.verifyTrue(issparse(preparedX.spatialWeights));
            testCase.verifyTrue(issparse(preparedY.spatialWeights));
            testCase.verifyEqual(source.numberOfPressureTerms,1);
            testCase.verifyEqual(source.numberOfVelocityXTerms,1);
            testCase.verifyEqual(source.numberOfVelocityYTerms,1);
            testCase.verifyGreaterThan(norm( ...
                source.evaluatePrepared(preparedP,0.2),'fro'),0);
            testCase.verifyGreaterThan(norm( ...
                source.evaluatePrepared(preparedX,0.2),'fro'),0);
            testCase.verifyGreaterThan(norm( ...
                source.evaluatePrepared(preparedY,0.2),'fro'),0);
        end

        function rejectsPressureStencilAtAnyBoundary(testCase)
            term = MultipoleTerm2D([0,0.5],[0,0],@(t) 1);
            mps = MultipoleSource2D(term,ApproximationOrder=4);

            testCase.verifyError(@() mps.discretize( ...
                testCase.PressureGrid,testCase.VelocityXGrid, ...
                testCase.VelocityYGrid), ...
                'MultipoleSource2D:PressureStencilAtBoundary');
        end

        function rejectsStencilOutsideGrid(testCase)
            term = MultipoleTerm2D([0.99,0.5],[2,0],@(t) 1, ...
                TargetField="velocityX");
            mps = MultipoleSource2D(term,ApproximationOrder=4);

            testCase.verifyError(@() mps.discretize( ...
                testCase.PressureGrid,testCase.VelocityXGrid, ...
                testCase.VelocityYGrid), ...
                'MPSappx:StencilOutsideGrid');
        end

        function multipoleSourceRunsInSolver(testCase)
            medium = AcousticMedium2D.constant( ...
                KappaValue=1,BetaValue=1);
            times = 0:0.005:0.02;
            solver = AcousticSolver2DStaggered( ...
                testCase.PressureGrid,times, ...
                testCase.PressureGrid,times(end),medium, ...
                EnforceCFL=false);

            term = MultipoleTerm2D([0.5,0.5],[0,0],@(t) 1);
            mps = MultipoleSource2D(term,ApproximationOrder=2);
            source = mps.discretize( ...
                testCase.PressureGrid,testCase.VelocityXGrid, ...
                testCase.VelocityYGrid);
            zero = @(x,y) zeros(size(x));

            [P,Vx,Vy] = solver.solve(zero,zero,zero,Source=source);
            testCase.verifyGreaterThan(max(abs(P.values),[],'all'),0);
            testCase.verifyGreaterThan(max(abs(Vx.values),[],'all'),0);
            testCase.verifyGreaterThan(max(abs(Vy.values),[],'all'),0);
        end
    end
end
