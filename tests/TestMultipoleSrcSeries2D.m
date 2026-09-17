classdef TestMultipoleSrcSeries2D < matlab.unittest.TestCase
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
        function termStoresComponentAndTimeFunction(testCase)
            component = MultipoleComponent2D( ...
                [0.4,0.6],[1,2],TargetField="velocityY");
            term = MultipoleSrcTerm2D(component,@(t) -2*(1+t));

            testCase.verifyEqual(term.component,component);
            testCase.verifyEqual(term.location,[0.4,0.6]);
            testCase.verifyEqual(term.derivativeOrder,[1,2]);
            testCase.verifyEqual(term.targetField,"velocityY");
            testCase.verifyEqual(term.evaluateTime(0.5),-3);
            testCase.verifyFalse(isprop(term,'amplitude'));
            testCase.verifyTrue(isa(term,'MultipoleSrcTerm'));
        end

        function rejectsInvalidTimeOutput(testCase)
            component = MultipoleComponent2D([0.5,0.5],[0,0]);
            term = MultipoleSrcTerm2D(component,@(t) Inf);

            testCase.verifyError(@() term.evaluateTime(0), ...
                'MultipoleSrcTerm:InvalidTimeOutput');
        end

        function singleTermBuildsAcousticSource(testCase)
            component = MultipoleComponent2D( ...
                [0.45,0.55],[1,0]);
            term = MultipoleSrcTerm2D(component,@(t) 1+t);

            source = term.discretize(testCase.PressureGrid, ...
                testCase.VelocityXGrid,testCase.VelocityYGrid,4);

            testCase.verifyTrue(isa(source,'AcousticSource'));
            testCase.verifyTrue(isa(source,'AcousticSource2D'));
            testCase.verifyEqual(source.numberOfPressureTerms,1);
        end

        function seriesExposesTermsComponentsAndLocations(testCase)
            c1 = MultipoleComponent2D([0.35,0.45],[0,0]);
            c2 = MultipoleComponent2D( ...
                [0.65,0.55],[1,0],TargetField="velocityX");
            terms = [ ...
                MultipoleSrcTerm2D(c1,@(t) 1), ...
                MultipoleSrcTerm2D(c2,@(t) 2-t)];

            series = MultipoleSrcSeries2D(terms);

            testCase.verifyEqual(series.numberOfTerms,2);
            testCase.verifyEqual(series.terms,terms);
            testCase.verifyEqual(series.components,[c1,c2]);
            testCase.verifyEqual(series.locations, ...
                [0.35,0.45;0.65,0.55]);
            testCase.verifyFalse( ...
                isprop(series,'approximationOrders'));
        end

        function discretizesAllTargetFields(testCase)
            pressureComponent = MultipoleComponent2D( ...
                [0.35,0.45],[0,0]);
            velocityXComponent = MultipoleComponent2D( ...
                [0.55,0.45],[1,0],TargetField="velocityX");
            velocityYComponent = MultipoleComponent2D( ...
                [0.45,0.55],[0,1],TargetField="velocityY");
            series = MultipoleSrcSeries2D([ ...
                MultipoleSrcTerm2D(pressureComponent,@(t) 1+t), ...
                MultipoleSrcTerm2D(velocityXComponent,@(t) 2-t), ...
                MultipoleSrcTerm2D(velocityYComponent,@(t) cos(t))]);

            source = series.discretize(testCase.PressureGrid, ...
                testCase.VelocityXGrid,testCase.VelocityYGrid,4);
            preparedP = source.preparePressure(testCase.PressureGrid);
            preparedX = source.prepareVelocityX(testCase.VelocityXGrid);
            preparedY = source.prepareVelocityY(testCase.VelocityYGrid);

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

        function acceptsPerTermApproximationOrders(testCase)
            c1 = MultipoleComponent2D([0.357,0.443],[0,0]);
            c2 = MultipoleComponent2D([0.643,0.557],[1,1]);
            series = MultipoleSrcSeries2D([ ...
                MultipoleSrcTerm2D(c1,@(t) 1), ...
                MultipoleSrcTerm2D(c2,@(t) 1)]);

            source = series.discretize(testCase.PressureGrid, ...
                testCase.VelocityXGrid,testCase.VelocityYGrid,[2,4]);
            prepared = source.preparePressure(testCase.PressureGrid);

            testCase.verifyEqual(nnz(prepared.spatialWeights(:,1)),4);
            testCase.verifyEqual(nnz(prepared.spatialWeights(:,2)),25);
        end

        function rejectsWrongNumberOfApproximationOrders(testCase)
            component = MultipoleComponent2D([0.5,0.5],[0,0]);
            term = MultipoleSrcTerm2D(component,@(t) 1);
            series = MultipoleSrcSeries2D([term,term]);

            testCase.verifyError(@() series.discretize( ...
                testCase.PressureGrid,testCase.VelocityXGrid, ...
                testCase.VelocityYGrid,[2,4,6]), ...
                'MultipoleSrcSeries:ApproximationOrderSizeMismatch');
        end

        function rejectsPressureStencilAtAnyBoundary(testCase)
            component = MultipoleComponent2D([0,0.5],[0,0]);
            term = MultipoleSrcTerm2D(component,@(t) 1);

            testCase.verifyError(@() term.discretize( ...
                testCase.PressureGrid,testCase.VelocityXGrid, ...
                testCase.VelocityYGrid,4), ...
                'MultipoleSrcTerm2D:PressureStencilAtBoundary');
        end

        function emptySeriesProducesZeroSource(testCase)
            series = MultipoleSrcSeries2D( ...
                MultipoleSrcTerm2D.empty(1,0));
            source = series.discretize(testCase.PressureGrid, ...
                testCase.VelocityXGrid,testCase.VelocityYGrid,4);

            testCase.verifyEqual(series.numberOfTerms,0);
            testCase.verifySize(series.locations,[0,2]);
            testCase.verifyEqual(source.numberOfPressureTerms,0);
            testCase.verifyEqual(source.numberOfVelocityXTerms,0);
            testCase.verifyEqual(source.numberOfVelocityYTerms,0);
        end
    end
end
