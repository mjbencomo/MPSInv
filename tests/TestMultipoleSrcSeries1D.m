classdef TestMultipoleSrcSeries1D < matlab.unittest.TestCase
    properties
        PressureGrid
        VelocityGrid
    end

    methods (TestMethodSetup)
        function createGrids(testCase)
            primal = Grid1DUnifPrimal([0,1],41);
            testCase.PressureGrid = GridSpace1D(primal);
            testCase.VelocityGrid = ...
                GridSpace1D(Grid1DUnifDual(primal));
        end
    end

    methods (Test)
        function termStoresComponentAndTimeFunction(testCase)
            component = MultipoleComponent1D( ...
                0.4,1,TargetField="velocity");
            term = MultipoleSrcTerm1D(component,@(t) -2*(1+t));

            testCase.verifyEqual(term.component,component);
            testCase.verifyEqual(term.location,0.4);
            testCase.verifyEqual(term.derivativeOrder,1);
            testCase.verifyEqual(term.targetField,"velocity");
            testCase.verifyEqual(term.evaluateTime(0.5),-3);
            testCase.verifyFalse(isprop(term,'amplitude'));
            testCase.verifyTrue(isa(term,'MultipoleSrcTerm'));
        end

        function rejectsInvalidTimeOutput(testCase)
            component = MultipoleComponent1D(0.5,0);
            term = MultipoleSrcTerm1D(component,@(t) [1,2]);

            testCase.verifyError(@() term.evaluateTime(0), ...
                'MultipoleSrcTerm:InvalidTimeOutput');
        end

        function singleTermBuildsAcousticSource(testCase)
            component = MultipoleComponent1D(0.45,1);
            term = MultipoleSrcTerm1D(component,@(t) 1+t);

            source = term.discretize(testCase.PressureGrid, ...
                testCase.VelocityGrid,4);

            testCase.verifyTrue(isa(source,'AcousticSource'));
            testCase.verifyTrue(isa(source,'AcousticSource1D'));
            testCase.verifyEqual(source.numberOfPressureTerms,1);
        end

        function seriesExposesTermsComponentsAndLocations(testCase)
            c1 = MultipoleComponent1D(0.35,0);
            c2 = MultipoleComponent1D( ...
                0.65,1,TargetField="velocity");
            terms = [ ...
                MultipoleSrcTerm1D(c1,@(t) 1), ...
                MultipoleSrcTerm1D(c2,@(t) 2-t)];

            series = MultipoleSrcSeries1D(terms);

            testCase.verifyEqual(series.numberOfTerms,2);
            testCase.verifyEqual(series.terms,terms);
            testCase.verifyEqual(series.components,[c1,c2]);
            testCase.verifyEqual(series.locations,[0.35;0.65]);
            testCase.verifyFalse( ...
                isprop(series,'approximationOrders'));
        end

        function discretizesPressureAndVelocityTerms(testCase)
            pressureComponent = MultipoleComponent1D(0.35,0);
            velocityComponent = MultipoleComponent1D( ...
                0.65,1,TargetField="velocity");
            series = MultipoleSrcSeries1D([ ...
                MultipoleSrcTerm1D(pressureComponent,@(t) 1+t), ...
                MultipoleSrcTerm1D(velocityComponent,@(t) 2-t)]);

            source = series.discretize(testCase.PressureGrid, ...
                testCase.VelocityGrid,4);
            preparedP = source.preparePressure(testCase.PressureGrid);
            preparedV = source.prepareVelocity(testCase.VelocityGrid);

            testCase.verifyEqual(source.numberOfPressureTerms,1);
            testCase.verifyEqual(source.numberOfVelocityTerms,1);
            testCase.verifyGreaterThan( ...
                norm(source.evaluatePrepared(preparedP,0.2)),0);
            testCase.verifyGreaterThan( ...
                norm(source.evaluatePrepared(preparedV,0.2)),0);
        end

        function acceptsPerTermApproximationOrders(testCase)
            c1 = MultipoleComponent1D(0.357,0);
            c2 = MultipoleComponent1D(0.643,1);
            series = MultipoleSrcSeries1D([ ...
                MultipoleSrcTerm1D(c1,@(t) 1), ...
                MultipoleSrcTerm1D(c2,@(t) 1)]);

            source = series.discretize(testCase.PressureGrid, ...
                testCase.VelocityGrid,[2,4]);
            prepared = source.preparePressure(testCase.PressureGrid);

            testCase.verifyEqual(nnz(prepared.spatialWeights(:,1)),2);
            testCase.verifyEqual(nnz(prepared.spatialWeights(:,2)),5);
        end

        function rejectsWrongNumberOfApproximationOrders(testCase)
            component = MultipoleComponent1D(0.5,0);
            term = MultipoleSrcTerm1D(component,@(t) 1);
            series = MultipoleSrcSeries1D([term,term]);

            testCase.verifyError(@() series.discretize( ...
                testCase.PressureGrid,testCase.VelocityGrid,[2,4,6]), ...
                'MultipoleSrcSeries:ApproximationOrderSizeMismatch');
        end

        function rejectsPressureStencilAtBoundary(testCase)
            component = MultipoleComponent1D(0,0);
            term = MultipoleSrcTerm1D(component,@(t) 1);

            testCase.verifyError(@() term.discretize( ...
                testCase.PressureGrid,testCase.VelocityGrid,4), ...
                'MultipoleSrcTerm1D:PressureStencilAtBoundary');
        end

        function emptySeriesProducesZeroSource(testCase)
            series = MultipoleSrcSeries1D( ...
                MultipoleSrcTerm1D.empty(1,0));
            source = series.discretize(testCase.PressureGrid, ...
                testCase.VelocityGrid,4);

            testCase.verifyEqual(series.numberOfTerms,0);
            testCase.verifySize(series.locations,[0,1]);
            testCase.verifyEqual(source.numberOfPressureTerms,0);
            testCase.verifyEqual(source.numberOfVelocityTerms,0);
        end
    end
end
