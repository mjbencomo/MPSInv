classdef TestAcousticSource2D < matlab.unittest.TestCase
    % Tests for separable AcousticSource2D components and grid sampling.

    properties
        PressureGrid
        VelocityXGrid
        VelocityYGrid
    end

    methods (TestMethodSetup)
        function createGrids(testCase)
            xPrimal = Grid1DUnifPrimal([-1,1],7);
            yPrimal = Grid1DUnifPrimal([-2,2],9);
            xDual = Grid1DUnifDual(xPrimal);
            yDual = Grid1DUnifDual(yPrimal);

            testCase.PressureGrid = GridSpace2D(xPrimal,yPrimal);
            testCase.VelocityXGrid = GridSpace2D(xDual,yPrimal);
            testCase.VelocityYGrid = GridSpace2D(xPrimal,yDual);
        end
    end

    methods (Test)
        function defaultSourceIsZero(testCase)
            source = AcousticSource2D();
            [X,Y] = testCase.PressureGrid.mesh();

            testCase.verifyEqual(source.evaluatePressure(X,Y,0.25),zeros(size(X)));
            testCase.verifyEqual(source.evaluateVelocityX(X,Y,0.25),zeros(size(X)));
            testCase.verifyEqual(source.evaluateVelocityY(X,Y,0.25),zeros(size(X)));
            testCase.verifyEqual(source.pressureUnits,"Pa/s");
            testCase.verifyEqual(source.velocityXUnits,"m/s^2");
            testCase.verifyEqual(source.velocityYUnits,"m/s^2");
        end

        function zeroFactoryReturnsZeroSource(testCase)
            source = AcousticSource2D.zero();
            [X,Y] = testCase.PressureGrid.mesh();
            testCase.verifyEqual(source.evaluatePressure(X,Y,2),zeros(size(X)));
        end

        function evaluatesAllSeparableComponents(testCase)
            source = AcousticSource2D( ...
                @(x,y) 1+x.^2+y,@(t) cos(2*t), ...
                @(x,y) x-y,@(t) exp(-t), ...
                @(x,y) x.*y,@(t) 1+t);
            [X,Y] = ndgrid([-1,0,2],[-2,1]);
            t = 0.3;

            testCase.verifyEqual(source.evaluatePressureSpatial(X,Y), ...
                1+X.^2+Y);
            testCase.verifyEqual(source.evaluatePressureTime(t),cos(2*t));
            testCase.verifyEqual(source.evaluatePressure(X,Y,t), ...
                (1+X.^2+Y)*cos(2*t),AbsTol=1e-14);
            testCase.verifyEqual(source.evaluateVelocityX(X,Y,t), ...
                (X-Y)*exp(-t),AbsTol=1e-14);
            testCase.verifyEqual(source.evaluateVelocityY(X,Y,t), ...
                X.*Y*(1+t),AbsTol=1e-14);
        end

        function scalarSpatialOutputsAreExpanded(testCase)
            source = AcousticSource2D(@(x,y) 3,@(t) 2);
            [X,Y] = ndgrid(1:2,1:3);

            testCase.verifyEqual(source.evaluatePressureSpatial(X,Y),3*ones(2,3));
            testCase.verifyEqual(source.evaluatePressure(X,Y,0),6*ones(2,3));
        end

        function samplesComponentsOnStaggeredGrids(testCase)
            source = AcousticSource2D( ...
                @(x,y) 1+x+y,@(t) 2, ...
                @(x,y) 3-x.^2+y,@(t) 4, ...
                @(x,y) 5+x-y.^2,@(t) 6);

            Gp = source.samplePressureSpatial(testCase.PressureGrid);
            Gx = source.sampleVelocityXSpatial(testCase.VelocityXGrid);
            Gy = source.sampleVelocityYSpatial(testCase.VelocityYGrid);
            [Xp,Yp] = testCase.PressureGrid.mesh();
            [Xx,Yx] = testCase.VelocityXGrid.mesh();
            [Xy,Yy] = testCase.VelocityYGrid.mesh();

            testCase.verifyClass(Gp,'GridFunction');
            testCase.verifyEqual(Gp.values,1+Xp+Yp,AbsTol=1e-14);
            testCase.verifyEqual(Gx.values,3-Xx.^2+Yx,AbsTol=1e-14);
            testCase.verifyEqual(Gy.values,5+Xy-Yy.^2,AbsTol=1e-14);
        end

        function samplingPropagatesCustomMetadata(testCase)
            source = AcousticSource2D( ...
                @(x,y) ones(size(x)),@(t) 1, ...
                @(x,y) 2*ones(size(x)),@(t) 1, ...
                @(x,y) 3*ones(size(x)),@(t) 1, ...
                PressureLabel="Injected pressure",PressureUnits="kPa/s", ...
                VelocityXLabel="Horizontal acceleration", ...
                VelocityXUnits="cm/s^2", ...
                VelocityYLabel="Vertical acceleration", ...
                VelocityYUnits="mm/s^2");

            Gp = source.samplePressureSpatial(testCase.PressureGrid);
            Gx = source.sampleVelocityXSpatial(testCase.VelocityXGrid);
            Gy = source.sampleVelocityYSpatial(testCase.VelocityYGrid);

            testCase.verifyEqual(Gp.label,"Injected pressure");
            testCase.verifyEqual(Gp.units,"kPa/s");
            testCase.verifyEqual(Gx.label,"Horizontal acceleration");
            testCase.verifyEqual(Gx.units,"cm/s^2");
            testCase.verifyEqual(Gy.label,"Vertical acceleration");
            testCase.verifyEqual(Gy.units,"mm/s^2");
        end

        function rejectsInvalidCoordinates(testCase)
            source = AcousticSource2D(@(x,y) x+y,@(t) 1);
            testCase.verifyError( ...
                @() source.evaluatePressureSpatial(zeros(2),zeros(3)), ...
                'AcousticSource2D:InvalidCoordinates');
        end

        function rejectsNonnumericSpatialOutput(testCase)
            source = AcousticSource2D(@(x,y) "invalid",@(t) 1);
            testCase.verifyError( ...
                @() source.evaluatePressureSpatial(zeros(2),zeros(2)), ...
                'AcousticSource2D:InvalidSpatialOutput');
        end

        function rejectsIncorrectSpatialOutputSize(testCase)
            source = AcousticSource2D(@(x,y) zeros(3,3),@(t) 1);
            testCase.verifyError( ...
                @() source.evaluatePressureSpatial(zeros(2),zeros(2)), ...
                'AcousticSource2D:InvalidSpatialSize');
        end

        function rejectsNonfiniteSpatialOutput(testCase)
            source = AcousticSource2D(@(x,y) NaN(size(x)),@(t) 1);
            testCase.verifyError( ...
                @() source.evaluatePressureSpatial(zeros(2),zeros(2)), ...
                'AcousticSource2D:NonfiniteSpatialOutput');
        end

        function rejectsNonscalarTimeOutput(testCase)
            source = AcousticSource2D(@(x,y) ones(size(x)),@(t) [1,2]);
            testCase.verifyError(@() source.evaluatePressureTime(0), ...
                'AcousticSource2D:InvalidTimeOutput');
        end

        function evaluatesSumOfSeparableTerms(testCase)
            source = AcousticSource2D.zero();
            source.addPressureTerm(@(x,y) 1+x+y,@(t) 2*t);
            source.addPressureTerm(@(x,y) x.*y,@(t) 1-t);
            [X,Y] = ndgrid([-1,0,2],[-2,1]);
            t = 0.25;

            expected = (1+X+Y)*(2*t)+X.*Y*(1-t);
            testCase.verifyEqual(source.evaluatePressure(X,Y,t),expected, ...
                AbsTol=1e-14);
            testCase.verifyEqual(source.numberOfPressureTerms,2);
            testCase.verifyError(@() source.evaluatePressureTime(t), ...
                'AcousticSource2D:NonseparableComponent');
        end

        function preparedSourceMatchesDirectEvaluation(testCase)
            source = AcousticSource2D.zero();
            source.addPressureTerm(@(x,y) 1+x+y,@(t) cos(t));
            source.addPressureTerm(@(x,y) x.*y,@(t) exp(-t));
            [X,Y] = testCase.PressureGrid.mesh();
            t = 0.3;

            actual = source.evaluatePrepared( ...
                source.preparePressure(testCase.PressureGrid),t);
            expected = source.evaluatePressure(X,Y,t);
            testCase.verifyEqual(actual,expected,AbsTol=1e-14);
            testCase.verifySize(actual,testCase.PressureGrid.N);
        end

        function localizedTermsUseSparseSpatialStorage(testCase)
            source = AcousticSource2D.zero();
            indices = [3,12];
            source.addLocalizedPressureTerm( ...
                testCase.PressureGrid,indices,[2,-1],@(t) 1+t);

            prepared = source.preparePressure(testCase.PressureGrid);
            values = source.evaluatePrepared(prepared,2);
            expected = zeros(testCase.PressureGrid.N);
            expected(indices) = 3*[2,-1];

            testCase.verifyTrue(issparse(prepared.spatialWeights));
            testCase.verifyEqual(values,expected);
        end

        function localizedTermsSupportBothVelocityComponents(testCase)
            source = AcousticSource2D.zero();
            source.addLocalizedVelocityXTerm( ...
                testCase.VelocityXGrid,4,2,@(t) t);
            source.addLocalizedVelocityYTerm( ...
                testCase.VelocityYGrid,5,-3,@(t) 1+t);

            vx = source.evaluatePrepared( ...
                source.prepareVelocityX(testCase.VelocityXGrid),2);
            vy = source.evaluatePrepared( ...
                source.prepareVelocityY(testCase.VelocityYGrid),2);

            testCase.verifyEqual(vx(4),4);
            testCase.verifyEqual(vy(5),-9);
        end

        function localizedTermsRequireConstructionGrid(testCase)
            source = AcousticSource2D.zero();
            source.addLocalizedPressureTerm(testCase.PressureGrid,4,1,@(t) 1);
            otherGrid = GridSpace2D( ...
                Grid1DUnifPrimal([-1,1],9),Grid1DUnifPrimal([-2,2],9));

            testCase.verifyError(@() source.preparePressure(otherGrid), ...
                'AcousticSource2D:LocalizedGridMismatch');
        end

        function rejectsLocalizedIndexOutsideGrid(testCase)
            source = AcousticSource2D.zero();
            badIndex = prod(testCase.PressureGrid.N)+1;
            testCase.verifyError(@() source.addLocalizedPressureTerm( ...
                testCase.PressureGrid,badIndex,1,@(t) 1), ...
                'AcousticSource2D:LocalizedIndexOutsideGrid');
        end

        function combinesIndependentSources(testCase)
            source1 = AcousticSource2D(@(x,y) x+y,@(t) 2);
            source2 = AcousticSource2D(@(x,y) 1+x.^2+y.^2,@(t) 3);
            combined = AcousticSource2D.combine({source1,source2});
            [X,Y] = testCase.PressureGrid.mesh();

            testCase.verifyEqual(combined.evaluatePressure(X,Y,0), ...
                source1.evaluatePressure(X,Y,0)+ ...
                source2.evaluatePressure(X,Y,0),AbsTol=1e-14);
        end
    end
end
