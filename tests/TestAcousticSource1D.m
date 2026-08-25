classdef TestAcousticSource1D < matlab.unittest.TestCase
    % Tests for separable AcousticSource1D source components and sampling.

    properties
        PrimalGrid
        DualGrid
    end

    methods (TestMethodSetup)
        function createGrids(testCase)
            primal = Grid1DUnifPrimal([-1,1],9);
            testCase.PrimalGrid = GridSpace1D(primal);
            testCase.DualGrid = GridSpace1D(Grid1DUnifDual(primal));
        end
    end

    methods (Test)
        function defaultSourceIsZero(testCase)
            source = AcousticSource1D();
            x = [-1,0,1];

            testCase.verifyEqual(source.evaluatePressure(x,0.25),zeros(size(x)));
            testCase.verifyEqual(source.evaluateVelocity(x,0.25),zeros(size(x)));
            testCase.verifyEqual(source.pressureUnits,"Pa/s");
            testCase.verifyEqual(source.velocityUnits,"m/s^2");
        end

        function zeroFactoryReturnsZeroSource(testCase)
            source = AcousticSource1D.zero();
            x = linspace(-1,1,5).';

            testCase.verifyEqual(source.evaluatePressure(x,2),zeros(5,1));
            testCase.verifyEqual(source.evaluateVelocity(x,2),zeros(5,1));
        end

        function evaluatesSeparablePressureSource(testCase)
            source = AcousticSource1D( ...
                @(x) 1+x.^2,@(t) cos(2*t));
            x = [-1,0,2];
            t = 0.3;

            testCase.verifyEqual( ...
                source.evaluatePressureSpatial(x),1+x.^2);
            testCase.verifyEqual( ...
                source.evaluatePressureTime(t),cos(2*t));
            testCase.verifyEqual( ...
                source.evaluatePressure(x,t),(1+x.^2)*cos(2*t), ...
                AbsTol=1e-14);
        end

        function evaluatesSeparableVelocitySource(testCase)
            source = AcousticSource1D( ...
                @(x) 0,@(t) 1, ...
                @(x) 2-x,@(t) exp(-t));
            x = [-1,0,1];
            t = 0.4;

            testCase.verifyEqual( ...
                source.evaluateVelocitySpatial(x),2-x);
            testCase.verifyEqual( ...
                source.evaluateVelocityTime(t),exp(-t));
            testCase.verifyEqual( ...
                source.evaluateVelocity(x,t),(2-x)*exp(-t), ...
                AbsTol=1e-14);
        end

        function scalarSpatialOutputsAreExpanded(testCase)
            source = AcousticSource1D(@(x) 3,@(t) 2);
            x = reshape(1:6,2,3);

            testCase.verifyEqual( ...
                source.evaluatePressureSpatial(x),3*ones(2,3));
            testCase.verifyEqual( ...
                source.evaluatePressure(x,0),6*ones(2,3));
        end

        function samplesSpatialFactorsOnSpecifiedGrids(testCase)
            source = AcousticSource1D( ...
                @(x) 1+x,@(t) 2, ...
                @(x) 3-x.^2,@(t) 4);

            Gp = source.samplePressureSpatial(testCase.PrimalGrid);
            Gv = source.sampleVelocitySpatial(testCase.DualGrid);

            testCase.verifyClass(Gp,'GridFunction');
            testCase.verifyClass(Gv,'GridFunction');
            testCase.verifyEqual(Gp.grid,testCase.PrimalGrid);
            testCase.verifyEqual(Gv.grid,testCase.DualGrid);
            testCase.verifyEqual( ...
                Gp.values,1+testCase.PrimalGrid.x.pts,AbsTol=1e-14);
            testCase.verifyEqual( ...
                Gv.values,3-testCase.DualGrid.x.pts.^2,AbsTol=1e-14);
        end

        function samplingPropagatesCustomMetadata(testCase)
            source = AcousticSource1D( ...
                @(x) 1+0*x,@(t) 1, ...
                @(x) 2+0*x,@(t) 1, ...
                PressureLabel="Injected pressure", ...
                PressureUnits="kPa/s", ...
                VelocityLabel="Injected acceleration", ...
                VelocityUnits="cm/s^2");

            Gp = source.samplePressureSpatial(testCase.PrimalGrid);
            Gv = source.sampleVelocitySpatial(testCase.DualGrid);

            testCase.verifyEqual(Gp.label,"Injected pressure");
            testCase.verifyEqual(Gp.units,"kPa/s");
            testCase.verifyEqual(Gv.label,"Injected acceleration");
            testCase.verifyEqual(Gv.units,"cm/s^2");
        end

        function rejectsNonnumericSpatialOutput(testCase)
            source = AcousticSource1D(@(x) "invalid",@(t) 1);

            testCase.verifyError( ...
                @() source.evaluatePressureSpatial([0,1]), ...
                'AcousticSource1D:InvalidSpatialOutput');
        end

        function rejectsIncorrectSpatialOutputSize(testCase)
            source = AcousticSource1D(@(x) ones(numel(x)+1,1),@(t) 1);

            testCase.verifyError( ...
                @() source.evaluatePressureSpatial([0;1]), ...
                'AcousticSource1D:InvalidSpatialSize');
        end

        function rejectsNonfiniteSpatialOutput(testCase)
            source = AcousticSource1D(@(x) NaN(size(x)),@(t) 1);

            testCase.verifyError( ...
                @() source.evaluatePressureSpatial([0,1]), ...
                'AcousticSource1D:NonfiniteSpatialOutput');
        end

        function rejectsNonscalarTimeOutput(testCase)
            source = AcousticSource1D(@(x) 1+0*x,@(t) [1,2]);

            testCase.verifyError( ...
                @() source.evaluatePressureTime(0), ...
                'AcousticSource1D:InvalidTimeOutput');
        end

        function evaluatesSumOfSeparableTerms(testCase)
            source = AcousticSource1D.zero();
            source.addPressureTerm(@(x) 1+x,@(t) 2*t);
            source.addPressureTerm(@(x) x.^2,@(t) 1-t);
            x = [-1,0,2];
            t = 0.25;

            expected = (1+x)*(2*t)+x.^2*(1-t);
            testCase.verifyEqual(source.evaluatePressure(x,t),expected, ...
                AbsTol=1e-14);
            testCase.verifyEqual(source.numberOfPressureTerms,2);
            testCase.verifyError(@() source.evaluatePressureTime(t), ...
                'AcousticSource1D:NonseparableComponent');
        end

        function preparedSourceMatchesDirectEvaluation(testCase)
            source = AcousticSource1D.zero();
            source.addPressureTerm(@(x) 1+x,@(t) cos(t));
            source.addPressureTerm(@(x) x.^2,@(t) exp(-t));
            t = 0.3;

            prepared = source.preparePressure(testCase.PrimalGrid);
            actual = source.evaluatePrepared(prepared,t);
            expected = source.evaluatePressure( ...
                testCase.PrimalGrid.x.pts,t);

            testCase.verifyEqual(actual,expected,AbsTol=1e-14);
        end

        function localizedTermsUseSparseSpatialStorage(testCase)
            source = AcousticSource1D.zero();
            source.addLocalizedPressureTerm( ...
                testCase.PrimalGrid,[3,5],[2,-1],@(t) 1+t);

            prepared = source.preparePressure(testCase.PrimalGrid);
            values = source.evaluatePrepared(prepared,2);
            expected = zeros(testCase.PrimalGrid.N,1);
            expected([3,5]) = 3*[2;-1];

            testCase.verifyTrue(issparse(prepared.spatialWeights));
            testCase.verifyEqual(values,expected);
        end

        function localizedTermsRequireTheirConstructionGrid(testCase)
            source = AcousticSource1D.zero();
            source.addLocalizedPressureTerm( ...
                testCase.PrimalGrid,4,1,@(t) 1);
            otherGrid = GridSpace1D(Grid1DUnifPrimal([-1,1],11));

            testCase.verifyError(@() source.preparePressure(otherGrid), ...
                'AcousticSource1D:LocalizedGridMismatch');
        end

        function combinesIndependentSources(testCase)
            source1 = AcousticSource1D(@(x) x,@(t) 2);
            source2 = AcousticSource1D(@(x) 1+x.^2,@(t) 3);
            combined = AcousticSource1D.combine({source1,source2});
            x = testCase.PrimalGrid.x.pts;

            testCase.verifyEqual(combined.evaluatePressure(x,0), ...
                source1.evaluatePressure(x,0)+ ...
                source2.evaluatePressure(x,0),AbsTol=1e-14);
        end
    end
end
