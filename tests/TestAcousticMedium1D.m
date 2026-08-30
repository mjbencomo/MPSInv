classdef TestAcousticMedium1D < matlab.unittest.TestCase
    % Tests for AcousticMedium1D and sampleMedium1D.

    properties
        Domain = [-1,1]
        NumPoints = 9
        PrimalGrid
    end

    methods (TestMethodSetup)
        function createPrimalGrid(testCase)
            x = Grid1DUnifPrimal(testCase.Domain,testCase.NumPoints);
            testCase.PrimalGrid = GridSpace1D(x);
        end
    end

    methods (Test)
        function customFunctionsEvaluateAtCoordinates(testCase)
            medium = AcousticMedium1D(@(x) 2+x.^2,@(x) 3-x);
            x = [-1,0,2];

            testCase.verifyEqual(medium.evaluateKappa(x),[3,2,6]);
            testCase.verifyEqual(medium.evaluateBeta(x),[4,3,1]);
        end

        function defaultMetadataUsesAcousticSIUnits(testCase)
            medium = AcousticMedium1D.constant();

            testCase.verifyEqual(medium.kappaLabel,"Bulk modulus");
            testCase.verifyEqual(medium.kappaUnits,"Pa");
            testCase.verifyEqual(medium.betaLabel,"Buoyancy");
            testCase.verifyEqual(medium.betaUnits,"m^3/kg");
            testCase.verifyEqual(medium.waveSpeedLabel,"Wave speed");
            testCase.verifyEqual(medium.waveSpeedUnits,"m/s");
        end

        function constructorAcceptsCustomMetadata(testCase)
            medium = AcousticMedium1D(@(x) 2+0*x,@(x) 3+0*x, ...
                KappaLabel="Compressibility parameter", ...
                KappaUnits="GPa", ...
                BetaLabel="Inverse density", ...
                BetaUnits="cm^3/g", ...
                WaveSpeedLabel="Sound speed", ...
                WaveSpeedUnits="km/s");

            testCase.verifyEqual( ...
                medium.kappaLabel,"Compressibility parameter");
            testCase.verifyEqual(medium.kappaUnits,"GPa");
            testCase.verifyEqual(medium.betaLabel,"Inverse density");
            testCase.verifyEqual(medium.betaUnits,"cm^3/g");
            testCase.verifyEqual(medium.waveSpeedLabel,"Sound speed");
            testCase.verifyEqual(medium.waveSpeedUnits,"km/s");
        end

        function constantModelReturnsConstantArrays(testCase)
            medium = AcousticMedium1D.constant( ...
                KappaValue=2.5,BetaValue=4.0);
            x = reshape(linspace(-1,1,6),2,3);

            testCase.verifyEqual(medium.evaluateKappa(x),2.5*ones(2,3));
            testCase.verifyEqual(medium.evaluateBeta(x),4.0*ones(2,3));
            testCase.verifyEqual( ...
                medium.evaluateWaveSpeed(x),sqrt(10)*ones(2,3), ...
                AbsTol=1e-14);
        end

        function linearModelHasSpecifiedEndpointValues(testCase)
            medium = AcousticMedium1D.linear( ...
                KappaRange=[2,6], ...
                BetaRange=[3,5], ...
                Domain=testCase.Domain);
            x = [-1,0,1];

            testCase.verifyEqual(medium.evaluateKappa(x),[2,4,6], ...
                AbsTol=1e-14);
            testCase.verifyEqual(medium.evaluateBeta(x),[3,4,5], ...
                AbsTol=1e-14);
        end

        function smoothModelHasExpectedQuarterPeriodValues(testCase)
            medium = AcousticMedium1D.smooth( ...
                KappaRange=[2,6], ...
                BetaRange=[4,8], ...
                Domain=[-1,3]);
            x = [-1,0,1,2,3];

            testCase.verifyEqual(medium.evaluateKappa(x),[4,6,4,2,4], ...
                AbsTol=1e-14);
            testCase.verifyEqual(medium.evaluateBeta(x),[6,8,6,4,6], ...
                AbsTol=1e-14);
        end

        function gaussianPeakHasCorrectCenterAndWidthValues(testCase)
            medium = AcousticMedium1D.gaussian( ...
                KappaInside=6, KappaOutside=2, ...
                BetaInside=5, BetaOutside=3, ...
                Center=0.25, Width=0.5);
            x = [0.25,0.75];

            testCase.verifyEqual(medium.evaluateKappa(x),[6,2+4/exp(1)], ...
                AbsTol=1e-14);
            testCase.verifyEqual(medium.evaluateBeta(x),[5,3+2/exp(1)], ...
                AbsTol=1e-14);
        end

        function gaussianDipHasCorrectCenterAndWidthValues(testCase)
            medium = AcousticMedium1D.gaussian( ...
                KappaInside=2, KappaOutside=6, ...
                BetaInside=3, BetaOutside=5, ...
                Center=0.25, Width=0.5);
            x = [0.25,0.75];

            testCase.verifyEqual(medium.evaluateKappa(x),[2,6-4/exp(1)], ...
                AbsTol=1e-14);
            testCase.verifyEqual(medium.evaluateBeta(x),[3,5-2/exp(1)], ...
                AbsTol=1e-14);
        end

        function layeredUsesSpecifiedRelativeWidths(testCase)
            medium = AcousticMedium1D.layered( ...
                KappaValues=[2,4,6], BetaValues=[3,5,7], ...
                RelativeWidths=[1,2,1], Domain=testCase.Domain);
            % domain = [-1,1] -> layer edges at [-1, -0.5, 0.5, 1]
            x = [-0.75, 0, 0.75];

            testCase.verifyEqual(medium.evaluateKappa(x), [2,4,6]);
            testCase.verifyEqual(medium.evaluateBeta(x), [3,5,7]);
        end
        
        function samplingOnPrimalGridUsesSameGrid(testCase)
            medium = AcousticMedium1D(@(x) 2+x,@(x) 4-x.^2);

            [K,B] = sampleMedium1D( ...
                testCase.PrimalGrid,medium,BetaGrid="primal");
            x = testCase.PrimalGrid.x.pts;

            testCase.verifyClass(K,'GridFunction');
            testCase.verifyClass(B,'GridFunction');
            testCase.verifyEqual(K.grid,testCase.PrimalGrid);
            testCase.verifyEqual(B.grid,testCase.PrimalGrid);
            testCase.verifyEqual(K.values,2+x);
            testCase.verifyEqual(B.values,4-x.^2);
        end

        function samplingOnDualGridUsesMidpoints(testCase)
            medium = AcousticMedium1D(@(x) 2+x,@(x) 4-x.^2);

            [K,B] = sampleMedium1D( ...
                testCase.PrimalGrid,medium,BetaGrid="dual");

            primalPoints = testCase.PrimalGrid.x.pts;
            dualPoints = 0.5*(primalPoints(1:end-1)+primalPoints(2:end));

            testCase.verifyClass(B.grid.x,'Grid1DUnifDual');
            testCase.verifyEqual(K.N,testCase.NumPoints);
            testCase.verifyEqual(B.N,testCase.NumPoints-1);
            testCase.verifyEqual(B.grid.x.pts,dualPoints,AbsTol=1e-14);
            testCase.verifyEqual(B.values,4-dualPoints.^2,AbsTol=1e-14);
        end

        function defaultBetaGridIsDual(testCase)
            medium = AcousticMedium1D.constant( ...
                KappaValue=2,BetaValue=3);

            [~,B] = sampleMedium1D(testCase.PrimalGrid,medium);

            testCase.verifyClass(B.grid.x,'Grid1DUnifDual');
            testCase.verifyEqual(B.N,testCase.NumPoints-1);
        end

        function samplingReturnsWaveSpeedAndMaximum(testCase)
            medium = AcousticMedium1D(@(x) 2+x,@(x) 4-x);

            [~,~,C,cMax] = sampleMedium1D( ...
                testCase.PrimalGrid,medium);
            x = testCase.PrimalGrid.x.pts;
            expected = sqrt((2+x).*(4-x));

            testCase.verifyClass(C,'GridFunction');
            testCase.verifyEqual(C.grid,testCase.PrimalGrid);
            testCase.verifyEqual(C.values,expected,AbsTol=1e-14);
            testCase.verifyEqual(cMax,max(expected),AbsTol=1e-14);
        end

        function samplingPropagatesEditedMetadata(testCase)
            medium = AcousticMedium1D.constant();
            medium.kappaLabel = "Kappa";
            medium.kappaUnits = "kPa";
            medium.betaLabel = "Reciprocal density";
            medium.betaUnits = "custom beta units";
            medium.waveSpeedLabel = "Sound speed";
            medium.waveSpeedUnits = "km/s";

            [K,B,C] = sampleMedium1D(testCase.PrimalGrid,medium);

            testCase.verifyEqual(K.label,"Kappa");
            testCase.verifyEqual(K.units,"kPa");
            testCase.verifyEqual(B.label,"Reciprocal density");
            testCase.verifyEqual(B.units,"custom beta units");
            testCase.verifyEqual(C.label,"Sound speed");
            testCase.verifyEqual(C.units,"km/s");
        end

        function dualSamplingRequiresPrimalInput(testCase)
            primal = Grid1DUnifPrimal(testCase.Domain,testCase.NumPoints);
            dualGrid = GridSpace1D(Grid1DUnifDual(primal));
            medium = AcousticMedium1D.constant();

            testCase.verifyError( ...
                @() sampleMedium1D(dualGrid,medium,BetaGrid="dual"), ...
                'sampleMedium1D:PrimalGridRequired');
        end

        function smoothRejectsDecreasingRange(testCase)
            testCase.verifyError( ...
                @() AcousticMedium1D.smooth(KappaRange=[2,1]), ...
                'AcousticMedium1D:InvalidRange');
        end

        function linearAllowsDecreasingRange(testCase)
            medium = AcousticMedium1D.linear( ...
                KappaRange=[6,2], ...
                BetaRange=[5,3], ...
                Domain=testCase.Domain);
            x = [-1,0,1];

            testCase.verifyEqual(medium.evaluateKappa(x),[6,4,2], ...
                AbsTol=1e-14);
            testCase.verifyEqual(medium.evaluateBeta(x),[5,4,3], ...
                AbsTol=1e-14);
        end
    end
end
