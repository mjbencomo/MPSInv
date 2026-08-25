classdef TestAcousticMedium2D < matlab.unittest.TestCase
    % Tests for AcousticMedium2D and sampleMedium2D.

    properties
        Domain = [0,2;0,1]
        XPrimal
        YPrimal
        Grid
    end

    methods (TestMethodSetup)
        function createGrid(testCase)
            testCase.XPrimal = Grid1DUnifPrimal([0,2],5);
            testCase.YPrimal = Grid1DUnifPrimal([0,1],3);
            testCase.Grid = GridSpace2D( ...
                testCase.XPrimal,testCase.YPrimal);
        end
    end

    methods (Test)
        function customFunctionsAndWaveSpeed(testCase)
            medium = AcousticMedium2D( ...
                @(x,y) 2+x+y,@(x,y) 4+x-y);
            x = [0,1;0.5,1.5];
            y = [0,0;1,1];
            kappa = 2+x+y;
            beta = 4+x-y;

            testCase.verifyEqual(medium.evaluateKappa(x,y),kappa);
            testCase.verifyEqual(medium.evaluateBeta(x,y),beta);
            testCase.verifyEqual( ...
                medium.evaluateWaveSpeed(x,y),sqrt(kappa.*beta), ...
                AbsTol=1e-14);
        end

        function defaultMetadataUsesAcousticSIUnits(testCase)
            medium = AcousticMedium2D.constant();

            testCase.verifyEqual(medium.kappaLabel,"Bulk modulus");
            testCase.verifyEqual(medium.kappaUnits,"Pa");
            testCase.verifyEqual(medium.betaLabel,"Buoyancy");
            testCase.verifyEqual(medium.betaUnits,"m^3/kg");
            testCase.verifyEqual(medium.waveSpeedLabel,"Wave speed");
            testCase.verifyEqual(medium.waveSpeedUnits,"m/s");
        end

        function constructorAcceptsCustomMetadata(testCase)
            medium = AcousticMedium2D( ...
                @(x,y) 2+0*x+0*y,@(x,y) 3+0*x+0*y, ...
                KappaLabel="Kappa",KappaUnits="GPa", ...
                BetaLabel="Inverse density",BetaUnits="cm^3/g", ...
                WaveSpeedLabel="Sound speed",WaveSpeedUnits="km/s");

            testCase.verifyEqual(medium.kappaLabel,"Kappa");
            testCase.verifyEqual(medium.kappaUnits,"GPa");
            testCase.verifyEqual(medium.betaLabel,"Inverse density");
            testCase.verifyEqual(medium.betaUnits,"cm^3/g");
            testCase.verifyEqual(medium.waveSpeedLabel,"Sound speed");
            testCase.verifyEqual(medium.waveSpeedUnits,"km/s");
        end

        function constantModelPreservesArraySize(testCase)
            medium = AcousticMedium2D.constant( ...
                KappaValue=2.5,BetaValue=4);
            [x,y] = ndgrid(0:2,0:3);

            testCase.verifyEqual( ...
                medium.evaluateKappa(x,y),2.5*ones(3,4));
            testCase.verifyEqual( ...
                medium.evaluateBeta(x,y),4*ones(3,4));
        end

        function linearXAllowsIncreasingAndDecreasingRanges(testCase)
            medium = AcousticMedium2D.linearX( ...
                KappaRange=[2,6],BetaRange=[5,3], ...
                Domain=testCase.Domain);
            x = [0,1,2];
            y = zeros(size(x));

            testCase.verifyEqual( ...
                medium.evaluateKappa(x,y),[2,4,6],AbsTol=1e-14);
            testCase.verifyEqual( ...
                medium.evaluateBeta(x,y),[5,4,3],AbsTol=1e-14);
        end

        function linearYAllowsIncreasingAndDecreasingRanges(testCase)
            medium = AcousticMedium2D.linearY( ...
                KappaRange=[2,6],BetaRange=[5,3], ...
                Domain=testCase.Domain);
            y = [0,0.5,1];
            x = zeros(size(y));

            testCase.verifyEqual( ...
                medium.evaluateKappa(x,y),[2,4,6],AbsTol=1e-14);
            testCase.verifyEqual( ...
                medium.evaluateBeta(x,y),[5,4,3],AbsTol=1e-14);
        end

        function smoothModelHasExpectedExtrema(testCase)
            medium = AcousticMedium2D.smooth( ...
                KappaRange=[2,6],BetaRange=[4,8], ...
                Domain=[0,2;-1,1]);
            x = [0.5,0.5,0];
            y = [-0.5,0.5,-0.5];

            testCase.verifyEqual( ...
                medium.evaluateKappa(x,y),[6,2,4],AbsTol=1e-14);
            testCase.verifyEqual( ...
                medium.evaluateBeta(x,y),[8,4,6],AbsTol=1e-14);
        end

        function gaussianPeakHasExpectedValues(testCase)
            medium = AcousticMedium2D.gaussianPeak( ...
                KappaRange=[2,6],BetaRange=[3,5], ...
                Center=[1,0.5],Width=[0.5,0.25]);
            x = [1,1.5];
            y = [0.5,0.5];

            testCase.verifyEqual( ...
                medium.evaluateKappa(x,y),[6,2+4/exp(1)], ...
                AbsTol=1e-14);
            testCase.verifyEqual( ...
                medium.evaluateBeta(x,y),[5,3+2/exp(1)], ...
                AbsTol=1e-14);
        end

        function gaussianDipHasExpectedValues(testCase)
            medium = AcousticMedium2D.gaussianDip( ...
                KappaRange=[2,6],BetaRange=[3,5], ...
                Center=[1,0.5],Width=[0.5,0.25]);
            x = [1,1.5];
            y = [0.5,0.5];

            testCase.verifyEqual( ...
                medium.evaluateKappa(x,y),[2,6-4/exp(1)], ...
                AbsTol=1e-14);
            testCase.verifyEqual( ...
                medium.evaluateBeta(x,y),[3,5-2/exp(1)], ...
                AbsTol=1e-14);
        end

        function layeredXUsesSpecifiedRelativeWidths(testCase)
            medium = AcousticMedium2D.layeredX( ...
                KappaValues=[2,4,6],BetaValues=[3,5,7], ...
                RelativeWidths=[1,2,1],Domain=testCase.Domain);
            x = [0.25,1,1.75];
            y = zeros(size(x));

            testCase.verifyEqual(medium.evaluateKappa(x,y),[2,4,6]);
            testCase.verifyEqual(medium.evaluateBeta(x,y),[3,5,7]);
        end

        function layeredYUsesSpecifiedRelativeWidths(testCase)
            medium = AcousticMedium2D.layeredY( ...
                KappaValues=[2,4,6],BetaValues=[3,5,7], ...
                RelativeWidths=[1,2,1],Domain=testCase.Domain);
            y = [0.125,0.5,0.875];
            x = zeros(size(y));

            testCase.verifyEqual(medium.evaluateKappa(x,y),[2,4,6]);
            testCase.verifyEqual(medium.evaluateBeta(x,y),[3,5,7]);
        end

        function circularModelDistinguishesInsideAndOutside(testCase)
            medium = AcousticMedium2D.circular( ...
                KappaInside=2,KappaOutside=6, ...
                BetaInside=3,BetaOutside=5, ...
                Center=[1,0.5],Radius=0.25);
            x = [1,1.25,1.5];
            y = [0.5,0.5,0.5];

            testCase.verifyEqual(medium.evaluateKappa(x,y),[2,2,6]);
            testCase.verifyEqual(medium.evaluateBeta(x,y),[3,3,5]);
        end

        function tabulatedModelReproducesBilinearData(testCase)
            xData = [0,1,2];
            yData = [0,0.5,1];
            [X,Y] = ndgrid(xData,yData);
            kappaData = 2+X+2*Y+X.*Y;
            betaData = 4+2*X+Y-X.*Y;
            medium = AcousticMedium2D.tabulated( ...
                XCoordinates=xData,YCoordinates=yData, ...
                KappaValues=kappaData,BetaValues=betaData);
            xq = [0.5,1.5];
            yq = [0.25,0.75];

            testCase.verifyEqual( ...
                medium.evaluateKappa(xq,yq), ...
                2+xq+2*yq+xq.*yq,AbsTol=1e-14);
            testCase.verifyEqual( ...
                medium.evaluateBeta(xq,yq), ...
                4+2*xq+yq-xq.*yq,AbsTol=1e-14);
        end

        function samplingOnStaggeredGrids(testCase)
            medium = AcousticMedium2D( ...
                @(x,y) 2+x+y,@(x,y) 4+x-y);

            [K,Bx,By,C,cMax] = sampleMedium2D(testCase.Grid,medium);
            [X,Y] = testCase.Grid.mesh();
            [Xbx,Ybx] = Bx.grid.mesh();
            [Xby,Yby] = By.grid.mesh();
            expectedC = sqrt((2+X+Y).*(4+X-Y));

            testCase.verifyEqual(K.N,[5,3]);
            testCase.verifyEqual(Bx.N,[4,3]);
            testCase.verifyEqual(By.N,[5,2]);
            testCase.verifyClass(Bx.grid.x,'Grid1DUnifDual');
            testCase.verifyClass(By.grid.y,'Grid1DUnifDual');
            testCase.verifyEqual(K.values,2+X+Y,AbsTol=1e-14);
            testCase.verifyEqual(Bx.values,4+Xbx-Ybx,AbsTol=1e-14);
            testCase.verifyEqual(By.values,4+Xby-Yby,AbsTol=1e-14);
            testCase.verifyEqual(C.values,expectedC,AbsTol=1e-14);
            testCase.verifyEqual(cMax,max(expectedC,[],'all'), ...
                AbsTol=1e-14);
        end

        function samplingBetaOnPrimalGrid(testCase)
            medium = AcousticMedium2D.constant( ...
                KappaValue=2,BetaValue=3);

            [~,Bx,By] = sampleMedium2D( ...
                testCase.Grid,medium,BetaGrid="primal");

            testCase.verifyEqual(Bx.grid.x.pts,testCase.Grid.x.pts);
            testCase.verifyEqual(Bx.grid.y.pts,testCase.Grid.y.pts);
            testCase.verifyEqual(By.grid.x.pts,testCase.Grid.x.pts);
            testCase.verifyEqual(By.grid.y.pts,testCase.Grid.y.pts);
            testCase.verifyEqual(Bx.values,3*ones(5,3));
            testCase.verifyEqual(By.values,3*ones(5,3));
        end

        function samplingPropagatesEditedMetadata(testCase)
            medium = AcousticMedium2D.constant();
            medium.kappaLabel = "Kappa";
            medium.kappaUnits = "kPa";
            medium.betaLabel = "Reciprocal density";
            medium.betaUnits = "custom beta units";
            medium.waveSpeedLabel = "Sound speed";
            medium.waveSpeedUnits = "km/s";

            [K,Bx,By,C] = sampleMedium2D(testCase.Grid,medium);

            testCase.verifyEqual(K.label,"Kappa");
            testCase.verifyEqual(K.units,"kPa");
            testCase.verifyEqual(Bx.label,"Reciprocal density");
            testCase.verifyEqual(Bx.units,"custom beta units");
            testCase.verifyEqual(By.label,"Reciprocal density");
            testCase.verifyEqual(By.units,"custom beta units");
            testCase.verifyEqual(C.label,"Sound speed");
            testCase.verifyEqual(C.units,"km/s");
        end

        function staggeredSamplingRequiresPrimalInput(testCase)
            xDual = Grid1DUnifDual(testCase.XPrimal);
            invalidGrid = GridSpace2D(xDual,testCase.YPrimal);
            medium = AcousticMedium2D.constant();

            testCase.verifyError( ...
                @() sampleMedium2D(invalidGrid,medium), ...
                'sampleMedium2D:PrimalGridRequired');
        end

        function layeredModelRejectsMismatchedData(testCase)
            testCase.verifyError( ...
                @() AcousticMedium2D.layeredX( ...
                    KappaValues=[2,3],BetaValues=[4,5,6], ...
                    RelativeWidths=[1,1]), ...
                'AcousticMedium2D:InvalidLayerData');
        end

        function tabulatedModelRejectsIncorrectArraySize(testCase)
            testCase.verifyError( ...
                @() AcousticMedium2D.tabulated( ...
                    XCoordinates=[0,1,2],YCoordinates=[0,1], ...
                    KappaValues=ones(2,3),BetaValues=ones(3,2)), ...
                'AcousticMedium2D:InvalidTabulatedSize');
        end

        function tabulatedModelRejectsUnorderedCoordinates(testCase)
            testCase.verifyError( ...
                @() AcousticMedium2D.tabulated( ...
                    XCoordinates=[0,2,1],YCoordinates=[0,1], ...
                    KappaValues=ones(3,2),BetaValues=ones(3,2)), ...
                'AcousticMedium2D:InvalidCoordinates');
        end

        function smoothModelRejectsDecreasingRange(testCase)
            testCase.verifyError( ...
                @() AcousticMedium2D.smooth(KappaRange=[2,1]), ...
                'AcousticMedium2D:InvalidRange');
        end

        function modelRejectsInvalidDomain(testCase)
            testCase.verifyError( ...
                @() AcousticMedium2D.linearX( ...
                    Domain=[1,0;0,1]), ...
                'AcousticMedium2D:InvalidDomain');
        end
    end
end
