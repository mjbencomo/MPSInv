classdef TestGridFunction < matlab.unittest.TestCase
    % Tests for GridFunction construction, assignment, and interpolation.

    properties
        XPrimal
        XDual
        YShort
        YPrimal
        Grid1DPrimal
        Grid1DDual
        Grid2DShort
        Grid2DPrimal
        Grid2DDualPrimal
    end

    methods (TestMethodSetup)
        function createGridFixtures(testCase)
            testCase.XPrimal = Grid1DUnifPrimal([0,1],5);
            testCase.XDual = Grid1DUnifDual(testCase.XPrimal);
            testCase.YShort = Grid1DUnifPrimal([0,2],3);
            testCase.YPrimal = Grid1DUnifPrimal([0,2],5);

            testCase.Grid1DPrimal = GridSpace1D(testCase.XPrimal);
            testCase.Grid1DDual = GridSpace1D(testCase.XDual);
            testCase.Grid2DShort = GridSpace2D( ...
                testCase.XPrimal,testCase.YShort);
            testCase.Grid2DPrimal = GridSpace2D( ...
                testCase.XPrimal,testCase.YPrimal);
            testCase.Grid2DDualPrimal = GridSpace2D( ...
                testCase.XDual,testCase.YPrimal);
        end
    end

    methods (Test)
        function zeroInitialization1D(testCase)
            u = GridFunction(testCase.Grid1DPrimal);

            testCase.verifyEqual(u.values,zeros(5,1));
            testCase.verifyEqual(u.N,5);
            testCase.verifyEqual(u.dim,1);
            testCase.verifyEqual(u.label,"Value");
            testCase.verifyEqual(u.units,"");
            testCase.verifyEqual(u.formattedLabel(),"Value");
        end

        function constructorSetsLabelAndUnits(testCase)
            u = GridFunction(testCase.Grid1DPrimal, ...
                Label="Pressure",Units="Pa");

            testCase.verifyEqual(u.label,"Pressure");
            testCase.verifyEqual(u.units,"Pa");
            testCase.verifyEqual(u.formattedLabel(),"Pressure (Pa)");
        end

        function metadataCanBeEdited(testCase)
            u = GridFunction(testCase.Grid1DPrimal);

            u.label = "Particle velocity";
            u.units = "m/s";

            testCase.verifyEqual( ...
                u.formattedLabel(),"Particle velocity (m/s)");
        end

        function setValuesArray1D(testCase)
            u = GridFunction(testCase.Grid1DPrimal);
            data = [1;2;3;4;5];

            u.setValues(data);

            testCase.verifyEqual(u.values,data);
        end

        function setValuesRowArray1D(testCase)
            u = GridFunction(testCase.Grid1DPrimal);
            data = [1,2,3,4,5];

            u.setValues(data);

            testCase.verifyEqual(u.values,data.');
        end

        function setValuesFunction1D(testCase)
            u = GridFunction(testCase.Grid1DPrimal);

            u.setValues(@(x) x.^2);

            expected = testCase.XPrimal.pts.^2;
            testCase.verifyEqual(u.values,expected,AbsTol=1e-14);
        end

        function zeroInitialization2D(testCase)
            u = GridFunction(testCase.Grid2DShort);

            testCase.verifyEqual(u.values,zeros(5,3));
            testCase.verifyEqual(u.N,[5,3]);
            testCase.verifyEqual(u.dim,2);
        end

        function setValuesArray2D(testCase)
            u = GridFunction(testCase.Grid2DShort);
            data = reshape(1:15,5,3);

            u.setValues(data);

            testCase.verifyEqual(u.values,data);
        end

        function setValuesFunction2D(testCase)
            grid = testCase.Grid2DShort;
            u = GridFunction(grid);

            u.setValues(@(x,y) x.^2+y);

            [X,Y] = grid.mesh();
            expected = X.^2+Y;
            testCase.verifyEqual(u.values,expected,AbsTol=1e-14);
        end

        function interpolationPrimalToDual1D(testCase)
            uP = GridFunction(testCase.Grid1DPrimal);
            uD = GridFunction(testCase.Grid1DDual);
            uP.setValues(@(x) 2*x+1);

            uD.interpolateFrom(uP);

            expected = 2*testCase.XDual.pts+1;
            testCase.verifyEqual(uD.values,expected,AbsTol=1e-14);
        end

        function interpolationPreservesTargetMetadata(testCase)
            source = GridFunction(testCase.Grid1DPrimal, ...
                Label="Pressure",Units="Pa");
            target = GridFunction(testCase.Grid1DDual, ...
                Label="Interpolated pressure",Units="kPa");
            source.setValues(@(x) 2*x+1);

            target.interpolateFrom(source);

            expected = 2*testCase.XDual.pts+1;
            testCase.verifyEqual(target.values,expected,AbsTol=1e-14);
            testCase.verifyEqual(target.label,"Interpolated pressure");
            testCase.verifyEqual(target.units,"kPa");
        end

        function plot1DUsesFormattedLabel(testCase)
            u = GridFunction(testCase.Grid1DPrimal, ...
                Label="Pressure",Units="Pa");

            fig = figure(Visible='off');
            cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
            ax = axes(fig);
            u.plot();

            testCase.verifyEqual(ax.YLabel.String,'Pressure (Pa)');
        end

        function plot2DLabelsColorbar(testCase)
            u = GridFunction(testCase.Grid2DShort, ...
                Label="Particle velocity",Units="m/s");

            fig = figure(Visible='off');
            cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
            axes(fig);
            u.plot();
            cb = findall(fig,Type='colorbar');

            testCase.verifyNumElements(cb,1);
            testCase.verifyEqual(cb.Label.String,'Particle velocity (m/s)');
        end

        function interpolationDualToPrimal1D(testCase)
            uD = GridFunction(testCase.Grid1DDual);
            uP = GridFunction(testCase.Grid1DPrimal);
            uD.setValues(@(x) 2*x+1);

            uP.interpolateFrom(uD);

            expectedInterior = 2*testCase.XPrimal.pts(2:end-1)+1;
            testCase.verifyEqual( ...
                uP.values(2:end-1),expectedInterior,AbsTol=1e-14);
            testCase.verifyTrue(isnan(uP.values(1)));
            testCase.verifyTrue(isnan(uP.values(end)));
        end

        function interpolation2D(testCase)
            uPP = GridFunction(testCase.Grid2DPrimal);
            uDP = GridFunction(testCase.Grid2DDualPrimal);
            uPP.setValues(@(x,y) 2*x-3*y+x.*y+1);

            uDP.interpolateFrom(uPP);

            [X,Y] = testCase.Grid2DDualPrimal.mesh();
            expected = 2*X-3*Y+X.*Y+1;
            testCase.verifyEqual(uDP.values,expected,AbsTol=1e-13);
        end

        function interpolationDimensionMismatch(testCase)
            grid2D = GridSpace2D(testCase.XPrimal,testCase.XPrimal);
            u1 = GridFunction(testCase.Grid1DPrimal);
            u2 = GridFunction(grid2D);

            testCase.verifyError( ...
                @() u1.interpolateFrom(u2), ...
                'GridFunction:DimensionMismatch');
        end

        function invalidArraySize1D(testCase)
            u = GridFunction(testCase.Grid1DPrimal);

            testCase.verifyError( ...
                @() u.setValues([1;2;3]), ...
                'GridFunction:InvalidSize');
        end

        function invalidArraySize2D(testCase)
            u = GridFunction(testCase.Grid2DShort);

            testCase.verifyError( ...
                @() u.setValues(zeros(3,5)), ...
                'GridFunction:InvalidSize');
        end

        function invalidInputType(testCase)
            u = GridFunction(testCase.Grid1DPrimal);

            testCase.verifyError( ...
                @() u.setValues("invalid"), ...
                'GridFunction:InvalidData');
        end
    end
end
