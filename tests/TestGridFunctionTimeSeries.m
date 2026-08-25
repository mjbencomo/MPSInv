classdef TestGridFunctionTimeSeries < matlab.unittest.TestCase
    % Tests for GridFunctionTimeSeries.

    properties
        XPrimal
        XDual
        YPrimal
        Grid1D
        Grid1DDual
        Grid2D
        Times = [0,0.25,1]
    end

    methods (TestMethodSetup)
        function createGrids(testCase)
            testCase.XPrimal = Grid1DUnifPrimal([0,1],5);
            testCase.XDual = Grid1DUnifDual(testCase.XPrimal);
            testCase.YPrimal = Grid1DUnifPrimal([0,2],3);
            testCase.Grid1D = GridSpace1D(testCase.XPrimal);
            testCase.Grid1DDual = GridSpace1D(testCase.XDual);
            testCase.Grid2D = GridSpace2D( ...
                testCase.XPrimal,testCase.YPrimal);
        end
    end

    methods (Test)
        function zeroInitialization1D(testCase)
            U = GridFunctionTimeSeries(testCase.Grid1D,testCase.Times);

            testCase.verifyEqual(U.values,zeros(5,3));
            testCase.verifyEqual(U.N,5);
            testCase.verifyEqual(U.Nt,3);
            testCase.verifyEqual(U.dim,1);
            testCase.verifyEqual(U.times,testCase.Times);
            testCase.verifyEqual(U.label,"Value");
            testCase.verifyEqual(U.units,"");
            testCase.verifyEqual(U.timeLabel,"Time");
            testCase.verifyEqual(U.timeUnits,"");
        end

        function constructorStoresMetadata(testCase)
            U = GridFunctionTimeSeries( ...
                testCase.Grid1D,testCase.Times, ...
                Label="Pressure",Units="Pa", ...
                TimeLabel="Time",TimeUnits="s");

            testCase.verifyEqual(U.label,"Pressure");
            testCase.verifyEqual(U.units,"Pa");
            testCase.verifyEqual(U.timeLabel,"Time");
            testCase.verifyEqual(U.timeUnits,"s");
            testCase.verifyEqual(U.formattedLabel(),"Pressure (Pa)");
            testCase.verifyEqual(U.formattedTime(3),"Time = 1 s");
        end

        function zeroInitialization2D(testCase)
            U = GridFunctionTimeSeries(testCase.Grid2D,testCase.Times);

            testCase.verifyEqual(U.values,zeros(5,3,3));
            testCase.verifyEqual(U.N,[5,3]);
            testCase.verifyEqual(U.Nt,3);
            testCase.verifyEqual(U.dim,2);
        end

        function constructorConvertsTimesToRow(testCase)
            U = GridFunctionTimeSeries(testCase.Grid1D,[0;0.5;1]);

            testCase.verifyEqual(U.times,[0,0.5,1]);
            testCase.verifySize(U.times,[1,3]);
        end

        function setNumericSnapshot1D(testCase)
            U = GridFunctionTimeSeries(testCase.Grid1D,testCase.Times);
            data = [1,2,3,4,5];

            U.setSnapshot(2,data);

            testCase.verifyEqual(U.values(:,2),data.');
            testCase.verifyEqual(U.values(:,1),zeros(5,1));
            testCase.verifyEqual(U.values(:,3),zeros(5,1));
        end

        function setNumericSnapshot2D(testCase)
            U = GridFunctionTimeSeries(testCase.Grid2D,testCase.Times);
            data = reshape(1:15,5,3);

            U.setSnapshot(3,data);

            testCase.verifyEqual(U.values(:,:,3),data);
            testCase.verifyEqual(U.values(:,:,1),zeros(5,3));
        end

        function snapshotReturnsGridFunction1D(testCase)
            U = GridFunctionTimeSeries( ...
                testCase.Grid1D,testCase.Times, ...
                Label="Pressure",Units="Pa");
            data = (1:5).';
            U.setSnapshot(2,data);

            u = U.snapshot(2);

            testCase.verifyClass(u,'GridFunction');
            testCase.verifyEqual(u.grid,testCase.Grid1D);
            testCase.verifyEqual(u.values,data);
            testCase.verifyEqual(u.label,"Pressure");
            testCase.verifyEqual(u.units,"Pa");
        end

        function snapshotReturnsGridFunction2D(testCase)
            U = GridFunctionTimeSeries(testCase.Grid2D,testCase.Times);
            data = reshape(1:15,5,3);
            U.setSnapshot(2,data);

            u = U.snapshot(2);

            testCase.verifyClass(u,'GridFunction');
            testCase.verifyEqual(u.grid,testCase.Grid2D);
            testCase.verifyEqual(u.values,data);
        end

        function animate1DUpdatesLineThroughFinalTime(testCase)
            U = GridFunctionTimeSeries(testCase.Grid1D,testCase.Times);
            for k = 1:U.Nt
                U.setSnapshot(k,k*testCase.XPrimal.pts);
            end

            fig = figure(Visible='off');
            cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
            ax = axes(fig);
            h = U.animate(Parent=ax,FramePause=0);

            testCase.verifyTrue(isgraphics(h,'line'));
            testCase.verifyEqual(h.YData(:),U.values(:,end));
            testCase.verifyEqual(string(ax.Title.String),"Time = 1");
        end

        function animate2DUpdatesImageThroughFinalTime(testCase)
            U = GridFunctionTimeSeries(testCase.Grid2D,testCase.Times);
            [X,Y] = testCase.Grid2D.mesh();
            for k = 1:U.Nt
                U.setSnapshot(k,X+k*Y);
            end

            fig = figure(Visible='off');
            cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
            ax = axes(fig);
            h = U.animate(Parent=ax,FramePause=0);

            testCase.verifyTrue(isgraphics(h,'image'));
            testCase.verifyEqual(h.CData,U.values(:,:,end).');
            testCase.verifyEqual(ax.YDir,'normal');
            testCase.verifyEqual(string(ax.Title.String),"Time = 1");
        end

        function animateRejectsInvalidParent(testCase)
            U = GridFunctionTimeSeries(testCase.Grid1D,testCase.Times);

            testCase.verifyError( ...
                @() U.animate(Parent=42,FramePause=0), ...
                'GridFunctionTimeSeries:InvalidParent');
        end

        function setSnapshotFromGridFunctionOnSameGrid(testCase)
            U = GridFunctionTimeSeries( ...
                testCase.Grid1D,testCase.Times, ...
                Label="Stored field",Units="stored units");
            u = GridFunction(testCase.Grid1D, ...
                Label="Input field",Units="input units");
            u.setValues(@(x) 2*x+1);

            U.setSnapshot(1,u);

            testCase.verifyEqual(U.values(:,1),u.values);
            testCase.verifyEqual(U.label,"Stored field");
            testCase.verifyEqual(U.units,"stored units");
        end

        function setSnapshotInterpolatesOntoOutputGrid(testCase)
            U = GridFunctionTimeSeries( ...
                testCase.Grid1DDual,testCase.Times);
            u = GridFunction(testCase.Grid1D);
            u.setValues(@(x) 2*x+1);

            U.setSnapshot(2,u);

            expected = 2*testCase.XDual.pts+1;
            testCase.verifyEqual(U.values(:,2),expected,AbsTol=1e-14);
        end

        function constructorRejectsEmptyTimes(testCase)
            testCase.verifyError( ...
                @() GridFunctionTimeSeries(testCase.Grid1D,[]), ...
                'GridFunctionTimeSeries:EmptyTimeGrid');
        end

        function constructorRejectsUnorderedTimes(testCase)
            testCase.verifyError( ...
                @() GridFunctionTimeSeries( ...
                    testCase.Grid1D,[0,0.5,0.25]), ...
                'GridFunctionTimeSeries:InvalidTimes');
        end

        function constructorRejectsNonvectorTimes(testCase)
            testCase.verifyError( ...
                @() GridFunctionTimeSeries(testCase.Grid1D,[0,1;2,3]), ...
                'GridFunctionTimeSeries:InvalidTimes');
        end

        function constructorRejectsRepeatedTimes(testCase)
            testCase.verifyError( ...
                @() GridFunctionTimeSeries(testCase.Grid1D,[0,0.5,0.5]), ...
                'GridFunctionTimeSeries:InvalidTimes');
        end

        function setSnapshotRejectsInvalidTimeIndex(testCase)
            U = GridFunctionTimeSeries(testCase.Grid1D,testCase.Times);

            testCase.verifyError( ...
                @() U.setSnapshot(4,zeros(5,1)), ...
                'GridFunctionTimeSeries:InvalidTimeIndex');
        end

        function snapshotRejectsInvalidTimeIndex(testCase)
            U = GridFunctionTimeSeries(testCase.Grid1D,testCase.Times);

            testCase.verifyError( ...
                @() U.snapshot(4), ...
                'GridFunctionTimeSeries:InvalidTimeIndex');
        end

        function setSnapshotRejectsInvalidSize1D(testCase)
            U = GridFunctionTimeSeries(testCase.Grid1D,testCase.Times);

            testCase.verifyError( ...
                @() U.setSnapshot(1,zeros(4,1)), ...
                'GridFunctionTimeSeries:InvalidSize');
        end

        function setSnapshotRejectsInvalidSize2D(testCase)
            U = GridFunctionTimeSeries(testCase.Grid2D,testCase.Times);

            testCase.verifyError( ...
                @() U.setSnapshot(1,zeros(3,5)), ...
                'GridFunctionTimeSeries:InvalidSize');
        end

        function setSnapshotRejectsInvalidDataType(testCase)
            U = GridFunctionTimeSeries(testCase.Grid1D,testCase.Times);

            testCase.verifyError( ...
                @() U.setSnapshot(1,"invalid"), ...
                'GridFunctionTimeSeries:InvalidData');
        end

        function setSnapshotRejectsDimensionMismatch(testCase)
            U = GridFunctionTimeSeries(testCase.Grid1D,testCase.Times);
            u2D = GridFunction(testCase.Grid2D);

            testCase.verifyError( ...
                @() U.setSnapshot(1,u2D), ...
                'GridFunctionTimeSeries:DimensionMismatch');
        end
    end
end
