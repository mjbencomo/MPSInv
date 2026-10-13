classdef TestPlotAcqGeoReceivers < matlab.unittest.TestCase
    methods (Test)
        function plotsOneDimensionalNumericLocations(testCase)
            grid = GridSpace1D(Grid1DUnifPrimal([0,1],5));
            background = GridFunction(grid);
            background.setValues(@(x) x);
            sourceX = [0.7;0.2];
            receiverX = [0.8;0.3;0.8];

            fig = figure(Visible='off');
            cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
            h = plotAcqGeo(background,sourceX,receiverX,Parent=axes(fig));

            testCase.verifyEqual(h.src_rcv(1).XData(:),sourceX);
            testCase.verifyEqual(h.src_rcv(2).XData(:),receiverX);
        end

        function plotsPairedTwoDimensionalNumericLocations(testCase)
            x = Grid1DUnifPrimal([0,1],5);
            background = GridFunction(GridSpace2D(x,x));
            background.setValues(@(X,Y) X+Y);
            sourceXY = [0.1 0.8;0.7 0.2];
            receiverXY = [0.6 0.1;0.2 0.9;0.6 0.1];

            fig = figure(Visible='off');
            cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
            h = plotAcqGeo(background,sourceXY,receiverXY, ...
                Parent=axes(fig));

            testCase.verifyEqual(h.src_rcv(1).XData(:),sourceXY(:,1));
            testCase.verifyEqual(h.src_rcv(1).YData(:),sourceXY(:,2));
            testCase.verifyEqual(h.src_rcv(2).XData(:),receiverXY(:,1));
            testCase.verifyEqual(h.src_rcv(2).YData(:),receiverXY(:,2));
            testCase.verifyEqual(numel(h.src_rcv(2).XData),3);
        end

        function rejectsSourceObjectAndComponentCell(testCase)
            x = Grid1DUnifPrimal([0,1],5);
            background = GridFunction(GridSpace2D(x,x));
            component = MultipoleComponent2D([0.2 0.7],[0,0]);

            testCase.verifyError(@() plotAcqGeo( ...
                background,component,[0.5 0.5]), ...
                'plotAcqGeo:InvalidSource');
            testCase.verifyError(@() plotAcqGeo( ...
                background,{component},[0.5 0.5]), ...
                'plotAcqGeo:InvalidSource');
        end

        function rejectsReceiverSpaceObject(testCase)
            x = Grid1DUnifPrimal([0,1],5);
            background = GridFunction(GridSpace2D(x,x));
            receivers = ReceiverSpace2D([0.2 0.8;0.7 0.3]);

            testCase.verifyError(@() plotAcqGeo( ...
                background,[0.5 0.5],receivers), ...
                'plotAcqGeo:InvalidReceivers');
        end

        function rejectsWrongShapeAndNonfiniteLocations(testCase)
            x = Grid1DUnifPrimal([0,1],5);
            background = GridFunction(GridSpace2D(x,x));

            testCase.verifyError(@() plotAcqGeo( ...
                background,[0.2;0.3],[0.5 0.5]), ...
                'plotAcqGeo:InvalidSource');
            testCase.verifyError(@() plotAcqGeo( ...
                background,[0.2 0.3],[NaN 0.5]), ...
                'plotAcqGeo:InvalidReceivers');
        end

        function allowsEmptySourceLocations(testCase)
            x = Grid1DUnifPrimal([0,1],5);
            background = GridFunction(GridSpace2D(x,x));

            fig = figure(Visible='off');
            cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
            h = plotAcqGeo(background,zeros(0,2),[0.3 0.7], ...
                Parent=axes(fig));

            testCase.verifyEqual(numel(h.src_rcv),2);
            testCase.verifyEmpty(h.src_rcv(1).XData);
            testCase.verifyEqual(h.src_rcv(2).XData,0.3);
            testCase.verifyEqual(h.src_rcv(2).YData,0.7);
        end

        function allowsEmptyReceiverLocations(testCase)
            grid = GridSpace1D(Grid1DUnifPrimal([0,1],5));
            background = GridFunction(grid);

            fig = figure(Visible='off');
            cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
            h = plotAcqGeo(background,[0.2;0.8],[],Parent=axes(fig));

            testCase.verifyEqual(numel(h.src_rcv),2);
            testCase.verifyEqual(h.src_rcv(1).XData(:),[0.2;0.8]);
            testCase.verifyEmpty(h.src_rcv(2).XData);
        end
    end
end
