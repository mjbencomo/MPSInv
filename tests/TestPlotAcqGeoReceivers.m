classdef TestPlotAcqGeoReceivers < matlab.unittest.TestCase
    methods (Test)
        function plotsOneDimensionalReceiverOrder(testCase)
            grid = GridSpace1D(Grid1DUnifPrimal([0,1],5));
            background = GridFunction(grid);
            background.setValues(@(x) x);
            source = MultipoleSrcTerm1D( ...
                MultipoleComponent1D(0.5,0),@(t) 1);
            receivers = ReceiverSpace1D([0.8;0.2;0.8]);

            fig = figure(Visible='off');
            cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
            ax = axes(fig);
            h = plotAcqGeo(background,source,receivers,Parent=ax);

            testCase.verifyEqual(h.src_rcv(1).XData,0.5);
            testCase.verifyEqual(h.src_rcv(2).XData(:),receivers.points);
        end

        function plotsPairedTwoDimensionalReceivers(testCase)
            x = Grid1DUnifPrimal([0,1],5);
            grid = GridSpace2D(x,x);
            background = GridFunction(grid);
            background.setValues(@(X,Y) X+Y);
            source = MultipoleSrcTerm2D( ...
                MultipoleComponent2D([0.5 0.5],[0,0]),@(t) 1);
            xy = [0.1 0.8;0.7 0.2;0.1 0.8];
            receivers = ReceiverSpace2D(xy);

            fig = figure(Visible='off');
            cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
            ax = axes(fig);
            h = plotAcqGeo(background,source,receivers,Parent=ax);

            testCase.verifyEqual(h.src_rcv(2).XData(:),xy(:,1));
            testCase.verifyEqual(h.src_rcv(2).YData(:),xy(:,2));
            testCase.verifyEqual(numel(h.src_rcv(2).XData),3);
        end

        function preservesCartesianReceiverPlotting(testCase)
            x = Grid1DUnifPrimal([0,1],5);
            grid = GridSpace2D(x,x);
            background = GridFunction(grid);
            source = MultipoleSrcTerm2D( ...
                MultipoleComponent2D([0.5 0.5],[0,0]),@(t) 1);
            receiverGrid = GridSpace2D( ...
                Grid1DUnifPrimal([0,1],3), ...
                Grid1DUnifPrimal([0,1],2));
            [X,Y] = receiverGrid.mesh();

            fig = figure(Visible='off');
            cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
            ax = axes(fig);
            h = plotAcqGeo(background,source,receiverGrid,Parent=ax);

            testCase.verifyEqual(h.src_rcv(2).XData(:),X(:));
            testCase.verifyEqual(h.src_rcv(2).YData(:),Y(:));
        end
    end
end
