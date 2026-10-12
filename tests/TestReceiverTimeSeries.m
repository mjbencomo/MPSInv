classdef TestReceiverTimeSeries < matlab.unittest.TestCase
    methods (Test)
        function storesOneDimensionalReceiversInInputOrder(testCase)
            space = ReceiverSpace1D([0.8;0.2;0.8]);
            U = GridFunctionTimeSeries(space,[0,0.5,1], ...
                Label="Pressure",Units="Pa");

            testCase.verifySize(U.values,[3,3]);
            testCase.verifyEqual(U.numPts,3);
            testCase.verifyEqual(U.dim,1);
            U.setSnapshot(2,[7,8,9]);
            testCase.verifyEqual(U.values(:,2),[7;8;9]);

            u = U.snapshot(2);
            testCase.verifyClass(u,'GridFunction');
            testCase.verifyEqual(u.grid,space);
            testCase.verifyEqual(u.values,[7;8;9]);
            testCase.verifyEqual(u.label,"Pressure");
            testCase.verifyEqual(u.units,"Pa");
        end

        function storesTwoDimensionalReceiversAsChannels(testCase)
            space = ReceiverSpace2D([0.8,0.2;0.1,0.9;0.8,0.2]);
            U = GridFunctionTimeSeries(space,[0,1]);

            testCase.verifySize(U.values,[3,2]);
            testCase.verifyEqual(U.dim,2);
            testCase.verifyEqual(U.numPts,3);
            U.setSnapshot(1,[3;4;5]);
            snapshot = U.snapshot(1);
            testCase.verifyEqual(snapshot.values,[3;4;5]);

            u = GridFunction(space);
            u.setValues(@(x,y) x+2*y);
            U.setSnapshot(2,u);
            testCase.verifyEqual(U.values(:,2),[1.2;1.9;1.2], ...
                AbsTol=1e-14);
        end

        function singleTwoDimensionalReceiver(testCase)
            space = ReceiverSpace2D([0.2,0.7]);
            U = GridFunctionTimeSeries(space,[0,1]);
            U.setSnapshot(2,4);
            testCase.verifySize(U.values,[1,2]);
            snapshot = U.snapshot(2);
            testCase.verifyEqual(snapshot.values,4);
        end

        function invalidReceiverSnapshotShape(testCase)
            space = ReceiverSpace2D([0,0;0,1;1,0;1,1]);
            U = GridFunctionTimeSeries(space,[0,1]);

            testCase.verifyError(@() U.setSnapshot(1,ones(2,2)), ...
                'GridFunctionTimeSeries:InvalidSize');
            u = GridFunction(space);
            testCase.verifyError(@() u.setValues(ones(2,2)), ...
                'GridFunction:InvalidSize');
        end

        function cartesianShapeIsStillAvailable(testCase)
            x = Grid1DUnifPrimal([0,1],5);
            y = Grid1DUnifPrimal([0,1],4);
            U = GridFunctionTimeSeries(GridSpace2D(x,y),[0,1]);

            testCase.verifyEqual(U.N,[5,4]);
            testCase.verifyEqual(U.numPts,20);
            testCase.verifySize(U.values,[5,4,2]);
        end

        function receiverGridFunctionFromOtherSpaceIsRejected(testCase)
            target = ReceiverSpace1D([0.2;0.7]);
            source = GridSpace1D(Grid1DUnifPrimal([0,1],5));
            U = GridFunctionTimeSeries(target,[0,1]);
            u = GridFunction(source);

            testCase.verifyError(@() U.setSnapshot(1,u), ...
                'GridFunctionTimeSeries:ReceiverInterpolationUnsupported');
        end
    end
end
