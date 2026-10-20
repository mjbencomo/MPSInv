classdef TestReceiverSpace < matlab.unittest.TestCase
    methods (Test)
        function oneDimensionalOrderAndRepeatedLocations(testCase)
            x = [0.8; 0.2; 0.8];
            receivers = ReceiverSpace1D(x);

            testCase.verifyTrue(isa(receivers,'GridSpace'));
            testCase.verifyEqual(receivers.dim,1);
            testCase.verifyEqual(receivers.numPts,3);
            testCase.verifyEqual(receivers.points,x);
        end

        function twoDimensionalPairsAndRepeatedLocations(testCase)
            xy = [0.8 0.2; 0.1 0.9; 0.8 0.2];
            receivers = ReceiverSpace2D(xy);

            testCase.verifyTrue(isa(receivers,'GridSpace'));
            testCase.verifyEqual(receivers.dim,2);
            testCase.verifyEqual(receivers.numPts,3);
            testCase.verifyEqual(receivers.points,xy);
        end

        function singleReceiver(testCase)
            testCase.verifyEqual(ReceiverSpace1D(0.5).points,0.5);
            testCase.verifyEqual(ReceiverSpace2D([0.5 0.7]).numPts,1);
        end

        function invalidCoordinates(testCase)
            testCase.verifyError(@() ReceiverSpace1D([]), ...
                'ReceiverSpace1D:InvalidPoints');
            testCase.verifyError(@() ReceiverSpace2D([1;2]), ...
                'ReceiverSpace2D:InvalidPoints');
            testCase.verifyError(@() ReceiverSpace2D(zeros(0,2)), ...
                'ReceiverSpace2D:InvalidPoints');
        end

        function cartesianPointCountIsPreserved(testCase)
            x = Grid1DUnifPrimal([0,1],5);
            y = Grid1DUnifPrimal([0,1],4);

            testCase.verifyEqual(GridSpace1D(x).numPts,5);
            testCase.verifyEqual(GridSpace2D(x,y).numPts,20);
        end
    end
end
