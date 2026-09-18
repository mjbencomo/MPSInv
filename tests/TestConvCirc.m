classdef TestConvCirc < matlab.unittest.TestCase
    % Unit tests for convCirc.

    methods (Test)
        function computesKnownRowVectorExample(testCase)
            kernel = [1,2,3,4];
            x = [5,-1,2,3];

            actual = convCirc(kernel,x);
            expected = TestConvCirc.directCircularConvolution(kernel,x);

            testCase.verifyEqual(actual,expected,'AbsTol',1e-12);
            testCase.verifySize(actual,size(x));
        end

        function computesColumnVectorExample(testCase)
            kernel = [2;-1;3;0.5];
            x = [1;4;-2;3];

            actual = convCirc(kernel,x);
            expected = TestConvCirc.directCircularConvolution(kernel,x);

            testCase.verifyEqual(actual,expected,'AbsTol',1e-12);
            testCase.verifySize(actual,size(x));
        end

        function impulseKernelReturnsInput(testCase)
            x = [3,-2,5,7,-1];
            kernel = [1,zeros(1,numel(x)-1)];

            testCase.verifyEqual(convCirc(kernel,x),x,'AbsTol',1e-12);
        end

        function circularConvolutionIsCommutative(testCase)
            kernel = [1.5,-2,0.25,4];
            x = [-3,1,2,0.5];

            testCase.verifyEqual(convCirc(kernel,x),convCirc(x,kernel), ...
                'AbsTol',1e-12);
        end

        function rejectsDifferentLengths(testCase)
            testCase.verifyError(@() convCirc([1,2,3],[4,5]), ...
                'convCirc:SizeMismatch');
        end

        function rejectsDifferentOrientations(testCase)
            kernel = [1,2,3];
            x = [4;5;6];

            testCase.verifyError(@() convCirc(kernel,x), ...
                'convCirc:SizeMismatch');
        end
    end

    methods (Static, Access = private)
        function y = directCircularConvolution(kernel,x)
            outputIsRow = isrow(x);
            kernel = kernel(:);
            x = x(:);
            n = numel(x);
            y = zeros(n,1,'like',kernel + x);

            for i = 1:n
                for j = 1:n
                    kernelIndex = mod(i-j,n) + 1;
                    y(i) = y(i) + kernel(kernelIndex)*x(j);
                end
            end

            if outputIsRow
                y = y.';
            end
        end
    end
end
