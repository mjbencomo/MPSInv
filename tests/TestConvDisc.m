classdef TestConvDisc < matlab.unittest.TestCase
    % Unit tests for convDisc.

    methods (Test)
        function computesKnownRowVectorExample(testCase)
            x = [2,-1,3];
            % Kernel order: [k(-2), k(-1), k(0), k(1), k(2)].
            kernel = [4,-2,1,3,5];

            actual = convDisc(kernel,x);
            expected = TestConvDisc.directDiscreteConvolution(kernel,x);

            testCase.verifyEqual(actual,expected,'AbsTol',1e-12);
            testCase.verifySize(actual,size(x));
        end

        function computesColumnVectorExample(testCase)
            x = [1;3;-2;4];
            kernel = [-1;2;0.5;3;-2;4;1];

            actual = convDisc(kernel,x);
            expected = TestConvDisc.directDiscreteConvolution(kernel,x);

            testCase.verifyEqual(actual,expected,'AbsTol',1e-12);
            testCase.verifySize(actual,size(x));
        end

        function acceptsOppositeInputOrientations(testCase)
            x = [1,2,3];
            kernel = [2;-1;4;3;0.5];

            actual = convDisc(kernel,x);
            expected = TestConvDisc.directDiscreteConvolution(kernel,x);

            testCase.verifyEqual(actual,expected,'AbsTol',1e-12);
            testCase.verifySize(actual,size(x));
        end

        function centeredImpulseKernelReturnsInput(testCase)
            x = [3;-1;4;2];
            n = numel(x);
            kernel = zeros(2*n-1,1);
            kernel(n) = 1;  % k(0) = 1

            testCase.verifyEqual(convDisc(kernel,x),x,'AbsTol',1e-12);
        end

        function handlesSingleEntryInputs(testCase)
            testCase.verifyEqual(convDisc(3,4),12,'AbsTol',1e-12);
        end

        function agreesWithDirectSummation(testCase)
            n = 8;
            kernel = cos((1:2*n-1).');
            x = sin((1:n).');

            actual = convDisc(kernel,x);
            expected = TestConvDisc.directDiscreteConvolution(kernel,x);

            testCase.verifyEqual(actual,expected,'AbsTol',1e-11);
        end

        function rejectsKernelWithWrongLength(testCase)
            x = [1,2,3];
            kernel = [1,2,3,4];

            testCase.verifyError(@() convDisc(kernel,x), ...
                'convDisc:SizeMismatch');
        end
    end

    methods (Static, Access = private)
        function y = directDiscreteConvolution(kernel,x)
            outputIsRow = isrow(x);
            kernel = kernel(:);
            x = x(:);
            n = numel(x);
            y = zeros(n,1,'like',kernel + x(1));

            % kernel(n) stores k(0), so kernel(i-j+n) stores k(i-j).
            for i = 1:n
                for j = 1:n
                    y(i) = y(i) + kernel(i-j+n)*x(j);
                end
            end

            if outputIsRow
                y = y.';
            end
        end
    end
end
