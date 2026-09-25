classdef TestConvMulti < matlab.unittest.TestCase
    % Unit tests for convMulti.

    methods (Test)
        function singleInputSingleOutputAgreesWithConvDisc(testCase)
            X = [2;-1;3;4];
            kernel = [1;-2;0.5;3;-1;4;2];

            actual = convMulti(kernel,X);
            expected = convDisc(kernel,X);

            testCase.verifyEqual(actual,expected,'AbsTol',1e-12);
            testCase.verifySize(actual,[4,1]);
        end

        function computesMultipleInputAndOutputChannels(testCase)
            K = 4;
            N = 2;
            M = 3;
            X = [1, 2,-1; ...
                 3,-2, 4; ...
                -1, 1, 2; ...
                 2, 0, 3];
            kernel = zeros(2*K-1,N,M);

            for n = 1:N
                for m = 1:M
                    kernel(:,n,m) = (n+2*m)*(1:2*K-1).' - (m-n);
                end
            end

            actual = convMulti(kernel,X);
            expected = TestConvMulti.directMultiConvolution(kernel,X);

            testCase.verifyEqual(actual,expected,'AbsTol',1e-11);
            testCase.verifySize(actual,[K,N]);
        end

        function centeredKernelsMixChannelsInstantaneously(testCase)
            K = 5;
            N = 2;
            M = 3;
            X = [1, 2, 3; ...
                 4,-1, 2; ...
                 0, 3,-2; ...
                 2, 1, 5; ...
                -1, 4, 1];
            gains = [2,-1,0.5; ...
                     0, 3, -2];
            kernel = zeros(2*K-1,N,M);

            for n = 1:N
                for m = 1:M
                    kernel(K,n,m) = gains(n,m);
                end
            end

            expected = X*gains.';
            testCase.verifyEqual(convMulti(kernel,X),expected, ...
                'AbsTol',1e-12);
        end

        function handlesOneOutputAndMultipleInputs(testCase)
            K = 3;
            M = 2;
            X = [1,4; 2,-1; 3,2];
            kernel = zeros(2*K-1,1,M);
            kernel(:,1,1) = [1;0;-1;2;3];
            kernel(:,1,2) = [-2;1;4;0;2];

            actual = convMulti(kernel,X);
            expected = TestConvMulti.directMultiConvolution(kernel,X);

            testCase.verifyEqual(actual,expected,'AbsTol',1e-12);
            testCase.verifySize(actual,[K,1]);
        end

        function handlesOneInputAndMultipleOutputs(testCase)
            K = 3;
            X = [2;-1;4];
            kernel = [1,-2; ...
                      0, 3; ...
                      2, 1; ...
                     -1, 4; ...
                      3, 0];

            actual = convMulti(kernel,X);
            expected = TestConvMulti.directMultiConvolution(kernel,X);

            testCase.verifyEqual(actual,expected,'AbsTol',1e-12);
            testCase.verifySize(actual,[K,2]);
        end

        function handlesSingleTimePoint(testCase)
            X = [2,3];
            kernel = reshape([1,2,4,-1],[1,2,2]);

            actual = convMulti(kernel,X);
            expected = TestConvMulti.directMultiConvolution(kernel,X);

            testCase.verifyEqual(actual,expected,'AbsTol',1e-12);
            testCase.verifySize(actual,[1,2]);
        end

        function preservesSinglePrecision(testCase)
            X = single([1;2;3]);
            kernel = single([1;0;2;-1;3]);

            actual = convMulti(kernel,X);

            testCase.verifyClass(actual,'single');
            testCase.verifyEqual(actual,convDisc(kernel,X), ...
                'AbsTol',single(1e-5));
        end

        function rejectsIncorrectKernelTimeDimension(testCase)
            X = zeros(3,2);
            kernel = zeros(4,2,2);

            testCase.verifyError(@() convMulti(kernel,X), ...
                'convMulti:TimeSizeMismatch');
        end

        function rejectsIncorrectNumberOfInputChannels(testCase)
            X = zeros(3,2);
            kernel = zeros(5,2,3);

            testCase.verifyError(@() convMulti(kernel,X), ...
                'convMulti:InputSizeMismatch');
        end

        function rejectsHigherDimensionalInput(testCase)
            X = zeros(3,2,2);
            kernel = zeros(5,2,2);

            testCase.verifyError(@() convMulti(kernel,X), ...
                'convMulti:InvalidInputDimensions');
        end

        function rejectsKernelWithMoreThanThreeDimensions(testCase)
            X = zeros(3,2);
            kernel = zeros(5,2,2,2);

            testCase.verifyError(@() convMulti(kernel,X), ...
                'convMulti:InvalidKernelDimensions');
        end
    end

    methods (Static, Access = private)
        function Y = directMultiConvolution(kernel,X)
            K = size(X,1);
            N = size(kernel,2);
            M = size(kernel,3);
            Y = zeros(K,N,'like',X);

            % kernel(K,n,m) stores the zero-lag value for channel (n,m).
            for n = 1:N
                for m = 1:M
                    for i = 1:K
                        for j = 1:K
                            lagIndex = i-j+K;
                            Y(i,n) = Y(i,n) + ...
                                kernel(lagIndex,n,m)*X(j,m);
                        end
                    end
                end
            end
        end
    end
end
