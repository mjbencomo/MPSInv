classdef TestCorrMulti < matlab.unittest.TestCase
    % Channel order: kernel(time,forwardOutput,forwardInput).
    % Dense block-matrix oracle is independent of all convolution routines.
    properties (TestParameter)
        geometry = struct('scalar',[1 1 1], 'oneTime',[1 2 3], ...
            'singleChannel',[5 1 1], 'oneOutput',[4 1 3], ...
            'oneInput',[5 3 1], 'rectangular',[8 2 3]); % [K N M]
        precision = {'double','single'};
    end
    methods (Test)
        function matchesBlockTranspose(t,geometry,precision)
            K = geometry(1); N = geometry(2); M = geometry(3);
            g = cast(reshape(sin(1:(2*K-1)*N*M),2*K-1,N,M),precision);
            Y = cast(reshape(cos(1:K*N),K,N),precision);
            A = TestCorrMulti.matrix(double(g),K,N,M);
            expected = cast(reshape(A.'*double(Y(:)),K,M),precision);
            actual = corrMulti(g,Y);
            t.verifySize(actual,[K M]);
            t.verifyClass(actual,precision);
            t.verifyTrue(isreal(actual));
            t.verifyEqual(actual,expected,'AbsTol', ...
                cast(100*K*N*eps(precision)*max(1,norm(double(expected(:)),inf)),precision));
        end
        function satisfiesAdjointIdentity(t,geometry)
            K = geometry(1); N = geometry(2); M = geometry(3);
            g = reshape(cos(1:(2*K-1)*N*M),2*K-1,N,M);
            X = reshape(sin(1:K*M),K,M);
            Y = reshape(cos(2*(1:K*N)),K,N);
            AX = convMulti(g,X);
            ATY = corrMulti(g,Y);
            scale = max([1,norm(AX(:))*norm(Y(:)),norm(X(:))*norm(ATY(:))]);
            t.verifyEqual(AX(:).'*Y(:),X(:).'*ATY(:), ...
                'AbsTol',100*K*max(N,M)*eps*scale);
        end
        function centeredKernelsApplyTransposedChannelMap(t)
            K = 4;
            gains = [2 -1 0.5;0 3 -2]; % forward map N-by-M
            g = zeros(2*K-1,2,3);
            g(K,:,:) = reshape(gains,1,2,3);
            Y = [1 2;3 -2;-1 1;2 0];
            t.verifyEqual(corrMulti(g,Y),Y*gains,'AbsTol',1e-12);
        end
        function isolatedChannelAndLagDoNotLeak(t)
            g = zeros(7,2,3);
            g(5,2,3) = 2; % only output channel 2 -> input channel 3
            Y = [100 2;200 -3;300 5;400 7];
            expected = zeros(4,3);
            expected(:,3) = [-6;10;14;0];
            t.verifyEqual(corrMulti(g,Y),expected,'AbsTol',1e-12);
        end
        function zeroInputs(t)
            g = reshape(sin(1:42),7,2,3);
            Y = reshape(cos(1:8),4,2);
            t.verifyEqual(corrMulti(0*g,Y),zeros(4,3),'AbsTol',1e-12);
            t.verifyEqual(corrMulti(g,0*Y),zeros(4,3),'AbsTol',1e-12);
        end
        function rejectsWrongTimeDimension(t)
            t.verifyError(@() corrMulti(zeros(4,2,3),zeros(3,2)), ...
                'corrMulti:TimeSizeMismatch');
        end
        function rejectsWrongChannelDimension(t)
            % Y has M columns instead of N; unequal N and M expose swaps.
            t.verifyError(@() corrMulti(zeros(5,2,3),zeros(3,3)), ...
                'corrMulti:InputSizeMismatch');
        end
        function rejectsHigherDimensionalData(t)
            t.verifyError(@() corrMulti(zeros(5,2,3),zeros(3,2,2)), ...
                'corrMulti:InvalidInputDimensions');
        end
        function rejectsHigherDimensionalKernel(t)
            t.verifyError(@() corrMulti(zeros(5,2,3,2),zeros(3,2)), ...
                'corrMulti:InvalidKernelDimensions');
        end
        function rejectsInvalidArguments(t)
            g = ones(5,2,3); Y = ones(3,2);
            invalid = {[],1+1i,{1}};
            for k = 1:numel(invalid)
                bad = invalid{k};
                t.verifyTrue(TestCorrMulti.throws(@() corrMulti(bad,Y)));
                t.verifyTrue(TestCorrMulti.throws(@() corrMulti(g,bad)));
            end
        end
    end
    methods (Static, Access = private)
        function A = matrix(g,K,N,M)
            % vec(convMulti(g,X)) = A*vec(X), column-major ordering.
            A = zeros(K*N,K*M);
            for n = 1:N
                for m = 1:M
                    for i = 1:K
                        for j = 1:K
                            A((n-1)*K+i,(m-1)*K+j) = g(K+i-j,n,m);
                        end
                    end
                end
            end
        end
        function tf = throws(f)
            tf = false;
            try
                f();
            catch
                tf = true;
            end
        end
    end
end
