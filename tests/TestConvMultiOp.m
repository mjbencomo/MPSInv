classdef TestConvMultiOp < matlab.unittest.TestCase
    % Tests for ConvMultiOp and the LinearOp interface.

    methods (Test)
        function constructorSetsDimensionsAndChannelCounts(testCase)
            K = 4;
            N = 2;
            M = 3;
            kernels = zeros(2*K-1,N,M);

            op = ConvMultiOp(kernels,Scale=0.25);

            testCase.verifyTrue(isa(op,'LinearOp'));
            testCase.verifyEqual(op.dimDomain,K*M);
            testCase.verifyEqual(op.dimRange,K*N);
            testCase.verifyEqual(op.N_time,K);
            testCase.verifyEqual(op.N_input,M);
            testCase.verifyEqual(op.N_output,N);
            testCase.verifyEqual(op.scale,0.25);
        end

        function forwardActionAgreesWithConvMulti(testCase)
            K = 4;
            N = 2;
            M = 3;
            scale = 0.2;
            X = [1, 2,-1; ...
                 3,-2, 4; ...
                -1, 1, 2; ...
                 2, 0, 3];
            kernels = TestConvMultiOp.makeKernels(K,N,M);
            op = ConvMultiOp(kernels,Scale=scale);

            actual = op.fwd(X(:));
            expected = scale*convMulti(kernels,X);

            testCase.verifyEqual(actual,expected(:),'AbsTol',1e-11);
        end

        function adjointSatisfiesDotProductIdentity(testCase)
            rng(7);
            K = 5;
            N = 3;
            M = 2;
            kernels = randn(2*K-1,N,M);
            op = ConvMultiOp(kernels,Scale=0.125);
            x = randn(op.dimDomain,1);
            y = randn(op.dimRange,1);

            lhs = dot(op.fwd(x),y);
            rhs = dot(x,op.adj(y));

            tolerance = 1e3*eps*max([1,abs(lhs),abs(rhs)]);
            testCase.verifyEqual(lhs,rhs,'AbsTol',tolerance);
        end

        function matrixRepresentationMatchesActions(testCase)
            rng(11);
            K = 4;
            N = 2;
            M = 3;
            kernels = randn(2*K-1,N,M);
            op = ConvMultiOp(kernels,Scale=0.05);
            A = op.matrixRep();
            x = randn(op.dimDomain,1);
            y = randn(op.dimRange,1);

            testCase.verifySize(A,[op.dimRange,op.dimDomain]);
            testCase.verifyEqual(op.fwd(x),A*x,'AbsTol',1e-12);
            testCase.verifyEqual(op.adj(y),A.'*y,'AbsTol',1e-12);
        end

        function matrixColumnsAreStandardBasisActions(testCase)
            K = 3;
            kernels = TestConvMultiOp.makeKernels(K,2,2);
            op = ConvMultiOp(kernels);
            A = op.matrixRep();

            for j = 1:op.dimDomain
                e = zeros(op.dimDomain,1);
                e(j) = 1;
                testCase.verifyEqual(A(:,j),op.fwd(e), ...
                    'AbsTol',1e-12);
            end
        end

        function centeredKernelsGiveInstantaneousMixing(testCase)
            K = 4;
            gains = [2,-1,0.5; 0,3,-2];
            kernels = zeros(2*K-1,2,3);
            kernels(K,:,:) = reshape(gains,[1,2,3]);
            X = reshape(1:K*3,K,3);
            op = ConvMultiOp(kernels,Scale=0.5);

            expected = 0.5*X*gains.';
            actual = reshape(op.fwd(X(:)),K,2);

            testCase.verifyEqual(actual,expected,'AbsTol',1e-12);
        end

        function handlesSingleInputAndOutput(testCase)
            K = 4;
            kernel = (1:2*K-1).';
            x = [2;-1;3;4];
            op = ConvMultiOp(kernel,Scale=0.1);

            testCase.verifyEqual(op.fwd(x), ...
                0.1*convDisc(kernel,x),'AbsTol',1e-12);
            testCase.verifyEqual(op.dimDomain,K);
            testCase.verifyEqual(op.dimRange,K);
        end

        function rejectsEvenKernelTimeDimension(testCase)
            testCase.verifyError(@() ConvMultiOp(zeros(6,2,3)), ...
                'ConvMultiOp:InvalidKernelTimeSize');
        end

        function rejectsKernelWithTooManyDimensions(testCase)
            testCase.verifyError(@() ConvMultiOp(zeros(5,2,2,2)), ...
                'ConvMultiOp:InvalidKernelDimensions');
        end

        function rejectsIncorrectDomainVector(testCase)
            op = ConvMultiOp(zeros(5,2,2));

            testCase.verifyError(@() op.fwd(zeros(op.dimDomain-1,1)), ...
                'LinearOp:InvalidDomainVector');
            testCase.verifyError(@() op.fwd(zeros(1,op.dimDomain)), ...
                'LinearOp:InvalidDomainVector');
        end

        function rejectsIncorrectRangeVector(testCase)
            op = ConvMultiOp(zeros(5,2,2));

            testCase.verifyError(@() op.adj(zeros(op.dimRange-1,1)), ...
                'LinearOp:InvalidRangeVector');
            testCase.verifyError(@() op.adj(zeros(1,op.dimRange)), ...
                'LinearOp:InvalidRangeVector');
        end
    end

    methods (Static, Access = private)
        function kernels = makeKernels(K,N,M)
            kernels = zeros(2*K-1,N,M);
            for n = 1:N
                for m = 1:M
                    kernels(:,n,m) = ...
                        (n+2*m)*(1:2*K-1).'-(m-n);
                end
            end
        end
    end
end
