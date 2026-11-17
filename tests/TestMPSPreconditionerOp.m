classdef TestMPSPreconditionerOp < matlab.unittest.TestCase
    methods (Test)
        function metadata(testCase)
            s = [0,0;1,0;0,1;1,1];
            B = MPSPreconditionerOp(2,2,s,8,0.2);
            testCase.verifyClass(B,'MPSPreconditionerOp');
            testCase.verifyTrue(isa(B,'InvLinearOp'));
            testCase.verifyEqual(B.dimDomain,32);
            testCase.verifyEqual(B.dimRange,32);
            testCase.verifyEqual(B.s,s);
            testCase.verifyEqual(B.beta,0.5);
            testCase.verifyEqual(B.alpha,[0.5;1.5;1.5;2.5]);
        end

        function actionsMatchIndependentMatrices(testCase)
            % Build Q from products, independently of fractional functions.
            rng(7);
            Nt = 9; dt = 0.3; c = 1.7;
            for d = 1:3
                s = [zeros(1,d);eye(d);2*ones(1,d)];
                % Keep orders moderate for the inverse comparison.
                s(end,:) = 0; s(end,1) = 2;
                B = MPSPreconditionerOp(c,d,s,Nt,dt);
                Q = TestMPSPreconditionerOp.explicitQ(c,d,s,Nt,dt);
                A = Q'*Q;
                x = randn(size(A,1),1) + 1i*randn(size(A,1),1);
                expectedInverse = Q\(Q'\x);
                testCase.verifyEqual(B.fwd(x),A*x, ...
                    'RelTol',1e-11,'AbsTol',1e-10);
                testCase.verifyEqual(B.adj(x),A'*x, ...
                    'RelTol',1e-11,'AbsTol',1e-10);
                testCase.verifyEqual(B.inv(x),expectedInverse, ...
                    'RelTol',1e-10,'AbsTol',1e-10);
                testCase.verifyEqual(B.invAdj(x),expectedInverse, ...
                    'RelTol',1e-10,'AbsTol',1e-10);
                testCase.verifyEqual(B.matrixRep(),A, ...
                    'RelTol',1e-11,'AbsTol',1e-10);
                testCase.verifyEqual(B.inv(B.fwd(x)),x, ...
                    'RelTol',1e-8,'AbsTol',1e-8);
                testCase.verifyEqual(B.fwd(B.inv(x)),x, ...
                    'RelTol',1e-8,'AbsTol',1e-8);
                [~,flag] = chol(B.matrixRep());
                testCase.verifyEqual(flag,0);
            end
        end

        function adjointIdentity(testCase)
            rng(9);
            B = MPSPreconditionerOp(2,2,[0,0;1,0;0,1],10,0.2);
            x = randn(30,1)+1i*randn(30,1);
            y = randn(30,1)+1i*randn(30,1);
            testCase.verifyEqual(y'*B.fwd(x),B.adj(y)'*x, ...
                'RelTol',1e-12,'AbsTol',1e-10);
            testCase.verifyGreaterThan(real(x'*B.fwd(x)),0);
        end

        function zeroOrderAndSingleSample(testCase)
            B = MPSPreconditionerOp(3,1,[0;0],4,0.2);
            x = (1:8).';
            testCase.verifyEqual(B.fwd(x),x);
            testCase.verifyEqual(B.inv(x),x);
            B = MPSPreconditionerOp(2,2,[0,0;1,0],1,0.25);
            x = [3;-2];
            scale = [0.25^(-1);2^(-2)*0.25^(-3)];
            testCase.verifyEqual(B.fwd(x),scale.*x,'AbsTol',1e-12);
            testCase.verifyEqual(B.inv(x),x./scale,'AbsTol',1e-12);
        end

        function channelIsolation(testCase)
            B = MPSPreconditionerOp(2,2,[0,0;1,0;0,1],6,0.2);
            X = zeros(6,3); X(:,2) = (1:6).';
            Y = reshape(B.fwd(X(:)),6,3);
            Z = reshape(B.inv(X(:)),6,3);
            testCase.verifyEqual(Y(:,[1,3]),zeros(6,2));
            testCase.verifyEqual(Z(:,[1,3]),zeros(6,2));
        end

        function rejectsInvalidShapes(testCase)
            testCase.verifyError(@() MPSPreconditionerOp(2,2,[0;1],4,0.2), ...
                'MPSPreconditionerOp:InvalidMultiIndices');
            B = MPSPreconditionerOp(2,2,[0,0],4,0.2);
            testCase.verifyError(@() B.fwd(ones(1,4)), ...
                'LinearOp:InvalidDomainVector');
            testCase.verifyError(@() B.inv(ones(3,1)), ...
                'LinearOp:InvalidRangeVector');
        end
    end

    methods (Static, Access = private)
        function Q = explicitQ(c,d,s,Nt,dt)
            blocks = cell(1,size(s,1));
            for m = 1:size(s,1)
                r = sum(s(m,:));
                a = (d-1)/2+r;
                w = ones(Nt,1);
                for k = 1:Nt-1
                    w(k+1) = prod(((0:k-1)-a)./(1:k));
                end
                D = tril(toeplitz(w,[1,zeros(1,Nt-1)]));
                blocks{m} = c^(-r)*dt^(-a)*D;
            end
            Q = blkdiag(blocks{:});
        end
    end
end
