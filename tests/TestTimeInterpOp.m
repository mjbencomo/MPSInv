classdef TestTimeInterpOp < matlab.unittest.TestCase
    methods (Test)
        function interpolatesNonuniformTimesAndPreservesChannelOrder(tc)
            t = [0;1;3];
            q = [0;0.25;1;2;3];
            X = [0,10;2,20;6,40];
            P = TimeInterpOp(t,q,NumChannels=2);
            expected = [0,10;0.5,12.5;2,20;4,30;6,40];
            tc.verifyTrue(isa(P,'LinearOp'));
            tc.verifyEqual(P.dimDomain,6);
            tc.verifyEqual(P.dimRange,10);
            tc.verifyEqual(P.fwd(X(:)),expected(:),'AbsTol',1e-14);
        end

        function agreesWithInterp1ForComplexData(tc)
            t = [-1;0;0.3;2];
            q = [-1;0.1;0.8;2];
            X = [1,2;3,4;-2,5;0,-1]+1i*[0,2;1,0;3,1;-1,2];
            P = TimeInterpOp(t,q,NumChannels=2);
            expected = interp1(t,X,q,'linear');
            tc.verifyEqual(P.fwd(X(:)),expected(:),'AbsTol',1e-13);
        end

        function supportsUnorderedRepeatedQueries(tc)
            P = TimeInterpOp([0,1,3],[3,0.5,0.5,0]);
            tc.verifyEqual(P.fwd([0;2;6]),[6;1;1;0],'AbsTol',1e-14);
        end

        function identicalGridIsIdentity(tc)
            P = TimeInterpOp([0;0.2;1],[0;0.2;1],NumChannels=2);
            tc.verifyEqual(full(P.matrixRep()),eye(6));
        end

        function adjointAccumulatesWeights(tc)
            P = TimeInterpOp([0;1;2],[0.25;0.75;1.5]);
            tc.verifyEqual(P.adj([2;4;6]),[2.5;6.5;3],'AbsTol',1e-14);
        end

        function matrixAndComplexAdjointAgree(tc)
            P = TimeInterpOp([0;0.2;1],[0;0.1;0.7;1],NumChannels=3);
            x = (1:P.dimDomain).'+1i*(P.dimDomain:-1:1).';
            y = sin((1:P.dimRange).')+1i*cos((1:P.dimRange).');
            A = P.matrixRep();
            tc.verifyTrue(issparse(A));
            tc.verifySize(A,[12,9]);
            tc.verifyEqual(P.fwd(x),A*x,'AbsTol',1e-13);
            tc.verifyEqual(P.adj(y),A'*y,'AbsTol',1e-13);
            tc.verifyEqual(dot(P.fwd(x),y),dot(x,P.adj(y)), ...
                'AbsTol',1e-12);
        end

        function rejectsOutsideTimesByDefault(tc)
            tc.verifyError(@() TimeInterpOp([0;1],[-eps;0.5]), ...
                'TimeInterpOp:OutsideDomain');
            tc.verifyError(@() TimeInterpOp([0;1],[0.5;1+eps]), ...
                'TimeInterpOp:OutsideDomain');
        end

        function zeroExtensionHasZeroAdjointContributions(tc)
            P = TimeInterpOp([0;1],[-1;0;0.5;1;2],OutsideDomain="zero");
            tc.verifyEqual(P.fwd([2;4]),[0;2;3;4;0]);
            tc.verifyEqual(P.adj([100;0;0;0;200]),[0;0]);
            Q = TimeInterpOp([0;1],[-2;2],OutsideDomain="zero");
            tc.verifyEqual(Q.fwd([2;4]),[0;0]);
        end

        function singletonInputAndOutput(tc)
            P = TimeInterpOp(2,2,NumChannels=2);
            tc.verifyEqual(P.fwd([3;4]),[3;4]);
            Q = TimeInterpOp(2,[1;2;3;2],OutsideDomain="zero");
            tc.verifyEqual(Q.fwd(7),[0;7;0;7]);
            tc.verifyEqual(Q.adj([1;2;3;4]),6);
            R = TimeInterpOp([0;2],1);
            tc.verifyEqual(R.fwd([2;4]),3);
            tc.verifyError(@() TimeInterpOp(2,3), ...
                'TimeInterpOp:OutsideDomain');
        end

        function scalarActionsReturnFullVectors(tc)
            % Sparse times scalar can produce sparse results in MATLAB.
            % Exercise the forward and adjoint scalar cases explicitly.
            P = TimeInterpOp(2,[1;2;3;2],OutsideDomain="zero");
            tc.verifyFalse(issparse(P.fwd(7)));
            Q = TimeInterpOp([0;2],1);
            tc.verifyEqual(Q.adj(6),[3;3]);
            tc.verifyFalse(issparse(Q.adj(6)));
            R = TimeInterpOp(2,2,NumChannels=2);
            tc.verifyFalse(issparse(R.fwd([3;4])));
            tc.verifyFalse(issparse(R.adj([3;4])));
            tc.verifyTrue(issparse(R.weights));
            tc.verifyTrue(issparse(R.matrixRep()));
        end

        function rejectsInvalidGridsAndVectorSizes(tc)
            tc.verifyError(@() TimeInterpOp([0;0;1],0.5), ...
                'TimeInterpOp:InvalidInputTimes');
            tc.verifyError(@() TimeInterpOp([1;0],0.5), ...
                'TimeInterpOp:InvalidInputTimes');
            tc.verifyError(@() TimeInterpOp(ones(2),0.5), ...
                'TimeInterpOp:InvalidTimeVector');
            tc.verifyError(@() TimeInterpOp([0;1],ones(2)), ...
                'TimeInterpOp:InvalidTimeVector');
            P = TimeInterpOp([0;1],[0;0.5;1]);
            tc.verifyError(@() P.fwd([1,2]),'LinearOp:InvalidDomainVector');
            tc.verifyError(@() P.fwd(1),'LinearOp:InvalidDomainVector');
            tc.verifyError(@() P.adj([1,2,3]),'LinearOp:InvalidRangeVector');
            tc.verifyError(@() P.adj(1),'LinearOp:InvalidRangeVector');
        end

        function composesBeforeAndAfterConvolution(tc)
            K = 4;
            t = (0:K-1).';
            G = ConvMultiOp(reshape(sin(1:42),7,3,2),Scale=0.2);
            Pin = TimeInterpOp([0;1.5;3],t,NumChannels=G.N_input);
            Pout = TimeInterpOp(t,[0;0.5;2;3],NumChannels=G.N_output);
            x = cos((1:Pin.dimDomain).');
            y = sin((1:Pout.dimRange).');
            fwd = @(v) Pout.fwd(G.fwd(Pin.fwd(v)));
            adj = @(v) Pin.adj(G.adj(Pout.adj(v)));
            A = Pout.matrixRep()*G.matrixRep()*Pin.matrixRep();
            tc.verifyEqual(fwd(x),A*x,'AbsTol',1e-12);
            tc.verifyEqual(adj(y),A'*y,'AbsTol',1e-12);
            tc.verifyEqual(dot(fwd(x),y),dot(x,adj(y)),'AbsTol',1e-12);
        end
    end
end
