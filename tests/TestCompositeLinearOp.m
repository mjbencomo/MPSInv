classdef TestCompositeLinearOp < matlab.unittest.TestCase
    methods (Test)
        function rectangularCompositionHasCorrectOrderAndDimensions(tc)
            A = [1,2,-1;0,3,4];
            B = [2,0,1,3;1,-1,2,0;0,4,1,-2];
            D = [1,2;0,1;3,-1;2,0];
            C = CompositeLinearOp({MatrixLinearOp(A), ...
                MatrixLinearOp(B),MatrixLinearOp(D)});
            x = [2;-3];
            y = [1;4];
            tc.verifyTrue(isa(C,'LinearOp'));
            tc.verifyEqual(C.dimDomain,2);
            tc.verifyEqual(C.dimRange,2);
            tc.verifyEqual(C.numOperators,3);
            tc.verifyEqual(C.fwd(x),A*(B*(D*x)));
            tc.verifyEqual(C.adj(y),D'*(B'*(A'*y)));
            tc.verifyEqual(C.matrixRep(),A*B*D);
        end

        function complexAdjointAndMatrixAgree(tc)
            A = [1+2i,3,-1i;2,1i,4];
            B = [1,2i;3-1i,4;2,1+1i];
            C = CompositeLinearOp({MatrixLinearOp(A),MatrixLinearOp(B)});
            x = [2+1i;-3+2i];
            y = [1-1i;4+2i];
            M = C.matrixRep();
            tc.verifyEqual(M,A*B,'AbsTol',1e-12);
            tc.verifyEqual(C.fwd(x),M*x,'AbsTol',1e-12);
            tc.verifyEqual(C.adj(y),M'*y,'AbsTol',1e-12);
            tc.verifyEqual(dot(C.fwd(x),y),dot(x,C.adj(y)), ...
                'AbsTol',1e-12);
        end

        function singleOperatorIsSupported(tc)
            A = MatrixLinearOp([1,2,3;4,5,6]);
            C = CompositeLinearOp({A});
            tc.verifyEqual(C.dimDomain,3);
            tc.verifyEqual(C.dimRange,2);
            tc.verifyEqual(C.fwd([1;2;3]),A.fwd([1;2;3]));
            tc.verifyEqual(C.adj([2;3]),A.adj([2;3]));
            tc.verifyEqual(C.matrixRep(),A.matrixRep());
        end

        function columnCellVectorIsNormalized(tc)
            C = CompositeLinearOp({MatrixLinearOp([1,2]); ...
                MatrixLinearOp([3;4])});
            tc.verifySize(C.operators,[1,2]);
            tc.verifyEqual(C.fwd(2),22);
        end

        function homogeneousObjectArrayIsSupported(tc)
            A = MatrixLinearOp([1,2]);
            B = MatrixLinearOp([3;4]);
            C = CompositeLinearOp([A,B]);
            tc.verifyEqual(C.matrixRep(),11);
        end

        function nestedCompositionIsSupported(tc)
            A = MatrixLinearOp([1,2;3,4]);
            B = MatrixLinearOp([2,0;1,3]);
            D = MatrixLinearOp([1;2]);
            C = CompositeLinearOp({CompositeLinearOp({A,B}),D});
            tc.verifyEqual(C.matrixRep(),A.matrixRep()*B.matrixRep()*D.matrixRep());
            tc.verifyEqual(C.fwd(2),A.fwd(B.fwd(D.fwd(2))));
            tc.verifyEqual(C.adj([1;2]),D.adj(B.adj(A.adj([1;2]))));
        end

        function rejectsIncompatibleOperators(tc)
            tc.verifyError(@() CompositeLinearOp({ ...
                MatrixLinearOp(ones(2,3)),MatrixLinearOp(ones(4,5))}), ...
                'CompositeLinearOp:IncompatibleDimensions');
            % Ensure later adjacent pairs are checked too.
            tc.verifyError(@() CompositeLinearOp({ ...
                MatrixLinearOp(ones(2,3)),MatrixLinearOp(ones(3,4)), ...
                MatrixLinearOp(ones(5,2))}), ...
                'CompositeLinearOp:IncompatibleDimensions');
        end

        function rejectsInvalidContainersAndEntries(tc)
            A = MatrixLinearOp(1);
            tc.verifyError(@() CompositeLinearOp({}), ...
                'CompositeLinearOp:InvalidOperators');
            tc.verifyError(@() CompositeLinearOp([1,2]), ...
                'CompositeLinearOp:InvalidOperators');
            tc.verifyError(@() CompositeLinearOp({A,A;A,A}), ...
                'CompositeLinearOp:InvalidOperators');
            tc.verifyError(@() CompositeLinearOp({A,3}), ...
                'CompositeLinearOp:InvalidOperator');
            tc.verifyError(@() CompositeLinearOp({A,{A}}), ...
                'CompositeLinearOp:InvalidOperator');
            tc.verifyError(@() CompositeLinearOp({[A,A]}), ...
                'CompositeLinearOp:InvalidOperator');
        end

        function validatesExternalVectors(tc)
            C = CompositeLinearOp({MatrixLinearOp(ones(2,3))});
            tc.verifyError(@() C.fwd(zeros(1,3)), ...
                'LinearOp:InvalidDomainVector');
            tc.verifyError(@() C.fwd(zeros(2,1)), ...
                'LinearOp:InvalidDomainVector');
            tc.verifyError(@() C.adj(zeros(1,2)), ...
                'LinearOp:InvalidRangeVector');
            tc.verifyError(@() C.adj(zeros(3,1)), ...
                'LinearOp:InvalidRangeVector');
        end

        function interpolationAndConvolutionCompose(tc)
            t = (0:3).';
            G = ConvMultiOp(reshape(sin(1:42),7,3,2),Scale=0.2);
            Pin = TimeInterpOp([0;1.5;3],t,NumChannels=G.N_input);
            Pout = TimeInterpOp(t,[0;0.5;2],NumChannels=G.N_output);
            C = CompositeLinearOp({Pout,G,Pin});
            x = cos((1:Pin.dimDomain).');
            y = sin((1:Pout.dimRange).');
            M = Pout.matrixRep()*G.matrixRep()*Pin.matrixRep();
            tc.verifyEqual(C.dimDomain,6);
            tc.verifyEqual(C.dimRange,9);
            tc.verifyEqual(C.fwd(x),M*x,'AbsTol',1e-12);
            tc.verifyEqual(C.adj(y),M'*y,'AbsTol',1e-12);
            tc.verifyEqual(C.matrixRep(),M,'AbsTol',1e-12);
            tc.verifyEqual(dot(C.fwd(x),y),dot(x,C.adj(y)), ...
                'AbsTol',1e-12);
        end
    end
end
