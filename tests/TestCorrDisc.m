classdef TestCorrDisc < matlab.unittest.TestCase
    % Independent dense references; no FFT or convolution in the oracle.
    % Run with the project's src tree on the MATLAB path.
    properties (TestParameter)
        vectorLength = {1, 2, 5, 8, 17};
        rowOutput = {false, true};
        precision = {'double', 'single'};
    end
    methods (Test)
        function matchesDenseTranspose(t,vectorLength,rowOutput,precision)
            K = vectorLength;
            g = cast(cos((1:2*K-1).') + (1:2*K-1).'/7,precision);
            y = cast(sin((1:K).') - 0.3,precision);
            A = TestCorrDisc.matrix(double(g),K);
            expected = cast(A.'*double(y),precision);
            if rowOutput
                y = y.';
                expected = expected.';
                % corrDisc accepts either kernel orientation.
            end
            actual = corrDisc(g,y);
            t.verifySize(actual,size(y));
            t.verifyClass(actual,precision);
            t.verifyTrue(isreal(actual));
            t.verifyEqual(actual,expected,'AbsTol', ...
                cast(100*K*eps(precision)*max(1,norm(double(expected),inf)),precision));
        end
        function satisfiesAdjointIdentity(t,vectorLength)
            K = vectorLength;
            g = cos((1:2*K-1).');
            x = sin((1:K).');
            y = cos(2*(1:K).');
            Ax = convDisc(g,x);
            ATy = corrDisc(g,y);
            scale = max([1,norm(Ax)*norm(y),norm(x)*norm(ATy)]);
            t.verifyEqual(Ax.'*y,x.'*ATy,'AbsTol',100*K*eps*scale);
        end
        function impulsePinsLagDirection(t)
            y = [2;-3;5;7];
            g = zeros(7,1);
            g(5) = 1; % lag +1: correlation advances y
            t.verifyEqual(corrDisc(g,y),[y(2:end);0],'AbsTol',1e-12);
            g(:) = 0;
            g(3) = 1; % lag -1
            t.verifyEqual(corrDisc(g,y),[0;y(1:end-1)],'AbsTol',1e-12);
        end
        function zeroLagAndZeroInputs(t)
            y = [2;-3;5;7];
            g = zeros(7,1);
            g(4) = 1;
            t.verifyEqual(corrDisc(g,y),y,'AbsTol',1e-12);
            t.verifyEqual(corrDisc(0*g,y),zeros(4,1),'AbsTol',1e-12);
            t.verifyEqual(corrDisc(g,0*y),zeros(4,1),'AbsTol',1e-12);
        end
        function rejectsWrongSize(t)
            t.verifyError(@() corrDisc([1;2],[1;2;3]),'corrDisc:SizeMismatch');
        end
        function acceptsAllOrientationPairs(t)
            g = [4;-2;1;3;5];
            y = [2;-1;3];
            expected = [14;4;13]; % [1 -2 4;3 1 -2;5 3 1].'*y
            for kernelRow = [false,true]
                for dataRow = [false,true]
                    gg = g; yy = y; ee = expected;
                    if kernelRow, gg = gg.'; end
                    if dataRow, yy = yy.'; ee = ee.'; end
                    actual = corrDisc(gg,yy);
                    t.verifySize(actual,size(yy));
                    t.verifyEqual(actual,ee,'AbsTol',1e-12);
                end
            end
        end
        function extremeLagsDoNotWrap(t)
            y = [2;-3;5;7];
            g = zeros(7,1); g(7) = 1;
            t.verifyEqual(corrDisc(g,y),[7;0;0;0],'AbsTol',1e-12);
            g(:) = 0; g(1) = 1;
            t.verifyEqual(corrDisc(g,y),[0;0;0;2],'AbsTol',1e-12);
        end
        function rejectsInvalidArguments(t)
            g = ones(5,1); y = ones(3,1);
            invalid = {[],ones(2,2),[1;2;3]+1i,{1}};
            for k = 1:numel(invalid)
                bad = invalid{k};
                t.verifyTrue(TestCorrDisc.throws(@() corrDisc(bad,y)));
                t.verifyTrue(TestCorrDisc.throws(@() corrDisc(g,bad)));
            end
        end
    end
    methods (Static, Access = private)
        function A = matrix(g,K)
            A = zeros(K);
            for i = 1:K
                for j = 1:K
                    A(i,j) = g(K+i-j);
                end
            end
        end
        function tf = throws(f)
            % Built-in argument-validation IDs can vary across MATLAB releases.
            tf = false;
            try
                f();
            catch
                tf = true;
            end
        end
    end
end
