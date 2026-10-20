classdef TestCorrCirc < matlab.unittest.TestCase
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
            g = cast(cos((1:K).') + (1:K).'/7,precision);
            y = cast(sin((1:K).') - 0.3,precision);
            A = TestCorrCirc.matrix(double(g),K);
            expected = cast(A.'*double(y),precision);
            if rowOutput
                y = y.';
                expected = expected.';
                g = g.';
            end
            actual = corrCirc(g,y);
            t.verifySize(actual,size(y));
            t.verifyClass(actual,precision);
            t.verifyTrue(isreal(actual));
            t.verifyEqual(actual,expected,'AbsTol', ...
                cast(100*K*eps(precision)*max(1,norm(double(expected),inf)),precision));
        end
        function satisfiesAdjointIdentity(t,vectorLength)
            K = vectorLength;
            g = cos((1:K).');
            x = sin((1:K).');
            y = cos(2*(1:K).');
            Ax = convCirc(g,x);
            ATy = corrCirc(g,y);
            scale = max([1,norm(Ax)*norm(y),norm(x)*norm(ATy)]);
            t.verifyEqual(Ax.'*y,x.'*ATy,'AbsTol',100*K*eps*scale);
        end
        function impulsePinsLagDirection(t)
            y = [2;-3;5;7];
            g = zeros(4,1);
            g(2) = 1; % lag +1: correlation advances y
            t.verifyEqual(corrCirc(g,y),[y(2:end);y(1)],'AbsTol',1e-12);
            g(:) = 0;
            g(4) = 1; % lag -1
            t.verifyEqual(corrCirc(g,y),[y(end);y(1:end-1)],'AbsTol',1e-12);
        end
        function zeroLagAndZeroInputs(t)
            y = [2;-3;5;7];
            g = zeros(4,1);
            g(1) = 1;
            t.verifyEqual(corrCirc(g,y),y,'AbsTol',1e-12);
            t.verifyEqual(corrCirc(0*g,y),zeros(4,1),'AbsTol',1e-12);
            t.verifyEqual(corrCirc(g,0*y),zeros(4,1),'AbsTol',1e-12);
        end
        function rejectsWrongSize(t)
            t.verifyError(@() corrCirc([1;2],[1;2;3]),'corrCirc:SizeMismatch');
        end
        function rejectsMixedOrientations(t)
            t.verifyError(@() corrCirc([1,2,3],[1;2;3]),'corrCirc:SizeMismatch');
        end
        function rejectsInvalidArguments(t)
            g = ones(3,1); y = ones(3,1);
            invalid = {[],ones(2,2),[1;2;3]+1i,{1}};
            for k = 1:numel(invalid)
                bad = invalid{k};
                t.verifyTrue(TestCorrCirc.throws(@() corrCirc(bad,y)));
                t.verifyTrue(TestCorrCirc.throws(@() corrCirc(g,bad)));
            end
        end
    end
    methods (Static, Access = private)
        function A = matrix(g,K)
            A = zeros(K);
            for i = 1:K
                for j = 1:K
                    A(i,j) = g(mod(i-j,K)+1);
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
