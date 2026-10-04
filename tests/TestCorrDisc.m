classdef TestCorrDisc < matlab.unittest.TestCase
    methods (Test)
        function isAdjointOfDiscreteConvolution(testCase)
            x = rand(8,1);
            y = rand(8,1);
            g = rand(15,1);

            Lx = convDisc(g,x);
            Ladjy = corrDisc(g,y);

            dot_x = dot(Lx,y);
            dot_y = dot(x,Ladjy);

            testCase.verifyEqual(dot_x,dot_y,'RelTol',10*eps);
        end
    end
end
