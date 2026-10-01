classdef TestCorrCirc < matlab.unittest.TestCase
    methods (Test)
        function isAdjointOfCircularConvolution(testCase)
            x = rand(8,1);
            y = rand(8,1);
            g = rand(8,1);

            Lx = convCirc(g,x);
            Ladjy = corrCirc(g,y);

            dot_x = dot(Lx,y)
            dot_y = dot(x,Ladjy)

            testCase.verifyEqual(dot_x,dot_y,'AbsTol',10*eps);
        end
    end
end
