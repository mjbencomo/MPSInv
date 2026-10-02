classdef TestCorrMulti < matlab.unittest.TestCase
    methods (Test)
        function isAdjointOfMultiConvolution(testCase)
            x = rand(8,3);
            y = rand(8,5);
            g = rand(15,5,3);

            Lx = convMulti(g,x);
            Ladjy = corrMulti(g,y);

            dot_x = sum(dot(Lx,y))
            dot_y = sum(dot(x,Ladjy))

            testCase.verifyEqual(dot_x,dot_y,'RelTol',10*eps);
        end
    end
end
