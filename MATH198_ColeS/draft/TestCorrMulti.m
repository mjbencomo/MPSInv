classdef TestCorrMulti < matlab.unittest.TestCase
    methods (Test)
        function isAdjointOfMultiConvolution(testCase)
            x = rand(8,1);
            y = rand(8,1);
            g = rand(15,1);

            Lx = convMulti(g,x);
            Ladjy = corrMulti(g,y);

            dot_x = dot(Lx,y)
            dot_y = dot(x,Ladjy)

            testCase.verifyEqual(dot_x,dot_y,'AbsTol',10*eps);
        end
    end
end
