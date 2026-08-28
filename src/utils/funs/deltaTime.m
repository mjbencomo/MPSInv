function vals = deltaTime(t,tc,dt)
%DELTATIME Hat-function approximation of the Dirac delta.
%   vals = deltaTime(t, tc, dt) evaluates a piecewise-linear hat function
%   centered at tc with support [tc-dt, tc+dt]. The function has unit area.
%
%   The maximum value is 1/dt at t = tc.
arguments
    t (:,1) double {mustBeVector}
    tc (1,1) double {mustBeReal}
    dt (1,1) double {mustBePositive}
end

vals = max(1 - abs(t - tc)./dt, 0)./dt;

end