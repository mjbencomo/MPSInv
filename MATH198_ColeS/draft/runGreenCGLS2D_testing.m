clear; clc;
%% Experiment 1
rcv_loc = [0.2, 0.5; 0.5, 0.5; 0.8, 0.5];
src_location0 = [0.33,0.43];
src_loc = [src_location0; src_location0; src_location0];
derivativeOrders = [0,0; 1,0; 0,1];
tc0 = 0.35; tc1 = 0.35; tc2 = 0.35;
w = {@(t) exp(-((t-tc0)/0.1).^2), ...
    @(t) exp(-((t-tc1)/0.1).^2), ...
    @(t) exp(-((t-tc2)/0.1).^2)};

Experiment1 = runGreenCGLS2D(rcv_loc,src_loc,derivativeOrders,w);
plotGreenCGLS2D(Experiment1);

%% Experiment 2
