clear; clc;
%% Experiment 1
% This experiment is a control of sorts. Standard CGLS with noise. 3
% receivers, 1 multipole source consisting of a monopole and 2 dipoles.

rcv_loc = [0.2, 0.5; 0.5, 0.5; 0.8, 0.5];
src_location0 = [0.5,0.45];
src_loc = [src_location0; src_location0; src_location0];
derivativeOrders = [0,0; 1,0; 0,1];
tc0 = 0.35; tc1 = 0.35; tc2 = 0.35;
w = {@(t) exp(-((t-tc0)/0.1).^2), ...
    @(t) exp(-((t-tc1)/0.1).^2), ...
    @(t) exp(-((t-tc2)/0.1).^2)};

Experiment1 = runGreenCGLS2D(rcv_loc,src_loc,derivativeOrders,w);
plotGreenCGLS2D(Experiment1);

%% Experiment 2a
% Here we will change the grid of receivers into a rectangular 
% perimeter around the source. We also change the source into different
% cases for each type of multipole.

%Creating receiver grid
x_rcv = linspace(0.25,0.75,10)';
y_rcv = linspace(0.25,0.65,8)';

rcv_loc = [
    x_rcv, 0.25*ones(numel(x_rcv),1);
    x_rcv, 0.65*ones(numel(x_rcv),1);
    0.25*ones(numel(y_rcv)-2,1), y_rcv(2:end-1);
    0.75*ones(numel(y_rcv)-2,1), y_rcv(2:end-1)];
src_location0 = [0.5,0.45];
src_loc = [src_location0];
derivativeOrders = [0,0];
tc0 = 0.35;
w = {@(t) exp(-((t-tc0)/0.1).^2)};
kappa = 1;
beta = 1;

Experiment2a = runGreenCGLS2D(rcv_loc,src_loc,derivativeOrders,w, ...
    KappaValue=kappa, BetaValue=beta);
plotGreenCGLS2D(Experiment2a);

%% Experiment 2b
% This is the same experiment as 2a but with dipoles and no monopole.

%Creating receiver grid
x_rcv = linspace(0.25,0.75,10)';
y_rcv = linspace(0.25,0.65,8)';

rcv_loc = [
    x_rcv, 0.25*ones(numel(x_rcv),1);
    x_rcv, 0.65*ones(numel(x_rcv),1);
    0.25*ones(numel(y_rcv)-2,1), y_rcv(2:end-1);
    0.75*ones(numel(y_rcv)-2,1), y_rcv(2:end-1)];
src_location0 = [0.5,0.45];
src_loc = [src_location0; src_location0];
derivativeOrders = [1,0; 0,1];
tc0 = 0.35; tc1 = 0.35;
w = {@(t) exp(-((t-tc0)/0.1).^2), ...
    @(t) exp(-((t-tc1)/0.1).^2),};
kappa = 1;
beta = 1;

Experiment2b = runGreenCGLS2D(rcv_loc,src_loc,derivativeOrders,w, ...
    KappaValue=kappa, BetaValue=beta);
plotGreenCGLS2D(Experiment2b);

%% Experiment 2c
% This is the same experiment as 2a but a multipole consisting of a
% monopole and a dipole.

%Creating receiver grid
x_rcv = linspace(0.25,0.75,10)';
y_rcv = linspace(0.25,0.65,8)';

rcv_loc = [
x_rcv, 0.25*ones(numel(x_rcv),1);
x_rcv, 0.65*ones(numel(x_rcv),1);
0.25*ones(numel(y_rcv)-2,1), y_rcv(2:end-1);
0.75*ones(numel(y_rcv)-2,1), y_rcv(2:end-1)];
src_location0 = [0.5,0.45];
src_loc = [src_location0; src_location0; src_location0];
derivativeOrders = [0,0; 1,0; 0,1];
tc0 = 0.35; tc1 = 0.35;
w = {@(t) exp(-((t-tc0)/0.1).^2), ...
    @(t) exp(-((t-tc1)/0.1).^2), ...
    @(t) exp(-((t-tc2)/0.1).^2)};
kappa = 1;
beta = 1;

Experiment2c = runGreenCGLS2D(rcv_loc,src_loc,derivativeOrders,w, ...
KappaValue=kappa, BetaValue=beta);
plotGreenCGLS2D(Experiment2c);
%% Experiment 3
% Here we will adjust our bulk modulus so the speed of sound through the
% medium is 1.5 km/s. Keep the Buoyancy at 1 cm^3/g. It follows that the
% bulk modulus is 0.00225 Pa. We will be running this on with a monopole
% source to simplify the experiment

%Creating receiver grid
x_rcv = linspace(0.25,0.75,10)';
y_rcv = linspace(0.25,0.65,8)';

rcv_loc = [
    x_rcv, 0.25*ones(numel(x_rcv),1);
    x_rcv, 0.65*ones(numel(x_rcv),1);
    0.25*ones(numel(y_rcv)-2,1), y_rcv(2:end-1);
    0.75*ones(numel(y_rcv)-2,1), y_rcv(2:end-1)];
src_location0 = [0.5,0.45];
src_loc = [src_location0];
derivativeOrders = [0,0];
tc0 = 0.35;
w = {@(t) exp(-((t-tc0)/0.1).^2)};
kappa = 0.00225;
beta = 1;

Experiment3 = runGreenCGLS2D(rcv_loc,src_loc,derivativeOrders,w, ...
    KappaValue=kappa, BetaValue=beta);
plotGreenCGLS2D(Experiment3);