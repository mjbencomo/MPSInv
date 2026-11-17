% Quick script to visualize fractional derivatives and integrals.

Nt = 101;
t_grid = linspace(0,1,Nt);
dt = t_grid(2) - t_grid(1);

tc = 0.3;
fpeak = 5;
w = RickerWave(t_grid,tc,fpeak);

figure
plot(t_grid,w)
title("Ricker Wavelet")

%% Fractional derivatives

Dw_1 = fractionalDerivative(w,dt,1);
Dw_half = fractionalDerivative(w,dt,1/2);
Dw_2 = fractionalDerivative(w,dt,2);
Dw_3half = fractionalDerivative(w,dt,3/2);

figure
plot(t_grid,Dw_half)
title("Derivative Order 1/2")

figure
plot(t_grid,Dw_1)
title("Derivative Order 1")

figure
plot(t_grid,Dw_3half)
title("Derivative Order 3/2")

figure
plot(t_grid,Dw_2)
title("Derivative Order 2")

%% Fractional integrals

Iw_1 = fractionalIntegral(w,dt,1);
Iw_half = fractionalIntegral(w,dt,1/2);
Iw_2 = fractionalIntegral(w,dt,2);
Iw_3half = fractionalIntegral(w,dt,3/2);

figure
plot(t_grid,Iw_half)
title("Integral Order 1/2")

figure
plot(t_grid,Iw_1)
title("Integral Order 1")

figure
plot(t_grid,Iw_3half)
title("Integral Order 3/2")

figure
plot(t_grid,Iw_2)
title("Integral Order 2")
