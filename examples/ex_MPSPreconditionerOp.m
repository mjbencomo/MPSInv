% Example of initializing MPSPreconditionerOp

c = 2;
d = 2;
s = [0,0; 1,0; 0,1]; % Monopole, x dipole, y dipole.
Nt = 301;
dt = 0.01;
t = (0:Nt-1).'*dt;

B = MPSPreconditionerOp(c,d,s,Nt,dt);
W = [sin(pi*t).^2, t.^2.*(1-t).^2, sin(2*pi*t).^2];
w = W(:); % Complete time series for each channel, stacked in order.

Bw = B.fwd(w);
wRecovered = B.inv(Bw);
fprintf('Relative recovery error: %.3e\n',norm(wRecovered-w)/norm(w));