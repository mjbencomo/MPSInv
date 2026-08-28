function vals = RickerWave(t,tc,fpeak)
% RICKERWAVE Ricker wavelet function definition from seismology.
% 
% vals = RickerWave(t,tc,fpeak)
%
% Inputs:
%     t - time points
%    tc - time center of wavelet
% fpeak - peak frequency of wavelet

arguments
    t (:,1) double {mustBeVector}
    tc (1,1) double {mustBeReal}
    fpeak (1,1) double {mustBePositive}
end
t = t(:);
tmp = ((t-tc)*pi*fpeak).^2;
vals = (1-2*tmp).* exp( -tmp );

end