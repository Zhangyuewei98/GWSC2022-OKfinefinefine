function sigVec = crcbgensinsig(dataX,snr,sinCoefs)
% Generate a sinusoidal signal
% S = CRCBGENSINSIG(X,SNR,C)
% Generates a sinusoidal signal S. X is the vector of
% time stamps at which the samples of the signal are to be computed. SNR is
% the matched filtering signal-to-noise ratio of S and C is the vector of
% two coefficients [f0, phi0] that parametrize the phase of the signal:
% 2*pi*f0*t+phi0. 

%Yuewei Zhang, February 2022

phaseVec = 2*pi*sinCoefs(1)*dataX + sinCoefs(2);
sigVec = sin(phaseVec);
sigVec = snr*sigVec/norm(sigVec);


